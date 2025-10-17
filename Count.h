#include "../utils/communication.h"
#include "../blogel/BVertex.h"
#include "../blogel/Block.h"
#include "../blogel/BWorker.h"
#include "../blogel/BGlobal.h"
#include "../blogel/BType.h"
#include <queue>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <unordered_map>
#include <string.h>
#include <map>
#include <numeric>
#include <utility>
#include <time.h>
#include <deque>
#include <set>
#include <algorithm>
#include <bitset>
#include <vector>
#include <omp.h>
#include <malloc.h>

// 核心点：并行化树分解，如何确认顶点之间的互不影响性
// 结合分区和树分解的性质，进一步加速分解过程
#define MAXINT ((unsigned) 4294967295)

using namespace std;

struct EnumValue {
	vector<int> Local, Edge; 

    void set_Edge(int adj_id, int id1){

		if (v2part[id1] != v2part[adj_id])
			Edge.push_back(adj_id);
		else
			Local.push_back(adj_id);
    }

	void empty(){
		vector<int>().swap(Local);
		vector<int>().swap(Edge);
	}

    friend ibinstream & operator<<(ibinstream & m, const EnumValue & v){
		m<<v.Local;
		m<<v.Edge;

    	return m;
    }


    friend obinstream & operator>>(obinstream & m, EnumValue & v){
		m>>v.Local;
    	m>>v.Edge;
    	
		return m;
    }
};


struct MsgInf {
};


class LCRVertex : public BVertex<VertexID, EnumValue, MsgInf>{
public:
	virtual void compute(MessageContainer& messages){}

	// =======边界索引构建策略=========、
	void Bound_Compute(MessageContainer& messages){ // 按道理讲，这个过程不需要接收其他分区顶点的信息

		vote_to_halt();
	}

};


class LCRBlock : public Block<char, LCRVertex, MsgInf> {

};


class LCRBlockWorker : public BWorker<LCRBlock>{
public:
	vector<unsigned> *unsigned_lab, *unsigned_clab;
	
	vector<long> *long_lab, *long_clab;
	
	vector<int> *int_pos, *int_cpos, 
	            *long_pos, *long_cpos;

	vector<vector<unsigned> > con, conE;
	vector<vector<long> > conB;

	vector<vector<int> > rootlab;

    int BoundNum = 0, max_id_bit = 0, topk = 8;
	int id_bit = 0, dis_bit = 0, cnt_bit = 0;
	long long labelsize = 0, cntnum = 0;
	long mask;
	virtual LCRVertex *toVertex(char *line){
		LCRVertex *v;
		return v;
	}

	int countBits(unsigned value) {
        if (value == 0) return 1;
        if (value == 1) return 1;
        
        int bits = 0;
        unsigned temp = value;
        while (temp > 0) {
            bits++;
            temp >>= 1;
        }
        return bits;
    }

	unsigned merge32(unsigned id, unsigned cnt) {		
		return (id << (32 - max_id_bit)) | cnt; // 格式: [id: max_id_bit位][cnt: (32-max_id_bit)位]
	}

	void split32(unsigned combined, unsigned& id, unsigned& cnt) {
		unsigned id_mask = (1U << max_id_bit) - 1;
		unsigned cnt_mask = (1U << (32 - max_id_bit)) - 1;
		
		id = (combined >> (32 - max_id_bit)) & id_mask;  // 提取高max_id_bit位作为id
		cnt = combined & cnt_mask;                       // 提取低(32-max_id_bit)位作为cnt
	}

	long merge64(unsigned id, unsigned cnt) {
		// 确保 id 不超过 max_id_bit 位
		long masked_id = id & ((1ULL << max_id_bit) - 1);
		return (masked_id << (64 - max_id_bit)) | (static_cast<long>(cnt) & ((1ULL << (64 - max_id_bit)) - 1));
	}

	void split64(long combined, unsigned& id, unsigned& cnt) {
		id = static_cast<unsigned>((combined >> (64 - max_id_bit)) & ((1ULL << max_id_bit) - 1));
		cnt = static_cast<unsigned>(combined & ((1ULL << (64 - max_id_bit)) - 1));
	}

	long merge_edge(unsigned id, unsigned dis, unsigned cnt) {
		return ((long)(id & ((1UL << id_bit) - 1)) << (dis_bit + cnt_bit)) |
			((long)(dis & ((1UL << dis_bit) - 1)) << cnt_bit) |
			(cnt & ((1UL << cnt_bit) - 1));
	}

	// 分解函数：将long变量elem分解成id, dis, cnt三个unsigned变量
	void split_edge(long elem, unsigned& id, unsigned& dis, unsigned& cnt) {
		cnt = elem & ((1UL << cnt_bit) - 1);
		dis = (elem >> cnt_bit) & ((1UL << dis_bit) - 1);
		id = (elem >> (dis_bit + cnt_bit)) & ((1UL << id_bit) - 1);
	}

	void ParameterInitial(){
		unsigned aa = countBits(BoundNum);
		max_id_bit = all_max(aa);

		v2p.resize(totalV, -1);
		p2v.resize(totalV, -1);
		v2p_Global.resize(totalV, -1);
		p2v_Global.resize(totalV, -1);	
	}
	
    virtual void blockInit(VertexContainer &vertexes, BlockContainer &blocks){

		ParameterInitial();

		for (int i=0; i<v2degree.size(); ++i){
			int flg = v2degree[i] > 0? 1 : 0;
			int valB = flg * 1000000000 + v2degree_R[i] - v2degree[i]; // 优先copy，其次bound，
			int valC = flg * 1000000000 + v2degree_R[i];

			if (v2part[i] == _my_rank) sortList.push_back(make_pair(valB, -i));

			if (v2degree[i] > 0)       sortTotal.push_back(make_pair(valC, -i)); // pos就是id					
		}
		
		sort(sortList.rbegin(), sortList.rend());
		
		for (int i=0; i<sortList.size(); ++i){
			v2p[-sortList[i].second] = i;
			p2v[i] = -sortList[i].second; // new 2 old
		}
		
		sort(sortTotal.rbegin(), sortTotal.rend());
		
		for (int i=0; i<sortTotal.size(); ++i){
			v2p_Global[-sortTotal[i].second] = i;
			p2v_Global[i] = -sortTotal[i].second; // new 2 old
		}

		// ======================================
        con.resize(sortList.size());
		conB.resize(sortTotal.size());
		conE.resize(sortTotal.size());
		

		for (int i=0; i<sortList.size(); ++i){
            int ovid = p2v[i];
            vector<int>& local = vertexes[vert2place[ovid]]->value().Local;
            
			for( int p = 0; p < local.size(); ++p ) {
                int j = v2p[local[p]];

                con[j].push_back(i); // 内部存储的
		    }
			
			// ===================================
			vector<int>& edge = vertexes[vert2place[ovid]]->value().Edge;
			
			if (edge.size() > 0){
				
				int newPos = v2p_Global[ovid]; // 在边界图上的定位
				
				for( int p = 0; p < edge.size(); ++p ) {
					unsigned newVid = v2p_Global[edge[p]]; // 在边界图上的定位
					conE[newPos].push_back(newVid);
				}
			}

			vertexes[vert2place[ovid]]->value().empty();
        }

		for (auto* vertex : vertexes) 
			delete vertex;  // 删除每个对象
		vector<LCRVertex *>().swap(vertexes);

		vector<pair<int, int>>().swap(sortList);
		vector<pair<int, int>>().swap(sortTotal);	

		for (int i=0; i<conE.size(); i++){
			sort(conE[i].begin(), conE[i].end());
		}

	}


	void AddVertex(char *line, int va){
		// 三个元素  V_A, V_B, Label
		int vb, pa, pb;
		LCRVertex* v = new LCRVertex;
		v->id = va; v->bid = 0;
		load_vertex(v);
		vert2place[va] = vertexes.size()-1;
		char* s1 = strtok(line," ");

		while(s1){
			
			vb = atoi(s1) - 1;
			
			if ((va+vb)%10<10){
				v->value().set_Edge(vb, va);
				totalEdges += 1;	
			}
		
			s1=strtok(NULL," ");
		}
		
		v2degree[va] = v->value().Edge.size();
		v2degree_R[va] = v->value().Edge.size() + v->value().Local.size();

		if (v->value().Edge.size() > 0) BoundNum += 1;
	}

	int can_update(int v, int dis, char *nowdis) {
		
		unsigned flg = 2, d = 0, w, cc;
		
		if (int_pos[v][0] == 1 and nowdis[v] >= 0){
			if (nowdis[v] < dis) return 0;
			else if (nowdis[v] == dis) flg = 1;
		}
		
		while (d+1 < dis){
			for (int i=int_pos[v][d]; i<int_pos[v][d+1]; ++i){
				split32(unsigned_lab[v][i], w, cc);
				if( nowdis[w] >= 0 ){
					if (nowdis[w] + d + 1 < dis) return 0;
					else if (nowdis[w] + d + 1 == dis) flg = 1;
					else continue;
				}
			}

			for (int i=long_pos[v][d]; i<long_pos[v][d+1]; ++i){
				split64(long_lab[v][i], w, cc);
				if( nowdis[w] >= 0 ){
					if (nowdis[w] + d + 1 < dis) return 0;
					else if (nowdis[w] + d + 1 == dis) flg = 1;
					else continue;
				}
			}

			d += 1;
		}
		
		return flg;
	}


	int prune_by_root(int v, int u, int aaa){
		for (int i=0; i<topk; ++i){
			int dd = rootlab[v][i] + rootlab[u][i];

			if (dd < aaa) {
				// cout<<i<<"  "<<rootlab[v][i]<<"  "<<rootlab[u][i]<<endl;
				return 1;
			}
		}

		return 0;
	}


	void LabelInitial(int Num){
		unsigned_lab  = new vector<unsigned>[Num]; // id+cnt
		int_pos = new vector<int>[Num]; // 对应32位的label

		unsigned_clab = new vector<unsigned>[Num];
		int_cpos = new vector<int>[Num]; // 对应32位的clab

		long_lab  = new vector<long>[Num]; // id+cnt 64位
		long_pos = new vector<int>[Num]; // 对应64位的label

		long_clab = new vector<long>[Num];
		long_cpos = new vector<int>[Num]; // 对应64位的clab
	}

	void ValueUpdate(char *nowdis, unsigned w, int dis){

		unsigned v, cc, d = 0;

		if (int_pos[w][0] == 1){
			nowdis[w] = 0; // 这个label一定是32位的
		} 

		while (d+1 < dis){
			for (int i=int_pos[w][d]; i<int_pos[w][d+1]; ++i){
				split32(unsigned_lab[w][i], v, cc);
				nowdis[v] = d+1;
			} 
			
			for (int i=long_pos[w][d]; i<long_pos[w][d+1]; ++i){
				split64(long_lab[w][i], v, cc);
				nowdis[v] = d+1;
			}

			d += 1;
		}
	}

	void ValueReinitial(char *nowdis, unsigned w, int dis){

		unsigned v, cc, d = 0;

		if (int_pos[w][0] == 1) nowdis[w] = -1;

		while (d+1 < dis){
			for (int i=int_pos[w][d]; i<int_pos[w][d+1]; ++i){
				split32(unsigned_lab[w][i], v, cc);
				nowdis[v] = -1;
			}
			
			for (int i=long_pos[w][d]; i<long_pos[w][d+1]; ++i){
				split64(long_lab[w][i], v, cc);
				nowdis[v] = -1;
			}

			d += 1;
		}
	}

	void LabelUpdate(vector<unsigned>& L1, vector<int>& Pos1, vector<unsigned>& newL){
		L1.insert(L1.end(), newL.begin(), newL.end());
		vector<unsigned>(L1).swap(L1);
		vector<unsigned>().swap(newL);
		Pos1.push_back(L1.size());
	}

	void LabelUpdate(vector<long>& L1, vector<int>& Pos1, vector<long>& newL){
		L1.insert(L1.end(), newL.begin(), newL.end());
		vector<long>(L1).swap(L1);
		vector<long>().swap(newL);
		Pos1.push_back(L1.size());
	}

	void RedgeCount() {
		
		long long cnt = 0;
		int maxdis = int_pos[0].size()-1;
		char *nowdis = new char[BdV];
		memset( nowdis, -1, sizeof(char) * BdV);
		
		for (int u=0; u<BdV; ++u){
			
			ValueUpdate(nowdis, u, maxdis);

			for (int i=0; i<conB[u].size(); i++){
				unsigned v, d, c;
				split_edge(conB[u][i], v, d, c);

				if (d == 1 and v <= u) continue;

				int ans = can_update(v, d, nowdis);

				if (ans == 2) cnt += 1;
			}

			ValueReinitial(nowdis, u, maxdis);
		}
		
		cout<<"Redundant edges: "<<cnt<<endl;
		delete[] nowdis;
	}

	void Part2hop(){

		totalV = con.size();
		LabelInitial(totalV);

		omp_set_num_threads(threads);

		for(int i = 0; i < totalV; ++i){
			if ( i < BoundNum ){ 
				unsigned lab = merge32(i, 1); 
				unsigned_lab[i].push_back(lab);
				int_pos[i].push_back(1);
			}else{
				int_pos[i].push_back(0);
			}

			for (int j=0; j<con[i].size(); ++j){
				if (con[i][j] < i and con[i][j]<BoundNum){
					unsigned lab = merge32(con[i][j], 1);
					unsigned_lab[i].push_back(lab);
				}
			}
			
			int_pos[i].push_back(unsigned_lab[i].size());

			// === unsigned_clab，long_lab，long_clab 都是没有收录label的 ===
			int_cpos[i].push_back(0),  int_cpos[i].push_back(0);
			long_pos[i].push_back(0),  long_pos[i].push_back(0);
			long_cpos[i].push_back(0), long_cpos[i].push_back(0);
		}

		int dis = 2;

		for( long long cnt = 1; cnt > 0; ++dis ){
			
			cnt = 0;
			vector<unsigned> *label_un_new = new vector<unsigned>[totalV];	// 32位	 	
			vector<unsigned> *clab_un_new  = new vector<unsigned>[totalV];

			vector<long> *label_long_new = new vector<long>[totalV]; // 64位			
			vector<long> *clab_long_new  = new vector<long>[totalV];

			#pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				long long local_cnt = 0;

				vector<unsigned> cand, candp;
				vector<int> cand2count(totalV, -1); // 用于存储候选点和计数
				
				char *nowdis = new char[totalV];
				memset( nowdis, -1, sizeof(char) * totalV);
			
				for( int u = pid; u < totalV; u += np ){

					cand.clear(), candp.clear();

					for( int i = 0; i < con[u].size(); ++i ){

						unsigned w = con[u][i], v, cc;
						
						for( int j = int_pos[w][dis-2]; j < int_pos[w][dis-1]; ++j ) {
							// unsigned_lab and int_pos 32 位分解
							split32(unsigned_lab[w][j], v, cc);

							if( v >= u ) break;
							
							if (w < BoundNum) cc = 0; // 经过其他边界点时，意味着这部分路径可以被忽视

							if ( cand2count[v] == -1 ){
								cand.push_back(v);
								cand2count[v] = cc; // 初始化计数
							}else 
								cand2count[v] += cc; // 初始化计数 
						}

						for( int j = long_pos[w][dis-2]; j < long_pos[w][dis-1]; ++j ) {
							// unsigned_lab and int_pos 32 位分解
							
							split64(long_lab[w][j], v, cc);

							if( v >= u ) break;

							// 经过其他边界点时，意味着这部分路径的数量可以被忽视
							if (w < BoundNum) cc = 0; 

							if ( cand2count[v] == -1 ){
								cand.push_back(v);
								cand2count[v] = cc; // 初始化计数
							}else 
								cand2count[v] += cc; // 初始化计数 
						}
						
						// =====================
						if (w >= BoundNum){

							for( int j = int_cpos[w][dis-2]; j < int_cpos[w][dis-1]; ++j ) {
                        
								split32(unsigned_clab[w][j], v, cc);

								if( v >= u ) break;
								
								if (cc == 0) continue;

								if ( cand2count[v] == -1 ){
									cand.push_back(v);
									cand2count[v] = cc; // 初始化计数
								}else
									cand2count[v] += cc; // 初始化计数 
							}

							for( int j = long_cpos[w][dis-2]; j < long_cpos[w][dis-1]; ++j ) {
							
								split64(long_clab[w][j], v, cc);

								if( v >= u ) break;

								if (cc == 0) continue;
								
								if ( cand2count[v] == -1 ){
									cand.push_back(v);
									cand2count[v] = cc; // 初始化计数
								}else
									cand2count[v] += cc; // 初始化计数 
							}
						}						
					}

					// === 更新 nowdis 的值 ===
					if (cand.size() == 0) continue;

					ValueUpdate(nowdis, u, dis);
					
					int n_cand = 0;
					for( int i = 0; i < (int) cand.size(); ++i ){

						int ans = can_update(cand[i], dis, nowdis);

						if (ans == 2){
							cand[n_cand++] = cand[i]; // 有路径，需要考虑distance判断
						}else if (ans == 1){
							candp.push_back(cand[i]); // 有路径，但是不用考虑distance判断
						}else{
							cand2count[cand[i]] = -1; // 重置计数
						}
					}

					cand.resize(n_cand);
					sort(cand.begin(), cand.end());

					for( int i = 0; i < (int) cand.size(); ++i ) {
						// == 判断新的label是否大于32位 ==
						int bits = countBits(cand2count[cand[i]]);

						if (bits + max_id_bit <= 32){
							unsigned lab = merge32(cand[i], cand2count[cand[i]]);
							label_un_new[u].push_back(lab);
						}else{
							long lab = merge64(cand[i], cand2count[cand[i]]);
							label_long_new[u].push_back(lab);
 						}

						++local_cnt;
						cand2count[cand[i]] = -1; // 重置计数
					}

					sort(candp.begin(), candp.end());
					for( int i = 0; i < (int) candp.size(); ++i ){
						
						if (cand2count[candp[i]] == 0){
							cand2count[candp[i]] = -1;
							continue;
						}

						int bits = countBits(cand2count[candp[i]]);

						if (bits + max_id_bit <= 32){
							unsigned lab = merge32(candp[i], cand2count[candp[i]]);
							clab_un_new[u].push_back(lab); // 
						}else{
							long lab = merge64(candp[i], cand2count[candp[i]]);
							clab_long_new[u].push_back(lab);
						}

						local_cnt += 1;
						cand2count[candp[i]] = -1; // 重置计数
					}

					// ======== nowdis reinitialize ========
					ValueReinitial(nowdis, u, dis);					
				}
				
				#pragma omp critical
				{
					cnt += local_cnt;
				}
				
				delete[] nowdis;
			}


			#	pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				for( int u = pid; u < totalV; u += np ){

					LabelUpdate(unsigned_lab[u], int_pos[u], label_un_new[u]);
					
					LabelUpdate(unsigned_clab[u], int_cpos[u], clab_un_new[u]);

					LabelUpdate(long_lab[u], long_pos[u], label_long_new[u]);

					LabelUpdate(long_clab[u], long_cpos[u], clab_long_new[u]);
				}
			}

			delete[] label_un_new;
			delete[] clab_un_new;
			delete[] label_long_new;
			delete[] clab_long_new;
		}
	}

	void Core2hop(){

		rootlab.resize(conB.size(), vector<int>(topk,99));

		totalV = conB.size();
		
		max_id_bit = countBits(totalV);

		LabelInitial(totalV);
		
		omp_set_num_threads(threads);

		for(int i = 0; i < totalV; ++i){
			if ( i < totalV ){ 
				unsigned lab = merge32(i, 1); 
				unsigned_lab[i].push_back(lab);
				int_pos[i].push_back(1);
			}else{
				int_pos[i].push_back(0);
			}

			int_cpos[i].push_back(0);
			long_pos[i].push_back(0);
			long_cpos[i].push_back(0);

			if ( i < topk ) rootlab[i][i] = 0;	
		}

		int dis = 1;
		
		for( long long cnt = 1, cntt = 0, TotalCnt = 0; cnt+cntt > 0; ++dis ){
			
			cnt = 0, cntt = 0, TotalCnt = 0;

			vector<unsigned> *label_un_new = new vector<unsigned>[totalV];	// 32位	 	
			vector<unsigned> *clab_un_new  = new vector<unsigned>[totalV];

			vector<long> *label_long_new = new vector<long>[totalV]; // 64位			
			vector<long> *clab_long_new  = new vector<long>[totalV];


			#pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				long long local_cnt = 0, aaa = 0, bbb = 0;

				vector<unsigned> cand, candp;
				vector<pair<unsigned, unsigned> > cand_delete;
				vector<int> cand2count(totalV, -1);

				char *nowdis = new char[totalV];
				memset( nowdis, -1, sizeof(char) * totalV);

				for( int u = pid; u < totalV; u += np ){

					cand.clear(), candp.clear(), cand_delete.clear(); 

					unsigned w, d, c;

					for( int i = 0; i < conB[u].size(); ++i ){
						
						split_edge(conB[u][i], w, d, c);

						if (d > dis or c == 0) 
							continue; // c=0意味着内部最短不是全局最短

						if (d == dis and dis > 2) {
							// 这一部分不受排序函数影响
							cand_delete.push_back(make_pair(w, i)); // 确认这条边是不是不需要的
						}
						
						unsigned v, cc;

						for( int j = d==dis? 0:int_pos[w][dis-d-1]; j < int_pos[w][dis-d]; ++j ) {

							split32(unsigned_lab[w][j], v, cc);
							
							if( v >= u ) break;

							if ( cand2count[v] == -1 ){
								cand.push_back(v);
								cand2count[v] = c*cc; // 初始化计数
							}else 
								cand2count[v] += c*cc; // 初始化计数 
						}

						for( int j = d==dis? 0:long_pos[w][dis-d-1]; j < long_pos[w][dis-d]; ++j ) {
							
							split64(long_lab[w][j], v, cc);

							if( v >= u ) break;

							if ( cand2count[v] == -1 ){
								cand.push_back(v);
								cand2count[v] = c*cc; // 初始化计数
							}else 
								cand2count[v] += c*cc; // 初始化计数 
						}
					
						for( int j = d==dis? 0:int_cpos[w][dis-d-1]; j < int_cpos[w][dis-d]; ++j ) {

							split32(unsigned_clab[w][j], v, cc);
							
							if( v >= u ) break;

							if ( cand2count[v] == -1 ){
								cand.push_back(v);
								cand2count[v] = c*cc; // 初始化计数
							}else 
								cand2count[v] += c*cc; // 初始化计数 
						}
						
						for( int j = d==dis? 0:long_cpos[w][dis-d-1]; j < long_cpos[w][dis-d]; ++j ) {
							
							split64(long_clab[w][j], v, cc);

							if( v >= u ) break;

							if ( cand2count[v] == -1 ){
								cand.push_back(v);
								cand2count[v] = c*cc; // 初始化计数
							}else 
								cand2count[v] += c*cc; // 初始化计数 
						}
					} 


					if (cand.size() == 0 and cand_delete.size() == 0) continue;

					ValueUpdate(nowdis, u, dis);
					
					int n_cand = 0;
					for( int i = 0; i < (int) cand_delete.size(); ++i ){
						unsigned w = cand_delete[i].first;
						int ans = can_update(w, dis, nowdis);
						if (ans == 0){ // cand_delete中的元素要在后续删除
							cand_delete[n_cand++] = cand_delete[i];
						}
					}
					cand_delete.resize(n_cand);


					n_cand = 0;

					for( int i = 0; i < (int) cand.size(); ++i ){

						RepeatNum[u][cand[i]] += 1;

						int aas = prune_by_root(u, cand[i], dis);

						if (aas == 1){
							cand2count[cand[i]] = -1;
							continue;
						} 

						int ans = can_update(cand[i], dis, nowdis);

						if (ans == 2){

							cand[n_cand++] = cand[i]; // 有路径，需要考虑distance判断

							if (cand[i] < topk) {
								rootlab[u][cand[i]] = dis;
							}

						}else if (ans == 1){
							candp.push_back(cand[i]); // 有路径，但是不用考虑distance判断
						}else{
							cand2count[cand[i]] = -1; // 重置计数
						}
					}

					cand.resize(n_cand);
					sort(cand.begin(), cand.end());

									 
					for( int i = 0; i < (int) cand.size(); ++i ) {
						// == 判断新的label是否大于32位 ==
						int bits = countBits(cand2count[cand[i]]);

						if (bits + id_bit <= 32){
							unsigned lab = merge32(cand[i], cand2count[cand[i]]);
							label_un_new[u].push_back(lab);
						}else{
							long lab = merge64(cand[i], cand2count[cand[i]]);
							label_long_new[u].push_back(lab);
 						}

						++local_cnt;
						cand2count[cand[i]] = -1; // 重置计数
					}

					sort(candp.begin(), candp.end());
					for( int i = 0; i < (int) candp.size(); ++i ){

						int bits = countBits(cand2count[candp[i]]);

						if (bits + id_bit <= 32){
							unsigned lab = merge32(candp[i], cand2count[candp[i]]);
							clab_un_new[u].push_back(lab); // 
						}else{
							long lab = merge64(candp[i], cand2count[candp[i]]);
							clab_long_new[u].push_back(lab);
						}

						local_cnt += 1;
						cand2count[candp[i]] = -1; // 重置计数
					}

					for (int i=0; i<cand_delete.size(); i++){
						conB[u][cand_delete[i].second] = conB[u][cand_delete[i].second] & mask;
					}
					// ======== nowdis reinitialize ========
					ValueReinitial(nowdis, u, dis);
				}

				#pragma omp critical
				{
					cnt += local_cnt;
				}
				delete[] nowdis;
			}

			#	pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				for( int u = pid; u < totalV; u += np ){

					LabelUpdate(unsigned_lab[u], int_pos[u], label_un_new[u]);
					
					LabelUpdate(unsigned_clab[u], int_cpos[u], clab_un_new[u]);

					LabelUpdate(long_lab[u], long_pos[u], label_long_new[u]);

					LabelUpdate(long_clab[u], long_cpos[u], clab_long_new[u]);
				}
			}

			// if (_my_rank == 0)
			cout<<"Distance: "<<dis<<"   Cnt: "<<cnt<<endl;


			delete[] label_un_new; 
			delete[] clab_un_new; 
			delete[] label_long_new; 
			delete[] clab_long_new;
		}
	}


	void Global_Init(){ 
		vector<vector<int> > V2D(_num_workers), V2DR(_num_workers);
		for(int i=0; i<_num_workers; ++i){
			if (i == _my_rank) 
				continue;
			V2D[i] = v2degree; // edge.size()
			V2DR[i] = v2degree_R; // edge.size() + local.size()
		}

		all_to_all_cat(V2D, V2DR);
		worker_barrier();

		for (int i=0; i<_num_workers; ++i){
			for (int j=0; j<V2D[i].size(); ++j){
				v2degree[j]   += V2D[i][j];
				v2degree_R[j] += V2DR[i][j]; 
			}
		}

		BdV = all_sum(BoundNum);
		totalV = all_sum(vertexes.size());

		vector<vector<int> >().swap(V2D);
		vector<vector<int> >().swap(V2DR);
	}


	void Lab2Edge(unsigned self_id, vector<unsigned>& lab, vector<int>& pos){
		
		unsigned vid, cc, d = 0;

		while (d < pos.size()-1){
			
			for (int i=pos[d]; i<pos[d+1]; ++i){

				split32(lab[i], vid, cc);

				if (cc == 0) continue; // 该标签只用于距离判断

				int ovid = p2v[vid], ob_vid = v2p_Global[ovid];

				if (self_id == ob_vid) continue;

				long elem1 = merge_edge(ob_vid, d+1, cc), 
				     elem2 = merge_edge(self_id, d+1, cc);

				unsigned vv, dd, cc;
				split_edge(elem1, vv, dd, cc);

				conB[self_id].push_back(elem1);
				conB[ob_vid].push_back(elem2);
			}
			d = d + 1;
		}
	
		vector<unsigned>().swap(lab);
		vector<int>().swap(pos);
	}

	void Lab2Edge(unsigned self_id, vector<long>& lab, vector<int>& pos){
		
		unsigned vid, cc, d = 0;

		while (d < pos.size()-1){
			for (int i=pos[d]; i<pos[d+1]; ++i){

				split64(lab[i], vid, cc);

				if (cc == 0) continue; // 该标签只用于距离判断

				int ovid = p2v[vid], ob_vid = v2p_Global[ovid];

				if (self_id == ob_vid) continue;

				long elem1 = merge_edge(ob_vid, d+1, cc), elem2 = merge_edge(self_id, d+1, cc);

				conB[self_id].push_back(elem1);
				conB[ob_vid].push_back(elem2);
			}

			d = d + 1;
		}

		vector<long>().swap(lab);
		vector<int>().swap(pos);
	}

	void BoundGraph(){

		unsigned maxdis = 64; // all_max(dd);

		conB.resize(BdV);

		dis_bit = countBits(maxdis);
		id_bit = countBits(BdV);
		cnt_bit = 64 - id_bit - dis_bit;
		mask = ~((1UL << cnt_bit) - 1);

		for (int i=0; i<conE.size(); ++i){

			for (int j=0; j<conE[i].size(); ++j){
				unsigned vid = conE[i][j];
				long elem = merge_edge(vid, 1, 1);
				conB[i].push_back(elem);
			}
		}

		vector<vector<unsigned> >().swap(conE);

		for (int i=0; i<BoundNum; ++i){

			int oid = p2v[i], ob_id = v2p_Global[oid];
			
			Lab2Edge(ob_id, unsigned_lab[i], int_pos[i]);
			Lab2Edge(ob_id, unsigned_clab[i], int_cpos[i]);

			Lab2Edge(ob_id, long_lab[i], long_pos[i]);
			Lab2Edge(ob_id, long_clab[i], long_cpos[i]);
		}


		delete[] unsigned_lab; delete[] unsigned_clab;
		delete[] long_lab;     delete[] long_clab;
		delete[] int_pos;      delete[] int_cpos;
		delete[] long_pos;     delete[] long_cpos;

		#ifdef __linux__
			malloc_trim(0);
		#endif
		
		worker_barrier();
		vector<vector<vector<long> > > CB(_num_workers);
		int start = 0, batch = start+batch > BdV? BdV-start:10000;

		while (1){
			CB.resize(_num_workers);
			for (int i=start; i<start+batch; ++i){
				
				if (conB[i].size() == 0) continue;

				conB[i].push_back(i); 
				CB[0].push_back(conB[i]);
				vector<long>().swap(conB[i]);
			}
			
			all_to_all(CB);
			worker_barrier();

			for(int i=0; i<_num_workers; ++i){
				
				for (int j=0; j<CB[i].size(); ++j){
					
					vector<long>& inf = CB[i][j];
					int vid = inf[inf.size()-1];
					inf.pop_back();

					conB[vid].insert(conB[vid].end(), inf.begin(), inf.end());
					
					vector<long>().swap(inf);
				}
			}

			vector<vector<vector<long> > >().swap(CB);

			start += batch;

			int cntt = all_sum(start);

			// if (_my_rank == 0)
			// 	cout<<"cnt: "<<cntt<<endl;
	
			if (cntt == _num_workers*BdV) break;

			if (start+batch > BdV)
				batch = BdV-start;
		}

		if (_my_rank == 0){
			long long edge_cnt = 0;
			for (int i=0; i<conB.size(); ++i){
				sort(conB[i].begin(), conB[i].end());
				edge_cnt += conB[i].size();
			}
			cout<<"Edges: "<<BdV<<"   "<<edge_cnt<<endl;
		}

	}

	void DHI_Store(string new_filename){
		
		ofstream fout((new_filename).c_str(), ios::binary);

		// == 记录 BoundNum，用于后续的计算  ==
		fout.write(reinterpret_cast<const char*>(&BoundNum), sizeof(unsigned)); // 内部起始位置
		
		// == 存储 p2v ==
		int p2v_cnt = p2v.size();
		fout.write(reinterpret_cast<const char*>(&p2v_cnt), sizeof(unsigned)); 
		fout.write(reinterpret_cast<const char*>(p2v.data()), sizeof(unsigned)*p2v_cnt);

		int v2p_cnt = v2p.size();
		fout.write(reinterpret_cast<const char*>(&v2p_cnt), sizeof(unsigned)); 
		fout.write(reinterpret_cast<const char*>(v2p.data()), sizeof(unsigned)*v2p_cnt);
		
		

		// == 先写label[i]的size ==
		
		for (int i=BoundNum; i<totalV; ++i){
			unsigned vid, cc, c = 0;
			
			int maxd = int_pos[i].size();
			vector<int> change_pos(maxd, 0);

			for (int j=0; j<unsigned_lab[i].size(); ++j){
				
				split32(unsigned_lab[i][j], vid, cc);
				
				if (cc > 0)
					unsigned_lab[i][c++] = unsigned_lab[i][j];
				else{
					for (int kk=0; kk<maxd; ++kk){
						if (int_pos[i][kk] > j)
							change_pos[kk] += 1;
					}
				}
			}

			for (int kk=0; kk<maxd; ++kk){
				int_pos[i][kk] -= change_pos[kk];
			}
			unsigned_lab[i].resize(c);

			fout.write(reinterpret_cast<const char*>(&c), sizeof(unsigned)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(unsigned_lab[i].data()), c*sizeof(unsigned));
			fout.write(reinterpret_cast<const char*>(&maxd), sizeof(int)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(int_pos[i].data()), maxd*sizeof(int));
			
			// =======================
			c = unsigned_clab[i].size();

			fout.write(reinterpret_cast<const char*>(&c), sizeof(unsigned)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(unsigned_clab[i].data()), c*sizeof(unsigned));
			fout.write(reinterpret_cast<const char*>(&maxd), sizeof(int)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(int_cpos[i].data()), maxd*sizeof(int));
			
			
			// // =======================
			c = 0;
			vector<int> change_cpos(maxd, 0);

			for (int j=0; j<long_lab[i].size(); ++j){
				split64(long_lab[i][j], vid, cc);
				if (cc > 0) 
					long_lab[i][c++] = long_lab[i][j];
				else
					for (int kk=0; kk<maxd; ++kk){
						if (long_pos[i][kk] > j)
							change_cpos[kk] += 1;
					}
			}

			for (int kk=0; kk<maxd; ++kk){
				long_pos[i][kk] -= change_cpos[kk];
			}
			long_lab[i].resize(c);

			if (long_pos[i][maxd-1] != c) cout<<"error"<<endl;

			fout.write(reinterpret_cast<const char*>(&c), sizeof(long)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(long_lab[i].data()), c*sizeof(long));
			fout.write(reinterpret_cast<const char*>(&maxd), sizeof(int)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(long_pos[i].data()), maxd*sizeof(int));

			// // ====================
			c = long_clab[i].size();
			fout.write(reinterpret_cast<const char*>(&c), sizeof(long)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(long_clab[i].data()), c*sizeof(long));
			fout.write(reinterpret_cast<const char*>(&maxd), sizeof(int)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(long_cpos[i].data()), maxd*sizeof(int));
		}

		fout.close();
	}

	void DHBound_Store(string new_filename){
		
		ofstream fout((new_filename).c_str(), ios::binary);

		// == 记录 BoundNum，用于后续的计算  ==
		fout.write(reinterpret_cast<const char*>(&BdV), sizeof(unsigned)); // 内部起始位置
		
		// == 存储 p2v ==
		int p2v_cnt = p2v_Global.size();
		fout.write(reinterpret_cast<const char*>(&p2v_cnt), sizeof(unsigned)); 
		fout.write(reinterpret_cast<const char*>(p2v_Global.data()), sizeof(unsigned)*p2v_cnt);

		int v2p_cnt = v2p_Global.size();
		fout.write(reinterpret_cast<const char*>(&v2p_cnt), sizeof(unsigned)); 
		fout.write(reinterpret_cast<const char*>(v2p_Global.data()), sizeof(unsigned)*v2p_cnt);

		// == 先写label[i]的size ==
		int maxd = int_pos[0].size();
		
		for (int i=0; i<BdV; ++i){
			
			unsigned c = unsigned_lab[i].size();

			fout.write(reinterpret_cast<const char*>(&c), sizeof(unsigned)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(unsigned_lab[i].data()), c*sizeof(unsigned));
			fout.write(reinterpret_cast<const char*>(&maxd), sizeof(int)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(int_pos[i].data()), maxd*sizeof(int));
			
			// =======================
			c = unsigned_clab[i].size();

			fout.write(reinterpret_cast<const char*>(&c), sizeof(unsigned)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(unsigned_clab[i].data()), c*sizeof(unsigned));
			fout.write(reinterpret_cast<const char*>(&maxd), sizeof(int)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(int_cpos[i].data()), maxd*sizeof(int));
			
			
			// // =======================
			c = long_lab[i].size();

			fout.write(reinterpret_cast<const char*>(&c), sizeof(unsigned)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(long_lab[i].data()), c*sizeof(long));
			fout.write(reinterpret_cast<const char*>(&maxd), sizeof(int)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(long_pos[i].data()), maxd*sizeof(int));

			// // ====================
			c = long_clab[i].size();

			fout.write(reinterpret_cast<const char*>(&c), sizeof(unsigned)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(long_clab[i].data()), c*sizeof(long));
			fout.write(reinterpret_cast<const char*>(&maxd), sizeof(int)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(long_cpos[i].data()), maxd*sizeof(int));
		}

		fout.close();
	}

	void run_LCR(const WorkerParams& params){
		
		ifstream infile, infile1, infile2; // 先把分区文件读取
		
		threads = params.khop;

		const char *part_path = params.partition_path.c_str(); 
		infile.open(part_path);
		string s;
		if(!infile.is_open()){
			cout<<"No such file!"<<endl;
			exit(-1);
		}

		int nnn = 0;
		while(getline(infile, s)){
			char* part = new char[strlen(s.c_str())+1];
			strcpy(part, s.c_str());
			v2part.push_back( atoi(part) ); // 
			v2degree.push_back(0); 
			delete part;
			nnn += 1;
		}
		v2degree_R.resize(v2degree.size());
		
		const char *filepath = params.input_path.c_str(); //"data//Amazon_New.txt";

		// ======================
		infile1.open(filepath);
		if(!infile1.is_open()){
			cout<<"No such file!"<<endl;
			exit(-1);
		}

		int xxx = 0;
		while(getline(infile1, s)){
			if (xxx > 0 and v2part[xxx-1] == _my_rank){
				char* strc = new char[strlen(s.c_str())+1];
				strcpy(strc, s.c_str());
				AddVertex(strc, xxx-1);
				delete strc;
			}
			xxx += 1;
		}

		Global_Init();

		blockInit(vertexes, blocks);

		float t = omp_get_wtime();

		Part2hop();

		// string new_filename = params.output_path + "In_"+ to_String(_my_rank);
		// DHI_Store(new_filename);

		for (int i=BoundNum; i<totalV; ++i){
			unsigned vid, cc;
			labelsize += unsigned_clab[i].size()*4;
			labelsize += long_clab[i].size()*8;

			for (unsigned& elem:unsigned_lab[i]){
				split32(elem, vid, cc);
				if (cc > 0) labelsize += 4;
			}
			for (long& elem:long_lab[i]){
				split64(elem, vid, cc);
				if (cc > 0) labelsize += 8;
			}

			vector<unsigned>().swap(unsigned_lab[i]);
			vector<unsigned>().swap(unsigned_clab[i]);
			vector<long>().swap(long_lab[i]);
			vector<long>().swap(long_clab[i]);
		}

		#ifdef __linux__
			malloc_trim(0);
		#endif

		BoundGraph();
		worker_barrier();
		float t1 = omp_get_wtime();
		cout<<_my_rank<<"  "<<t1-t<<"  s"<<endl;

		long long total_size = all_sum_LL(labelsize);
		
		if (_my_rank == 0){
			cout<<float(total_size)/(1024*1024)<<" MB"<<endl;
			float t2 = omp_get_wtime();
			Core2hop();
			float t3 = omp_get_wtime();
			cout<<t3-t2<<"  s"<<endl;

			// new_filename = params.output_path + "Bd";
			// DHBound_Store(new_filename);

			for (int i=0; i<totalV; ++i){
				total_size += unsigned_lab[i].size()*4;
				total_size += unsigned_clab[i].size()*4;
				total_size += long_lab[i].size()*8;
				total_size += long_clab[i].size()*8;
			}

			cout<<"Label size: "<<float(total_size)/(1024*1024)<<" MB"<<endl;
		}
			
	
	}
};



void SetLCR_Construction(string in_path, string partition_path, string out_path, int src, int dst, int khop){
	WorkerParams param;
	param.input_path=in_path;
	param.partition_path = partition_path;
	param.output_path=out_path;
	param.src = src;
	param.dst = dst;
	param.khop = khop;
	param.force_write=true;
	param.native_dispatcher=false;
	LCRBlockWorker worker;
	worker.set_compute_mode(LCRBlockWorker::VB_COMP);

	worker.run_LCR(param);
};

