#include "../utils/communication.h"
#include "../blogel/BVertex.h"
#include "../blogel/Block.h"
#include "../blogel/BWorker.h"
#include "../blogel/BGlobal.h"
#include "../blogel/BType.h"
#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <ext/hash_map>
#include <string.h>
#include <map>
#include <numeric>
#include <utility>
#include <time.h>
#include <deque>
#include <set>
#include <algorithm>
#include <bitset>

using namespace std;


struct STLCRVertexValue {
	vector<pair<unsigned, unsigned> > indL, indE; // 两种index
	int src, dst;

    // void set_Edge(unsigned val, int flg){

	// 	if (flg == 1)
	// 		indE.push_back(val);
	// 	else
	// 		indL.push_back(val);
    // }

	void Initial(){
		src = 1000, dst = -1000;
	}

	void empty(){
		vector<pair<unsigned, unsigned> >().swap(indL);
		vector<pair<unsigned, unsigned> >().swap(indE);
	}

    friend ibinstream & operator<<(ibinstream & m, const STLCRVertexValue & v){
		m<<v.indL;
		m<<v.indE;
		m<<v.src;
		m<<v.dst;

    	return m;
    }


    friend obinstream & operator>>(obinstream & m, STLCRVertexValue & v){
		m>>v.indL;
    	m>>v.indE;
		m>>v.src;
		m>>v.dst;

		return m;
    }
};


struct MsgInf {
	int dis;

	MsgInf(){dis = -1;}
    
	MsgInf(int d){dis = d;}

    friend ibinstream& operator<<(ibinstream& m, const MsgInf& idm)
    {
        m << idm.dis;

        return m;
    }

    friend obinstream& operator>>(obinstream& m, MsgInf& idm)
    {
        m >> idm.dis;

        return m;
    }
};


class LCRVertex : public BVertex<VertexID, STLCRVertexValue, MsgInf>{
public:
	virtual void compute(MessageContainer& messages){}

	void Inner_Compute(MessageContainer& messages){
		if (step_num() == 1){ // 添加元素，剩下的在block里面做
			if (id == src) 
				value().src = 0;
			if (id == dst) 
				value().dst = 0;
		}else{ // 第二个超步，将message中的顶点添加到send中就完成
			for (MessageIter it = messages.begin(); it != messages.end(); it++) {
				if ((*it).dis > 0 and value().src > (*it).dis)
					value().src = (*it).dis;
				
				if ((*it).dis < 0 and value().dst < (*it).dis)
					value().dst = (*it).dis;
			}
			
			if (value().src < 1000 and value().dst > -1000){
				if (khop > value().src-value().dst)
					khop = value().src-value().dst;
			}
		}

		vote_to_halt();
    }
};


class LCRBlock : public Block<char, LCRVertex, MsgInf> {
public:
	//
	virtual void compute(MessageContainer &messages, VertexContainer &vertexes){}
		// == 先执行 内部的传递 ==

};


class LCRBlockWorker : public BWorker<LCRBlock>{
public:
	vector<vector<vector<unsigned> > > Lab1, Lab2;
	vector<int> D1, D2;

	vector<unsigned> *unsigned_lab, *unsigned_clab;
	vector<long> *long_lab, *long_clab;
	vector<int> *int_pos, *int_cpos, *long_pos, *long_cpos;

	vector<unsigned> v2p, p2v, v2p_Global, p2v_Global;
	vector<pair<int, int> > QueryPairs;

	vector<int> ActiveV;

	vector<pair<unsigned, unsigned> > Psrc, Pdst, Csrc, Cdst;
	long long InSize = 0, BoundSize = 0;
	int id_bit = 0, dis_bit = 0, cnt_bit = 0, max_id_bit = 0;

	int BoundNum = 0, BdV = 0, candvalue = 0;
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

	virtual void blockInit(VertexContainer& vertexList, BlockContainer& blockList){
		
		vector<int> srcL, dstL;
		int start = BoundNum, mid = 10000;
		for (int i=start; i<start+mid; ++i){
			if (_my_rank % 2 == 1)  srcL.push_back(p2v[i]);
			else                    dstL.push_back(p2v[i]);
		}

		// start = 0;
		// for (int i=start; i<start+mid; ++i){
		// 	if (_my_rank % 2 == 0)  dstL.push_back(p2v[i]);
		// }

		vector<vector<int> > SrcL(_num_workers), DstL(_num_workers);
		for (int i=0; i<_num_workers; ++i){
			SrcL[i] = srcL; DstL[i] = dstL;
		}

		all_to_all_cat(SrcL, DstL);
		worker_barrier();

		srcL.clear(), dstL.clear();

		for (int i=0; i<_num_workers; ++i){
			srcL.insert(srcL.end(), SrcL[i].begin(), SrcL[i].end());
			dstL.insert(dstL.end(), DstL[i].begin(), DstL[i].end());
		}

		vector<vector<int> >().swap(SrcL);
		vector<vector<int> >().swap(DstL);

		sort(srcL.begin(), srcL.end());
		sort(dstL.begin(), dstL.end());

		for (int i=0; i<srcL.size(); ++i){
			if (srcL[i] != dstL[i] and v2part[srcL[i]] != v2part[dstL[i]]){
				pair<int, int> elem(srcL[i], dstL[i]);
				QueryPairs.push_back(elem);
			}
		}
		// =================


	}


	void active_LCR_vcompute(){ }

	void all_LCR_vcompute(){ }


	void all_bcompute_LCR(){
		active_bcount = 0;
		BMessageBufT* mbuf = (BMessageBufT*)get_bmessage_buffer();
		vector<BMessageContainerT>& b_msgbufs = mbuf->get_b_msg_bufs();
		for (int i = 0; i < blocks.size(); i++){
			blocks[i]->activate();
			blocks[i]->compute(b_msgbufs[i], vertexes);
			b_msgbufs[i].clear(); //clear used msgs
			if (blocks[i]->is_active())
				active_bcount++;
		}
	};

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

	void InLRead(string indexName){
		
		int p2v_cnt = 0;

		// ==  ActiveV[_my_rank] 代表本分区顶点的数量  ==
		LabelInitial(ActiveV[_my_rank]); 
		
		ifstream fin((indexName).c_str(), ios::binary);

		fin.read(reinterpret_cast<char*>(&BoundNum), sizeof(unsigned));

		// ===== 读取 p2v =====
		fin.read(reinterpret_cast<char*>(&p2v_cnt), sizeof(unsigned));
		p2v.resize(p2v_cnt);
		fin.read(reinterpret_cast<char*>(p2v.data()), p2v_cnt * sizeof(unsigned)); // 读内层内容

		fin.read(reinterpret_cast<char*>(&p2v_cnt), sizeof(unsigned));
		v2p.resize(p2v_cnt);
		fin.read(reinterpret_cast<char*>(v2p.data()), p2v_cnt * sizeof(unsigned)); // 读内层内容


		// =====  load labels  ===== 
		for (int i=BoundNum; i<ActiveV[_my_rank]; i++){
			unsigned labcnt = 0;
			int poscnt = 0;
			
			fin.read(reinterpret_cast<char*>(&labcnt), sizeof(unsigned));
			unsigned_lab[i].resize(labcnt);
			fin.read(reinterpret_cast<char*>(unsigned_lab[i].data()), labcnt * sizeof(unsigned)); // 读内层内容
			
			fin.read(reinterpret_cast<char*>(&poscnt), sizeof(int));
			int_pos[i].resize(poscnt);
			fin.read(reinterpret_cast<char*>(int_pos[i].data()), poscnt * sizeof(int)); // 读内层内容

			InSize += labcnt*4;
			// ==========================

			fin.read(reinterpret_cast<char*>(&labcnt), sizeof(unsigned));
			unsigned_clab[i].resize(labcnt);
			fin.read(reinterpret_cast<char*>(unsigned_clab[i].data()), labcnt * sizeof(unsigned)); // 读内层内容
			
			fin.read(reinterpret_cast<char*>(&poscnt), sizeof(int));
			int_cpos[i].resize(poscnt);
			fin.read(reinterpret_cast<char*>(int_cpos[i].data()), poscnt * sizeof(int)); // 读内层内容

			InSize += labcnt*4;

			// ==========================
			fin.read(reinterpret_cast<char*>(&labcnt), sizeof(long));
			long_lab[i].resize(labcnt);
			fin.read(reinterpret_cast<char*>(long_lab[i].data()), labcnt * sizeof(long)); // 读内层内容
			
			fin.read(reinterpret_cast<char*>(&poscnt), sizeof(int));
			long_pos[i].resize(poscnt);
			fin.read(reinterpret_cast<char*>(long_pos[i].data()), poscnt * sizeof(int)); // 读内层内容

			InSize += labcnt*8;

			// ==========================
			fin.read(reinterpret_cast<char*>(&labcnt), sizeof(long));
			long_clab[i].resize(labcnt);
			fin.read(reinterpret_cast<char*>(long_clab[i].data()), labcnt * sizeof(long)); // 读内层内容
			
			fin.read(reinterpret_cast<char*>(&poscnt), sizeof(int));
			long_cpos[i].resize(poscnt);
			fin.read(reinterpret_cast<char*>(long_cpos[i].data()), poscnt * sizeof(int)); // 读内层内容

			InSize += labcnt*8;
		}

		fin.close();
	}
	
	void skipVectorData(ifstream& fin, size_t elementSize) {
		unsigned count = 0;
		fin.read(reinterpret_cast<char*>(&count), sizeof(unsigned));
		fin.seekg(count * elementSize, ios::cur);
	}

	void BdLRead(string indexName){
		int p2v_cnt = 0;
		
		ifstream fin((indexName).c_str(), ios::binary);

		fin.read(reinterpret_cast<char*>(&BdV), sizeof(unsigned));

		// ===== 读取 p2v =====
		fin.read(reinterpret_cast<char*>(&p2v_cnt), sizeof(unsigned));
		p2v_Global.resize(p2v_cnt);
		fin.read(reinterpret_cast<char*>(p2v_Global.data()), p2v_cnt * sizeof(unsigned)); // 读内层内容

		fin.read(reinterpret_cast<char*>(&p2v_cnt), sizeof(unsigned));
		v2p_Global.resize(p2v_cnt);
		fin.read(reinterpret_cast<char*>(v2p_Global.data()), p2v_cnt * sizeof(unsigned)); // 读内层内容

		// =====  load labels  ===== ActiveV[_my_rank]
		for (int i=0; i<BdV; i++){
			unsigned labcnt = 0, poscnt = 0;
			
			int orig = p2v_Global[i], pos = v2p[orig]; // pos对应在子图内部的位置

			// === 不满足要求的，直接不读取 ==

			if (v2part[orig] != _my_rank){
				skipVectorData(fin, sizeof(unsigned));  // unsigned_lab
				skipVectorData(fin, sizeof(int));       // int_pos

				skipVectorData(fin, sizeof(unsigned));  // unsigned_lab
				skipVectorData(fin, sizeof(int));       // int_pos

				skipVectorData(fin, sizeof(long));  // unsigned_lab
				skipVectorData(fin, sizeof(int));       // int_pos

				skipVectorData(fin, sizeof(long));  // unsigned_lab
				skipVectorData(fin, sizeof(int));       // int_pos
			}else{
				fin.read(reinterpret_cast<char*>(&labcnt), sizeof(unsigned));
				unsigned_lab[pos].resize(labcnt);
				fin.read(reinterpret_cast<char*>(unsigned_lab[pos].data()), labcnt * sizeof(unsigned)); // 读内层内容
				
				InSize += labcnt*4;

				fin.read(reinterpret_cast<char*>(&poscnt), sizeof(int));
				int_pos[pos].resize(poscnt);
				fin.read(reinterpret_cast<char*>(int_pos[pos].data()), poscnt * sizeof(int)); // 读内层内容

				// ====
				fin.read(reinterpret_cast<char*>(&labcnt), sizeof(unsigned));
				unsigned_clab[pos].resize(labcnt);
				fin.read(reinterpret_cast<char*>(unsigned_clab[pos].data()), labcnt * sizeof(unsigned)); // 读内层内容
				
				InSize += labcnt*4;

				fin.read(reinterpret_cast<char*>(&poscnt), sizeof(int));
				int_cpos[pos].resize(poscnt);
				fin.read(reinterpret_cast<char*>(int_cpos[pos].data()), poscnt * sizeof(int)); // 读内层内容

				//  ===========
				fin.read(reinterpret_cast<char*>(&labcnt), sizeof(unsigned));
				long_lab[pos].resize(labcnt);
				fin.read(reinterpret_cast<char*>(long_lab[pos].data()), labcnt * sizeof(long)); // 读内层内容
				
				InSize += labcnt*8;

				fin.read(reinterpret_cast<char*>(&poscnt), sizeof(int));
				long_pos[pos].resize(poscnt);
				fin.read(reinterpret_cast<char*>(long_pos[pos].data()), poscnt * sizeof(int)); // 读内层内容

				// =============

				fin.read(reinterpret_cast<char*>(&labcnt), sizeof(unsigned));
				long_clab[pos].resize(labcnt);
				fin.read(reinterpret_cast<char*>(long_clab[pos].data()), labcnt * sizeof(long)); // 读内层内容
				
				InSize += labcnt*8;

				fin.read(reinterpret_cast<char*>(&poscnt), sizeof(int));
				long_cpos[pos].resize(poscnt);
				fin.read(reinterpret_cast<char*>(long_cpos[pos].data()), poscnt * sizeof(int)); // 读内层内容

			}
			

			// // ==========================

			// fin.read(reinterpret_cast<char*>(&labcnt), sizeof(unsigned));
			// unsigned_clab[pos].resize(labcnt);
			// fin.read(reinterpret_cast<char*>(unsigned_clab[pos].data()), labcnt * sizeof(unsigned)); // 读内层内容
			
			// fin.read(reinterpret_cast<char*>(&poscnt), sizeof(int));
			// int_cpos[pos].resize(poscnt);
			// fin.read(reinterpret_cast<char*>(int_cpos[pos].data()), poscnt * sizeof(int)); // 读内层内容

			// InSize += labcnt*4;

			// // ==========================
			// fin.read(reinterpret_cast<char*>(&labcnt), sizeof(long));
			// long_lab[pos].resize(labcnt);
			// fin.read(reinterpret_cast<char*>(long_lab[pos].data()), labcnt * sizeof(long)); // 读内层内容
			
			// fin.read(reinterpret_cast<char*>(&poscnt), sizeof(int));
			// long_pos[pos].resize(poscnt);
			// fin.read(reinterpret_cast<char*>(long_pos[pos].data()), poscnt * sizeof(int)); // 读内层内容

			// InSize += labcnt*8;

			// // ==========================
			// fin.read(reinterpret_cast<char*>(&labcnt), sizeof(long));
			// long_clab[pos].resize(labcnt);
			// fin.read(reinterpret_cast<char*>(long_clab[pos].data()), labcnt * sizeof(long)); // 读内层内容
			
			// fin.read(reinterpret_cast<char*>(&poscnt), sizeof(int));
			// long_cpos[pos].resize(poscnt);
			// fin.read(reinterpret_cast<char*>(long_cpos[pos].data()), poscnt * sizeof(int)); // 读内层内容

			// InSize += labcnt*8;
		}

		fin.close();

	}

	void split32(unsigned combined, unsigned& id, unsigned& cnt) {
		unsigned id_mask = (1U << max_id_bit) - 1;
		unsigned cnt_mask = (1U << (32 - max_id_bit)) - 1;
		
		id = (combined >> (32 - max_id_bit)) & id_mask;  // 提取高max_id_bit位作为id
		cnt = combined & cnt_mask;                       // 提取低(32-max_id_bit)位作为cnt
	}
	
	void split64(long combined, unsigned& id, unsigned& cnt) {
		id = static_cast<unsigned>((combined >> (64 - max_id_bit)) & ((1ULL << max_id_bit) - 1));
		cnt = static_cast<unsigned>(combined & ((1ULL << (64 - max_id_bit)) - 1));
	}

	void initial(int src, int dst){
		
		if (v2part[src] == _my_rank){
			
			unsigned vsrc = v2p[src];

			if (vsrc >= BoundNum){
				
				unsigned w, cc;
				int dis = int_pos[vsrc].size(), d = 0;

				while (d+1 < dis){
					for (int i=int_pos[vsrc][d]; i<int_pos[vsrc][d+1]; ++i){
						split32(unsigned_lab[vsrc][i], w, cc);
						unsigned elem = w<<6 | (d+1);
						Psrc.push_back(make_pair(elem, cc));
					}

					for (int i=int_cpos[vsrc][d]; i<int_cpos[vsrc][d+1]; ++i){
						split32(unsigned_clab[vsrc][i], w, cc);
						unsigned elem = w<<6 | (d+1);
						Psrc.push_back(make_pair(elem, cc));
					}

					for (int i=long_pos[vsrc][d]; i<long_pos[vsrc][d+1]; ++i){
						split64(long_lab[vsrc][i], w, cc);
						unsigned elem = w<<6 | (d+1);
						Psrc.push_back(make_pair(elem, cc));
					}

					for (int i=long_cpos[vsrc][d]; i<long_cpos[vsrc][d+1]; ++i){
						split64(long_clab[vsrc][i], w, cc);
						unsigned elem = w<<6 | (d+1);
						Psrc.push_back(make_pair(elem, cc));
					}

					d += 1;
				}
			}else{
				unsigned elem = vsrc<<6 | 0;
				Psrc.push_back(make_pair(elem, 1)); // self 只有1条
			}

		}

		
		if (v2part[dst] == _my_rank){
			
			unsigned vdst = v2p[dst];
			
			if (vdst >= BoundNum){
				
				unsigned w, cc;
				int dis = int_pos[vdst].size(), d = 0;

				while (d+1 < dis){
					for (int i=int_pos[vdst][d]; i<int_pos[vdst][d+1]; ++i){
						split32(unsigned_lab[vdst][i], w, cc);
						unsigned elem = w<<6 | (d+1);
						Pdst.push_back(make_pair(elem, cc));
					}

					for (int i=int_cpos[vdst][d]; i<int_cpos[vdst][d+1]; ++i){
						split32(unsigned_clab[vdst][i], w, cc);
						unsigned elem = w<<6 | (d+1);
						Pdst.push_back(make_pair(elem, cc));
					}

					for (int i=long_pos[vdst][d]; i<long_pos[vdst][d+1]; ++i){
						split64(long_lab[vdst][i], w, cc);
						unsigned elem = w<<6 | (d+1);
						Pdst.push_back(make_pair(elem, cc));
					}

					for (int i=long_cpos[vdst][d]; i<long_cpos[vdst][d+1]; ++i){
						split64(long_clab[vdst][i], w, cc);
						unsigned elem = w<<6 | (d+1);
						Pdst.push_back(make_pair(elem, cc));
					}

					d += 1;
				}
			}else{
				unsigned elem = vdst<<6 | 0;
				Pdst.push_back(make_pair(elem, 1)); // self 只有1条
			}
		}

	}

	void update(unsigned w, unsigned newdis, unsigned newcnt, 
				vector<pair<unsigned, unsigned>>& InfList,
			    vector<int>& pos){

		unsigned elem = w<<6 | newdis;

		if (pos[w] == -1){
			InfList.push_back(make_pair(elem, newcnt));
			pos[w] = InfList.size()-1;
		}else{
			if (elem < InfList[pos[w]].first){
				InfList[pos[w]].first = elem;
				InfList[pos[w]].second = newcnt;
			}else if(elem == InfList[pos[w]].first){
				InfList[pos[w]].second += newcnt;
			}
		}
	}

	unsigned CountCompute(vector<pair<unsigned, unsigned>>& Csrc, 
					vector<pair<unsigned, unsigned>>& Cdst, 
					unsigned& min_distance_sum) {
		unsigned total_sum = 0;
		min_distance_sum = 10000;
		int i = 0, j = 0;
		
		while (i < Csrc.size() && j < Cdst.size()) {
			unsigned src_id = Csrc[i].first >> 6;
			unsigned dst_id = Cdst[j].first >> 6;
			
			if (src_id < dst_id) {
				i++;
			} else if (src_id > dst_id) {
				j++;
			} else {
				// ID相同，计算距离和
				unsigned src_dis = Csrc[i].first & 0x3F;
				unsigned dst_dis = Cdst[j].first & 0x3F;
				unsigned distance_sum = src_dis + dst_dis;
				unsigned product = Csrc[i].second * Cdst[j].second;
				
				if (distance_sum < min_distance_sum) {
					// 找到更小的距离和，替换total_sum
					min_distance_sum = distance_sum;
					total_sum = product;
				} else if (distance_sum == min_distance_sum) {
					// 距离和相等，累加到total_sum
					total_sum += product;
				}
				
				i++;
				j++;
			}
		}
		
		return total_sum;
	}

	void run_LCR(const WorkerParams& params){
		ifstream infile, infile1, infile2; // 先把分区文件读取
		
		src = params.src, dst = params.dst;
		khop = 10000;

		const char *part_path = params.partition_path.c_str(); 
		infile.open(part_path);
		string ss;
		if(!infile.is_open()){
			cout<<"No such file!"<<endl;
			exit(-1);
		}

		int nnn = 0;
		ActiveV.resize(_num_workers);

		while(getline(infile, ss)){
			char* part = new char[strlen(ss.c_str())+1];
			strcpy(part, ss.c_str());
			v2part.push_back( atoi(part) );
			ActiveV[atoi(part)] += 1;
			v2degree.push_back(0);
			delete part;
			nnn += 1;
		}

		totalV = nnn;
		Lab1.resize(_num_workers), Lab2.resize(_num_workers);
        
		string indexName = params.input_path+"In_"+ to_string(_my_rank);
		InLRead(indexName);

		indexName = params.input_path+"Bd";
		BdLRead(indexName);

		long long lab = 0;
		for (int i=0; i<ActiveV[_my_rank]; ++i){
			lab += unsigned_lab[i].size() * 8;
			lab += unsigned_clab[i].size() * 8;
			lab += long_lab[i].size() * 8;
			lab += long_clab[i].size() * 8;
		}

		long long labelsize = all_sum_LL(lab);
		if (_my_rank == 0){
			cout<<labelsize/float(1024*1024*1024)<<" GB"<<endl;
		}
        worker_barrier();

		blockInit(vertexes, blocks); // === 随机生成查询任务队列 ===
		unsigned aa = countBits(BoundNum);
		max_id_bit = all_max(aa); // 子图内部的算子

		// // =================================================================
		float TotalTime = 0;
		long long TotalMessage = 0;
		long long global_msg_num = 0;
		long long msg_num = 0;

		if (_my_rank == 0)
			cout<<"Starting Query..."<<endl;
	
		for (int ii=0; ii<QueryPairs.size(); ++ii){
			
			vector<vector<pair<unsigned, unsigned> > > MsgSrc(_num_workers), MsgDst(_num_workers);

			src = QueryPairs[ii].first, dst = QueryPairs[ii].second;

			vector<int> v2d2cnt(BdV, -1), v2d2cnt_reverse(BdV, -1);
			
			// == 初始化 ==
			initial(src, dst);

			float s_1 = clock();

			for (pair<unsigned, unsigned>& elem: Psrc){
				unsigned w, cc;
				unsigned vid = elem.first >> 6, // 这里是1-hop，还需要传输给vid的label
					     dd = elem.first & 0x3F, cnt = elem.second; 

				int dis = int_pos[vid].size(), d = 0;

				while (d+1 < dis){
					
					for (int i=int_pos[vid][d]; i<int_pos[vid][d+1]; ++i){
						
						split32(unsigned_lab[vid][i], w, cc);
						
						unsigned part = v2part[p2v_Global[w]];
						unsigned newdis = dd+d+1, newcnt = cnt*cc;

						update(w, newdis, newcnt, MsgSrc[part], v2d2cnt);
					}

					for (int i=int_cpos[vid][d]; i<int_cpos[vid][d+1]; ++i){
						split32(unsigned_clab[vid][i], w, cc);
						
						unsigned part = v2part[p2v_Global[w]];
						unsigned newdis = dd+d+1, newcnt = cnt*cc;

						update(w, newdis, newcnt, MsgSrc[part], v2d2cnt);
					}

					for (int i=long_pos[vid][d]; i<long_pos[vid][d+1]; ++i){
						split64(long_lab[vid][i], w, cc);
						unsigned part = v2part[p2v_Global[w]];
						unsigned newdis = dd+d+1, newcnt = cnt*cc;

						update(w, newdis, newcnt, MsgSrc[part], v2d2cnt);
					}

					for (int i=long_cpos[vid][d]; i<long_cpos[vid][d+1]; ++i){
						split64(long_clab[vid][i], w, cc);
						unsigned part = v2part[p2v_Global[w]];
						unsigned newdis = dd+d+1, newcnt = cnt*cc;

						update(w, newdis, newcnt, MsgSrc[part], v2d2cnt);
					}

					d += 1;
				}
			}

			for (pair<unsigned, unsigned>& elem: Pdst){
				unsigned w, cc;
				unsigned vid = elem.first >> 6, // 这里是1-hop，还需要传输给vid的label
					     dd = elem.first & 0x3F, cnt = elem.second; 

				int dis = int_pos[vid].size(), d = 0;

				while (d+1 < dis){
					
					for (int i=int_pos[vid][d]; i<int_pos[vid][d+1]; ++i){
						
						split32(unsigned_lab[vid][i], w, cc);
						unsigned part = v2part[p2v_Global[w]];
						unsigned newdis = dd+d+1, newcnt = cnt*cc;

						update(w, newdis, newcnt, MsgDst[part], v2d2cnt_reverse);
					}

					for (int i=int_cpos[vid][d]; i<int_cpos[vid][d+1]; ++i){
						split32(unsigned_clab[vid][i], w, cc);
						
						unsigned part = v2part[p2v_Global[w]];
						unsigned newdis = dd+d+1, newcnt = cnt*cc;

						update(w, newdis, newcnt, MsgDst[part], v2d2cnt_reverse);

					}

					for (int i=long_pos[vid][d]; i<long_pos[vid][d+1]; ++i){
						split64(long_lab[vid][i], w, cc);
						unsigned part = v2part[p2v_Global[w]];
						unsigned newdis = dd+d+1, newcnt = cnt*cc;

						update(w, newdis, newcnt, MsgDst[part], v2d2cnt_reverse);					
					}

					for (int i=long_cpos[vid][d]; i<long_cpos[vid][d+1]; ++i){
						split64(long_clab[vid][i], w, cc);
						unsigned part = v2part[p2v_Global[w]];
						unsigned newdis = dd+d+1, newcnt = cnt*cc;

						update(w, newdis, newcnt, MsgDst[part], v2d2cnt_reverse);
					}

					d += 1;
				}
			}
			
			all_to_all_cat(MsgSrc, MsgDst);
			worker_barrier();

			for (int i=0; i<_num_workers; ++i){
				Csrc.insert(Csrc.end(), MsgSrc[i].begin(), MsgSrc[i].end());
				if (i != _my_rank){
					msg_num += MsgSrc[i].size();
				}
				// MsgSrc[i].clear();
				Cdst.insert(Cdst.end(), MsgDst[i].begin(), MsgDst[i].end());
				if (i != _my_rank){
					msg_num += MsgDst[i].size();
				}
				// MsgDst[i].clear();
			}

			unsigned mindis = 100000, cnt = 0;
			if (Csrc.size() > 0 and Cdst.size() > 0){
				sort(Csrc.begin(), Csrc.end());
				sort(Cdst.begin(), Cdst.end());
				cnt = CountCompute(Csrc, Cdst, mindis);		
				// cout<<_my_rank<<" "<<mindis<<" "<<cnt<<endl;
			}

			worker_barrier();
			float s_2 = clock();
			

			TotalTime += (s_2-s_1);
			global_msg_num = all_sum_LL(msg_num);
			// float(TotalTime)/(ii*CLOCKS_PER_SEC)
			if (_my_rank == 0){
				cout<<" DDH "<<ii<<"  "<<totallab<<"  Average: "<<float(s_2-s_1)/(CLOCKS_PER_SEC)<<" s , Comm: "<<float(global_msg_num*6)/(ii*1024)<<endl;
			}


			Psrc.clear(), Pdst.clear();
			Csrc.clear(), Cdst.clear();

			vector<vector<pair<unsigned, unsigned> > >().swap(MsgSrc);
			vector<vector<pair<unsigned, unsigned> > >().swap(MsgDst);
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

