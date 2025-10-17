#ifndef GLOBAL_H
#define GLOBAL_H

#include <mpi.h>
#include <stddef.h>
#include <limits.h>
#include <string>
#include <ext/hash_set>
#include <ext/hash_map>
#include <map>
#define hash_map __gnu_cxx::hash_map
#define hash_set __gnu_cxx::hash_set
#include <assert.h> //for ease of debug
#include<iostream>
#include <utility>
#include <unordered_map>

using namespace std;

//============================
///worker info
#define MASTER_RANK 0

int _my_rank;
int _num_workers;
inline int get_worker_id()
{
    return _my_rank;
}
inline int get_num_workers()
{
    return _num_workers;
}

void init_workers()
{
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &_num_workers);
    MPI_Comm_rank(MPI_COMM_WORLD, &_my_rank);
}

void worker_finalize()
{
    MPI_Finalize();
}

void worker_barrier()
{
    MPI_Barrier(MPI_COMM_WORLD);
}

//------------------------
// worker parameters

struct WorkerParams {
    string input_path;
    string output_path;
    string partition_path;
    string index_path;
    int src;
    int dst;
    int khop;

    bool force_write;
    bool native_dispatcher; //true if input is the output of a previous blogel job

    WorkerParams()
    {
        force_write = true;
        native_dispatcher = false;
    }
};

struct MultiInputParams {
    vector<string> input_paths;
    string output_path;
    bool force_write;
    bool native_dispatcher; //true if input is the output of a previous blogel job

    MultiInputParams()
    {
        force_write = true;
        native_dispatcher = false;
    }

    void add_input_path(string path)
    {
        input_paths.push_back(path);
    }
};

//============================
//general types
typedef int VertexID;
typedef int WorkerID;
//============================
//global variables
int global_step_num;
inline int step_num()
{
    return global_step_num;
}

int global_phase_num;
inline int phase_num()
{
    return global_phase_num;
}

void* global_message_buffer = NULL;
inline void set_message_buffer(void* mb)
{
    global_message_buffer = mb;
}
inline void* get_message_buffer()
{
    return global_message_buffer;
}

void* global_combiner = NULL;
inline void set_combiner(void* cb)
{
    global_combiner = cb;
}
inline void* get_combiner()
{
    return global_combiner;
}

void* global_aggregator = NULL;
inline void set_aggregator(void* ag)
{
    global_aggregator = ag;
}
inline void* get_aggregator()
{
    return global_aggregator;
}

void* global_agg = NULL; //for aggregator, FinalT of last round
inline void* getAgg()
{
    return global_agg;
}

int global_vnum = 0;
inline int& get_vnum()
{
    return global_vnum;
}

int global_edge = 0;
inline int& get_EdgeNum()
{
	return global_edge;
}


int global_active_vnum = 0;
inline int& active_vnum()
{
    return global_active_vnum;
}

enum BITS {
    HAS_MSG_ORBIT = 0,
    FORCE_TERMINATE_ORBIT = 1,
    WAKE_ALL_ORBIT = 2
};
//currently, only 3 bits are used, others can be defined by users
char global_bor_bitmap;

void clearBits()
{
    global_bor_bitmap = 0;
}

void setBit(int bit)
{
    global_bor_bitmap |= (2 << bit);
}

int getBit(int bit, char bitmap)
{
    return ((bitmap & (2 << bit)) == 0) ? 0 : 1;
}

void hasMsg()
{
    setBit(HAS_MSG_ORBIT);
}

void wakeAll()
{
    setBit(WAKE_ALL_ORBIT);
}

void forceTerminate()
{
    setBit(FORCE_TERMINATE_ORBIT);
}

//====================================================
//Ghost threshold
int global_ghost_threshold;

void set_ghost_threshold(int tau)
{
    global_ghost_threshold = tau;
}

//====================================================
string to_String(int n)
{
	if (n==0)
		return "0";
	int max = 100;
    int m = n;
    char s[max];
    char ss[max];
    int i=0,j=0;
    if (n < 0)// 处理负数
    {
        m = 0 - m;
        j = 1;
        ss[0] = '-';
    }
    while (m>0)
    {
        s[i++] = m % 10 + '0';
        m /= 10;
    }
    s[i] = '\0';
    i = i - 1;
    while (i >= 0)
    {
        ss[j++] = s[i--];
    }
    ss[j] = '\0';
    return ss;
}
//==============  ML2Hop 全局变量   自身分区所有顶点的入度和出度标签  ==============
//  { 源vid=标签vid, wid=imp }
vector<string> nOutLabel;
vector<string> nInLabel;

vector<vector<string> > NewOutL; // 每轮增加的标签， 这部分需要经过机器全通信
vector<vector<string> > NewInL;  //

//vector<vector<string> > Inner_Vertices_Inf; // 记录索引信息
vector<string> Inner_Vertices_Inf; // 记录索引信息


// 源id,  标签id>, 分区id
map<int, vector<int> > TotalOutL;   // 总的边界点标签，需要根据迭代更新的标签来进行检测
map<int, vector<int> > TotalInL;    //



int gloab_receive_size = 0;

set<int> ML_Top_Vertices;

unordered_map<int, int> vert2place, BVplace;  //  顶点id-顶点在vertexes中的位置                每个worker内，边界点在    vertexes 中的位置
//unordered_map
vector<string> QueryInf;
vector<vector<string> > QueryInfVector;

unordered_map<int, pair<int, int> > Bound2Imp2Wid;
vector<pair<int, int> > Imp2Vid; // 本地的排序列表
set<int> ActiveVSet; // 边界点计算的激活顶点集合
vector<map<int, vector<int> > > transIn;
set<int> newIn, newOut; // 针对单一顶点获得的新的索引信息
//========================测试结构===================================

set<int> srcSet;
set<int> dstSet;
set<int> dstWorker;

map<pair<int, int>, int> InitialResult; // 统计最后结果
vector<map<pair<int, int>, int> > ReachResults; // 可达的查询结果

//==========  公共查询函数    ==========
void SetInsert(int n, int num){
	for (int i=0; i<n; i++){
		if (i > n-num)
			dstSet.insert(i);
		if (i < num)
			srcSet.insert(i+100);
	}
}

void SetInsert1(int n1, int n2, int num1, int num2){
	for (int i=0; i<n1; i++){
		if (i > n1-num1 && i < n1-num1 + n2)
			dstSet.insert(i);
		if (i < num1 && i > num2)
			srcSet.insert(i+10);
	}
}
// batch query
void BatchQuerySet(){
	for(set<int>::iterator it1=srcSet.begin(); it1!=srcSet.end(); ++it1)
		for(set<int>::iterator it2=dstSet.begin(); it2!=dstSet.end(); ++it2){
			pair<int, int> keys(*it1, *it2);
			InitialResult[keys] = 0; // 查询结果全部初始化为0
		}
}

int spvalue = 1;

//========  索引更新算法全局变量    =============
int num_Insert = 5; // 每个分区随机选 num_Insert个顶点作为可能的更新点
map<int, int> InsertMap; //
vector<vector<string> > InsertVector;
vector<string> InsertInf;
// 第一个pair存储的是 源点的id和wid， 后面的pair，第一个是新加的出度点集， v1-w1-v2-w2  第二个是入度点集
map<int, vector<int> > Insert_Structure_Out;
map<int, vector<int> > Insert_Structure_In;
float percent = 1; // 加边的一个判断阈值

//---------
int minBoundImp = -1;
vector<int> minImpVect;
int UpdateFlag = 0;

//---------
map<int, int> Inner2Bound; // 所有从内部点变成边界点的id-imp的集合
set<int> InnerVert, BoundVert; //没有

set<int> virtual_id_set;
map<int, int> Source2Wid; // 只记录自身分区非循环的点
map<int, set<int> > Source2DestSet, AddElem;
set<int> StopElements; // 统计可以终止的源点信号
int global_reachnum = 0;

int bound_num = 0;


//---------
string ML2hop_index_name = "/vol6/home/stu_zyy/blogel/ML2hop/data/bound_data/";


// ===========  TVL 超参数     ============
int maxSize = 10; // 每个顶点所能获得的最大标签数量
int boundNum = 0;




// ===========  dsr 查询结构    ============
map<pair<int, int>, pair<int, int> > vid2Inf; // (vid, wid) — (if_In, if_Out)
map<int, set<int> > wId2vId;

vector<string> Related_Inf, Neig_Inf; // vid-wid-ifIn-ifOut, 出度信息=入度信息
vector<vector<string> > ALL_Inf, ALL_Neig_Inf;






















//==============    ECPG 定义变量    ======================
vector<double> vNum2Worker;// 记录每个进程上的顶点数
vector<double> edge2Worker;// 记录每个进程上的边的数量

double global_edges = 0;
double global_alpha_ECPG = 0.92;
double global_beta_ECPG = 0.0;
double totalVertex = 0;
double totalEdge = 0;


double delta = 0.15;
int global_changeNum = 0;
double global_degree = 0;

// ===========  IncKGGGP参数  =============
map<pair<int, int>, vector<int> > scNonBound; // 统计 <scvalue, id>, [part]
map<pair<int, int>, vector<int> > scBound; // 边界点总是先分配的，这个优先级更高
vector<vector<pair<int, int> > > topKpart;
vector<int> part_flag;
vector<int> part_bound;
vector<int> part_num;

map<pair<int, int>, int> LocalVert; // 所有分配到本分区的顶点 <sc, id>,需要频繁更新
map<int, int> topKplace; // 存放topK个对应的id和分区id, id具有唯一性
map<int, int> deleteVert; // 由于影响了balance，从而需要被删除的点
int topK = 20000; // 每次最多给topK个
double e_balance = 0.2;
int partition_flag = 0;
int maxNum, numCount;


set<int> Vr;




// ===========  flag  =================
int Insert_or_Delete = 2; //1-模拟加边；0-表示模拟减边；-1表示txt中读入streaming数据流

// ===========  新插入顶点  =================
int NewVert_Num = 10000; // 新插入顶点的数量
int newDegree = 200; // 新插入顶点的degree
int Insertion_Flag = 1;
set<int> Influence_Vertices;
map<pair<int, int>, pair<map<int, int>, map<int, int> > > Insert_Structure;

// ===========  删除顶点   ==================
int DeleteVert_Num = 100;
vector<vector<string> > DeleteVector; // d1=v1-v2-v3=v4-v5-v6  其中d1是被删除的点，v1到v3是相关的邻点
vector<string> DeleteInf;
map<int, pair<set<int>, set<int> > > Delete_Structure; // 第一个set是出度, 第二个set是入度


// ======= streaming情形下的参数设置  ===========
vector<map<int, int> > global_v2w;
map<int, int> global_vid2wid;




// ========= Fennel and LGD ==========
int if_Fennel = 0; // 0-对应 LGD 1-对应 Fennel
double beta = 0.1;
double alpha = 0.5, lamda = 2;
set<int> Vt_Set; // <vid, wid>





#define ROUND 11 //for PageRank

#endif
