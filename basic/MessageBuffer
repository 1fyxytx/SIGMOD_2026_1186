#ifndef MESSAGEBUFFER_H
#define MESSAGEBUFFER_H

#include <vector>
#include "../utils/global.h"
#include "../utils/Combiner.h"
#include "../utils/communication.h"
#include "../utils/vecs.h"
using namespace std;

template <class VertexT>
class MessageBuffer {
public:
    typedef typename VertexT::KeyType KeyT;
    typedef typename VertexT::MessageType MessageT;
    typedef typename VertexT::HashType HashT;
    typedef vector<MessageT> MessageContainerT;
    typedef hash_map<KeyT, int> Map; //int = position in v_msg_bufs //CHANGED FOR VADD
    typedef Vecs<KeyT, MessageT, HashT> VecsT;
    typedef typename VecsT::Vec Vec;
    typedef typename VecsT::VecGroup VecGroup;
    typedef typename Map::iterator MapIter;

    VecsT out_messages;
    Map in_messages;
    vector<VertexT*> to_add;
    vector<MessageContainerT> v_msg_bufs;
    HashT hash;

    void init(vector<VertexT*> vertexes)
    {
        v_msg_bufs.resize(vertexes.size()); // 针对自身进程内的每个顶点
        for (int i = 0; i < vertexes.size(); i++) {
            VertexT* v = vertexes[i];
            in_messages[v->id] = i; //也是针对自身进程的所有顶点
        }
    }

    void reinit(vector<VertexT*> vertexes)
    {
        v_msg_bufs.resize(vertexes.size());
        in_messages.clear();
        for (int i = 0; i < vertexes.size(); i++) {
            VertexT* v = vertexes[i];
            in_messages[v->id] = i; //CHANGED FOR VADD
        }
    }

    void add_message(const KeyT& id, const MessageT& msg)
    {
        hasMsg(); //cannot end yet even every vertex halts
        out_messages.append(id, msg);
    }

    void add_message_ECPG(const KeyT& id, const int wid, const MessageT& msg)
    {
        hasMsg(); //cannot end yet even every vertex halts
        out_messages.append_ECPG(id, wid, msg);
    }

    Map& get_messages()
    {
        return in_messages;
    }

    void combine()
    {
        //apply combiner
        Combiner<MessageT>* combiner = (Combiner<MessageT>*)get_combiner();
        if (combiner != NULL)
            out_messages.combine();
    }

    vector<VertexT*>& sync_messages()
    {	// to_add应该是要先自身把新添的那部分加进去，然后再分到对应的分区
        int np = get_num_workers();
        int me = get_worker_id();

        // get messages from remote
        vector<vector<VertexT*> > add_buf(_num_workers);
        for (int i = 0; i < to_add.size(); i++) {
            VertexT* v = to_add[i];
            add_buf[hash(v->id)].push_back(v);
        }
        // ================================================
        // exchange msgs
        // exchange vertices to add

        all_to_all_cat(out_messages.getBufs(), add_buf);


        //------------------------------------------------
        //delete sent vertices
        for (int i = 0; i < to_add.size(); i++) {
            VertexT* v = to_add[i];
            if (hash(v->id) != _my_rank)
                delete v;
        }
        // 这里的to_add记录的是一些要发送到相邻点的信息，若是相邻点不在自身的进程上，那么对于该进程是消耗内存的，所以需要删除
        to_add.clear();
        //collect new vertices to add
        for (int i = 0; i < np; i++) {
            to_add.insert(to_add.end(), add_buf[i].begin(), add_buf[i].end());
        }
        //================================================
        //Change of G33
        int oldsize = v_msg_bufs.size(); // 等于自身进程现有的顶点数量
        v_msg_bufs.resize(oldsize + to_add.size()); // to_add应该是外部迁移进来的顶点，或者说是动态图新加的顶点

        for (int i = 0; i < to_add.size(); i++) {
            int pos = oldsize + i;
            in_messages[to_add[i]->id] = pos; //CHANGED FOR VADD，in_messages 是个 map 格式
        }


        for (int i = 0; i < np; i++) {
            Vec& msgBuf = out_messages.getBuf(i); // msgBuf 是自身进程从进程i处获得的信息
            for (int i = 0; i < msgBuf.size(); i++) {
                MapIter it = in_messages.find(msgBuf[i].key); // it->first 直接对应的是 顶点的id, it->second是顶点在容器中的位置
                if (it != in_messages.end()){ //filter out msgs to non-existent vertices
                    v_msg_bufs[it->second].push_back(msgBuf[i].msg); //CHANGED FOR VADD
                }
            }
        }

        out_messages.clear();

        return to_add;
    }

    vector<VertexT*>& sync_messages_ECPG()
	{	// to_add应该是要先自身把新添的那部分加进去，然后再分到对应的分区
		int np = get_num_workers();

		for (int i = 0; i < np; i++) {
			Vec& msgBuf = out_messages.getBuf(i); // msgBuf 是自身进程从进程i处获得的信息
			for (int i = 0; i < msgBuf.size(); i++) {
				MapIter it = in_messages.find(msgBuf[i].key); // it->first 直接对应的是 顶点的id, it->second是顶点在容器中的位置
				if (it != in_messages.end()){ //filter out msgs to non-existent vertices
					v_msg_bufs[it->second].push_back(msgBuf[i].msg); //CHANGED FOR VADD
				}
			}
		}

		out_messages.clear();

		return to_add;
	}

    void add_vertex(VertexT* v)
    {
        hasMsg(); //cannot end yet even every vertex halts
        // 这个不是用来做顶点分区迁移的，或者说论文作者没想过这个功能
        to_add.push_back(v);
    }

    long long get_total_msg()
    {
        return out_messages.get_total_msg();
    }

    int get_total_vadd()
    {
        return to_add.size();
    }

    vector<MessageContainerT>& get_v_msg_bufs()
    {
        return v_msg_bufs;
    }
};

#endif
