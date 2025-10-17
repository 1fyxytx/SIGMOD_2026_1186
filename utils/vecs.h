#ifndef VECS_H_
#define VECS_H_

#include "serialization.h"
#include "Combiner.h"
#include "global.h"
#include <vector>
using namespace std;

template <class KeyT, class MessageT>
struct msgpair {
    KeyT key;
    MessageT msg;

    msgpair()
    {
    }

    msgpair(KeyT v1, MessageT v2)
    {
        key = v1;
        msg = v2;
    }

    inline bool operator<(const msgpair& rhs) const
    {
        return key < rhs.key;
    }
};

template <class KeyT, class MessageT>
ibinstream& operator<<(ibinstream& m, const msgpair<KeyT, MessageT>& v)
{
    m << v.key;
    m << v.msg;
    return m;
}

template <class KeyT, class MessageT>
obinstream& operator>>(obinstream& m, msgpair<KeyT, MessageT>& v)
{
    m >> v.key;
    m >> v.msg;
    return m;
}

//===============================================

template <class KeyT, class MessageT, class HashT>
class Vecs {
public:
    typedef vector<msgpair<KeyT, MessageT> > Vec;
    typedef vector<Vec> VecGroup;

    int np;
    VecGroup vecs;
    HashT hash;

    Vecs()
    {
        int np = _num_workers;
        this->np = np;
        vecs.resize(np);
    }

    void append(const KeyT key, const MessageT msg)
    {
        msgpair<KeyT, MessageT> item(key, msg);
        vecs[hash(key)].push_back(item);
    }

    void append_ECPG(const KeyT key, const int wid, const MessageT msg)
    {
        msgpair<KeyT, MessageT> item(key, msg);
        vecs[wid].push_back(item);
    }

    Vec& getBuf(int pos)
    {
        return vecs[pos];
    }

    VecGroup& getBufs()
    {
        return vecs;
    }

    void clear()
    {
        for (int i = 0; i < np; i++) {
            vecs[i].clear();
        }
    }

    //============================
    //apply combiner logic

    void combine()
    {
        Combiner<MessageT>* combiner = (Combiner<MessageT>*)get_combiner();
        for (int i = 0; i < np; i++) {
            sort(vecs[i].begin(), vecs[i].end()); // 这里排序很重要
            Vec newVec;
            int size = vecs[i].size(); // 发送给第i个进程的消息的数量
            if (size > 0) {
                newVec.push_back(vecs[i][0]);
                KeyT preKey = vecs[i][0].key; // 找到这个消息具体发送给哪个顶点
                for (int j = 1; j < size; j++) { // 拿到第一个消息之后，然后遍历下面的信息
                    msgpair<KeyT, MessageT>& cur = vecs[i][j];
                    if (cur.key != preKey) {
                        newVec.push_back(cur);
                        preKey = cur.key;
                    } else {
                        combiner->combine(newVec.back().msg, cur.msg);
                        // 针对发送给相同顶点的信息进行合并
                    }
                }
            }
            newVec.swap(vecs[i]); // 把每个进程合并后的新信息与原本的信息进行交换
        }
    }

    long long get_total_msg()
    {
        long long sum = 0;
        for (int i = 0; i < vecs.size(); i++) {
            sum += vecs[i].size();
        }
        return sum;
    }
};

#endif
