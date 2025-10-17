#ifndef WORKER_H
#define WORKER_H

#include <vector>
#include "../utils/global.h"
#include "MessageBuffer.h"
#include <string.h>
#include "../utils/communication.h"
//#include "../utils/ydhdfs.h"
#include "../utils/Combiner.h"
#include "../utils/Aggregator.h"
using namespace std;

template <class VertexT, class AggregatorT = DummyAgg> //user-defined VertexT
class Worker {
public:
    typedef vector<VertexT*> VertexContainer;
    typedef typename VertexContainer::iterator VertexIter;

    typedef typename VertexT::KeyType KeyT;
    typedef typename VertexT::MessageType MessageT;
    typedef typename VertexT::HashType HashT;

    typedef MessageBuffer<VertexT> MessageBufT;
    typedef typename MessageBufT::MessageContainerT MessageContainerT;
    typedef typename MessageBufT::Map Map;
    typedef typename MessageBufT::MapIter MapIter;

    typedef typename AggregatorT::PartialType PartialT;
    typedef typename AggregatorT::FinalType FinalT;


    map<int, int> vertNumMap;
    Worker()
    {
        //init_workers();//put to run.cpp
        message_buffer = new MessageBuffer<VertexT>;
        global_message_buffer = message_buffer;
        active_count = 0;
        combiner = NULL;
        global_combiner = NULL;
        aggregator = NULL;
        global_aggregator = NULL;
        global_agg = NULL;
    }

    void setCombiner(Combiner<MessageT>* cb)
    {
        combiner = cb;
        global_combiner = cb;
    }

    void setAggregator(AggregatorT* ag)
    {
        aggregator = ag;
        global_aggregator = ag;
        global_agg = new FinalT;
    }

    virtual ~Worker()
    {
        for (int i = 0; i < vertexes.size(); i++)
            delete vertexes[i];
        delete message_buffer;
        if (getAgg() != NULL)
            delete (FinalT*)global_agg;
        worker_barrier(); //newly added for ease of multi-job programming in run.cpp
    }

    //==================================
    //sub-functions
    void sync_graph()
    {
        vector<VertexContainer> _loaded_parts(_num_workers);
        for (int i = 0; i < vertexes.size(); i++) {
            VertexT* v = vertexes[i];
            _loaded_parts[hash(v->id)].push_back(v);
        }
        //exchange vertices to add
        all_to_all(_loaded_parts);
        //delete sent vertices
        for (int i = 0; i < vertexes.size(); i++) {
            VertexT* v = vertexes[i];
            if (hash(v->id) != _my_rank)
                delete v;
        }
        vertexes.clear();

        //collect vertices to add
        for (int i = 0; i < _num_workers; i++) {
            vertexes.insert(vertexes.end(), _loaded_parts[i].begin(), _loaded_parts[i].end());
        }

        _loaded_parts.clear();
    };

    void active_compute()
    {
        active_count = 0;
        MessageBufT* mbuf = (MessageBufT*)get_message_buffer();
        vector<MessageContainerT>& v_msgbufs = mbuf->get_v_msg_bufs();
        for (int i = 0; i < vertexes.size(); i++) {
            if (v_msgbufs[i].size() == 0) {
                if (vertexes[i]->is_active()) {
                    vertexes[i]->compute(v_msgbufs[i]);
                    AggregatorT* agg = (AggregatorT*)get_aggregator();
                    if (agg != NULL)
                        agg->stepPartial(vertexes[i]);
                    if (vertexes[i]->is_active())
                        active_count++;
                }
            } else {
                vertexes[i]->activate();
                vertexes[i]->compute(v_msgbufs[i]);
                v_msgbufs[i].clear(); //clear used msgs
                AggregatorT* agg = (AggregatorT*)get_aggregator();
                if (agg != NULL)
                    agg->stepPartial(vertexes[i]);
                if (vertexes[i]->is_active())
                    active_count++;
            }
        }
    }

    void all_compute()
    {
        active_count = 0;
        MessageBufT* mbuf = (MessageBufT*)get_message_buffer();
        vector<MessageContainerT>& v_msgbufs = mbuf->get_v_msg_bufs();
        for (int i = 0; i < vertexes.size(); i++) {
            vertexes[i]->activate();
            vertexes[i]->compute(v_msgbufs[i]);
            v_msgbufs[i].clear(); //clear used msgs
            AggregatorT* agg = (AggregatorT*)get_aggregator();
            if (agg != NULL)
                agg->stepPartial(vertexes[i]);
            if (vertexes[i]->is_active())
                active_count++;
        }
    }

    inline void add_vertex(VertexT* vertex)
    {
        vertexes.push_back(vertex);
        if (vertex->is_active())
            active_count++;
    }

    void agg_sync()
    {
        AggregatorT* agg = (AggregatorT*)get_aggregator();
        if (agg != NULL) {
            if (_my_rank != MASTER_RANK) { //send partialT to aggregator
                //gathering PartialT
                PartialT* part = agg->finishPartial();
                //------------------------ strategy choosing BEGIN ------------------------
                StartTimer(COMMUNICATION_TIMER);
                StartTimer(SERIALIZATION_TIMER);
                ibinstream m;
                m << part;
                int sendcount = m.size();
                StopTimer(SERIALIZATION_TIMER);
                int total = all_sum(sendcount);
                StopTimer(COMMUNICATION_TIMER);
                //------------------------ strategy choosing END ------------------------
                if (total <= AGGSWITCH)
                    slaveGather(*part);
                else {
                    send_ibinstream(m, MASTER_RANK);
                }
                //scattering FinalT
                slaveBcast(*((FinalT*)global_agg));
            } else {
                //------------------------ strategy choosing BEGIN ------------------------
                int total = all_sum(0);
                //------------------------ strategy choosing END ------------------------
                //gathering PartialT
                if (total <= AGGSWITCH) {
                    vector<PartialT*> parts(_num_workers);
                    masterGather(parts);
                    for (int i = 0; i < _num_workers; i++) {
                        if (i != MASTER_RANK) {
                            PartialT* part = parts[i];
                            agg->stepFinal(part);
                            delete part;
                        }
                    }
                } else {
                    for (int i = 0; i < _num_workers; i++) {
                        if (i != MASTER_RANK) {
                            obinstream um = recv_obinstream(i);
                            PartialT* part;
                            um >> part;
                            agg->stepFinal(part);
                            delete part;
                        }
                    }
                }
                //scattering FinalT
                FinalT* final = agg->finishFinal();
                //cannot set "global_agg=final" since MASTER_RANK works as a slave, and agg->finishFinal() may change
                *((FinalT*)global_agg) = *final; //deep copy
                masterBcast(*((FinalT*)global_agg));
            }
        }
    }

    //user-defined graphLoader ==============================
    virtual VertexT* toVertex(char* line) = 0; //this is what user specifies!!!!!!

    void load_vertex(VertexT* v)
    { //called by load_graph
        add_vertex(v);
    }


    HashT hash;
    VertexContainer vertexes;
    int active_count;

    MessageBuffer<VertexT>* message_buffer;
    Combiner<MessageT>* combiner;
    AggregatorT* aggregator;
};

#endif

