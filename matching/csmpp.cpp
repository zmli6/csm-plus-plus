#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <time.h>
#include <assert.h>
#include <bitset>

#include "utils/types.h"
#include "utils/globals.h"
#include "utils/utils.h"
#include "graph/graph.h"
#include "matching/csmpp.h"


csmpp::csmpp(Graph& data_graph, 
        Graph& query_graph,
        std::vector<Graph> & multiQueryGraph, 
        bool print_data,
        uint max_num_results,
        bool print_prep, 
        bool print_enum, 
        bool homo,
        std::string sampleCollectPath,
        bool LR,
        uint sampleNum,
        std::string SRP)
: matching(
    multiQueryGraph, query_graph, data_graph,  max_num_results, 
    print_prep, print_enum, homo)
{
    this->print_data = print_data;
    if(this->print_data == true){
        this->edgeNum = 0;
        this->graphNum = 0;
        this->updateTime = 0;
    }
    this->sample_ = SRP;
    this->Sample.setSample(sampleNum, this->sample_, "../data/sample2SEE.txt");
    this->FileSystem.init(this->sample_);
    this->sampleCollectPath = sampleCollectPath;
    this->LRButton = LR;
    this->intersectionResult.resize(100);
    this->combineStack.resize(100);
    this->headRecord.resize(100);
    this->stackSize.resize(100);
    this->type.resize(100);
    this->stackHead.resize(100);
    if(this->sampleCollectPath != "" && this->LRButton){
        this->FileSystem.writeSampleOpen();
        this->collectBotton = true;
        
        this->sampleNum = sampleNum;
    }
    
}

void csmpp::Preprocessing()
{
    std::cout << "make data init" << std::endl;
    this->data_.indexInit();
    for(int i = 0; i < this->queryVec.size(); i++){
        this->queryVec[i].indexInit();
        this->queryVec[i].MatchOrderInit();
        this->queryVec[i].matchOrderTypeSet();
        
        //this->queryVec[i].MatchOrderAllPrint();
    }
    this->EdgeFindQueryInit();
    this->initEdgeFindQueryInit();
    if(this->LRButton){
        this->decisionMakeSystem.init(this->queryVec, this->data_, 8, this->sample_);
        std::cout << "DMS init success" << std::endl;
        this->setLRSampleUpdate();
    }
    if(this->sampleCollectPath != "" && this->LRButton){
        int i = 0;
        while (!this->updates_.empty() && this->LRButton && this->collectBotton)
        {
            InsertUnit insert = this->updates_.front();
            this->updates_.pop();
            if (insert.type == 'e' && insert.is_add)
            {
                this->data_.AddEdge(insert.id1, insert.id2, insert.label);
                this->data_.indexUpdate(insert.id1, insert.id2, insert.label, true);
                if(this->sampleNum >= 0 ){
                    this->searchInitSampleCollect(insert.id1, insert.id2, insert.label, true);
                }
            }
            else if (insert.type == 'e' && !insert.is_add)
            {
                this->RemoveEdge(insert.id1, insert.id2);
            }
            i++;
        }
        this->FileSystem.writeSampleClose();
        this->decisionMakeSystem.SetPredict(this->Sample.getPosSampleCount(), this->Sample.getNegSampleCount());
        std::cout << "DMS setting is " << this->decisionMakeSystem.predictStatus << std::endl;
    }
    std::cout << "success preprocessing" << std::endl;
}

void csmpp::setLRSampleUpdate(){
    if (!io::file_exists(this->sampleCollectPath.c_str()))
    {
        std::cout << "Failed to open: " << this->sampleCollectPath << std::endl;
        exit(-1);
    }
    std::ifstream ifs(this->sampleCollectPath);

    std::string type;
    while (ifs >> type)
    {
        uint from_id, to_id, label;
        ifs >> from_id >> to_id >> label;
        this->updates_.emplace('e', type == "e", from_id, to_id, label);
    }
    ifs.close();
    std::cout << "# LR sample update size = " << this->updates_.size() << std::endl;
}

bool csmpp::indexCheck(uint data_v, uint query_v, uint queryID){
    const uint dataGraphVLabelMaxNum = this->data_.NumVLabels();
    const uint dataGraphELabelMaxNum = this->data_.NumELabels();
    const uint queryGraphVLabelMaxNum = this->queryVec[queryID].NumVLabels();
    const uint queryGraphELabelMaxNum = this->queryVec[queryID].NumELabels();
    const auto & dataGraph = this->data_.index[data_v];
    const auto & queryGraph = this->queryVec[queryID].index[query_v];
    for(int i = 0; i < queryGraphVLabelMaxNum; i++){
        if(i < dataGraphVLabelMaxNum){
            if(dataGraph[i] < queryGraph[i]){
                return false;
            }
        }
        else{
            if(queryGraph[i] > 0){
                return false;
            }
        }
    }
    for(int i = 0; i < queryGraphELabelMaxNum; i++){
        if(i < dataGraphELabelMaxNum){
            if(dataGraph[dataGraphVLabelMaxNum + i] < queryGraph[queryGraphVLabelMaxNum + i]){
                return false;
            }
        }
        else{
            if(queryGraph[queryGraphVLabelMaxNum + i] > 0){
                return false;
            }
        }
    }
    return true;
}

void csmpp::InitialMatching()
{
    for(int i = 0; i < this->data_.NumEdges(); ++i){
        const auto & edge = this->data_.GetEdge(i);
        uint v1 = edge.GetV1();
        uint v2 = edge.GetV2();
        uint elabel = edge.GeteLabel();
        //this->searchInit(v1, v2, elabel, init);
        if(reach_time_limit){
            break;
        }
    }
}

void csmpp::AddEdge(uint v1, uint v2, uint label)
{
    this->data_.AddEdge(v1, v2, label);
    this->data_.indexUpdate(v1, v2, label, true);
    this->searchInit(v1, v2, label, pos);
}

void csmpp::RemoveEdge(uint v1, uint v2)
{
    uint label = std::get<2>(this->data_.GetEdgeLabel(v1, v2));
    this->searchInit(v1, v2, label, neg);
    this->data_.RemoveEdge(v1, v2);
    this->data_.indexUpdate(v1, v2, label, false);
}

void csmpp::AddVertex(uint id, uint label)
{
    data_.AddVertex(id, label);
    
    visited_.resize(id + 1, false);
}

void csmpp::RemoveVertex(uint id)
{
    data_.RemoveVertex(id);
}

void csmpp::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0ul;
    num_vertices = 0ul;
}

void csmpp::EdgeFindQueryInit(){
    for(uint i = 0; i < this->queryVec.size(); i++){
        Graph query = this->queryVec[i];
        for(auto & edge : query.GetEdge()){
            uint label1 = edge.GetV1Label();
            uint label2 = edge.GetV2Label();
            if(label1 > label2)  {
                std::swap(label1, label2);
            }
            Edge triggerEdge(label1, label2, edge.GeteLabel());
            this->updateEdgeFindQuery[triggerEdge].push_back(std::make_pair(i, edge.GetIndex()));
        }
    }
}

void csmpp::initEdgeFindQueryInit(){
    for(uint i = 0; i < this->queryVec.size(); i++){
        const Edge & edge = this->queryVec[i].GetEdge(0);
        uint label1 = edge.GetV1Label();
        uint label2 = edge.GetV2Label();
        if(label1 > label2){
            std::swap(label1, label2);
        }
        Edge triggerEdge(label1, label2, edge.GeteLabel());
        this->initEdgeFindQuery[triggerEdge].push_back(std::make_pair(i, 0));
    }
}

std::vector<std::pair<uint, uint>> csmpp::UpdateEdgeFineQuery(uint v1, uint v2, uint edgeLabel, searchType type){
    uint v1Label = this->data_.GetVertexLabel(v1);
    uint v2Label = this->data_.GetVertexLabel(v2);
    if(v1Label > v2Label){
        std::swap(v1Label, v2Label);
    }
    Edge edge(v1Label, v2Label, edgeLabel);
    if(type == pos || type == neg){
        if(this->updateEdgeFindQuery.find(edge) != this->updateEdgeFindQuery.end()){
            return this->updateEdgeFindQuery[edge];
        }
        else{
            return {};
        }
    }
    else{
        if(this->initEdgeFindQuery.find(edge) != this->initEdgeFindQuery.end()){
            return this->initEdgeFindQuery[edge];
        }
        else{
            return {};
        }
    }
}

void csmpp::searchInit(uint v1, uint v2, uint label, searchType type){
    const auto &  querySeries = this->UpdateEdgeFineQuery(v1, v2, label, type);

    std::vector<std::pair<uint, uint>> queryCandidate;

    for(int index = 0; index < querySeries.size(); index++){
        uint queryGraph = querySeries[index].first;
        uint queryEdge = querySeries[index].second;
        Edge _queryEdge = this->queryVec[queryGraph].GetEdge(queryEdge);
        uint v1label = this->data_.GetVertexLabel(v1);
        uint v2label = this->data_.GetVertexLabel(v2);
        uint v1Query = _queryEdge.GetV1();
        uint v2Query = _queryEdge.GetV2();
        if(v1label == v2label){
            if((this->indexCheck(v1, v1Query, queryGraph) && this->indexCheck(v2, v2Query, queryGraph))||
            (this->indexCheck(v1, v2Query, queryGraph) && this->indexCheck(v2, v1Query, queryGraph))){
                queryCandidate.push_back(querySeries[index]);
            }
        }
        else{
            if(v1label == _queryEdge.GetV1Label()){
                if(this->indexCheck(v1, v1Query, queryGraph) &&
                    this->indexCheck(v2, v2Query, queryGraph)
                ) {
                    queryCandidate.push_back(querySeries[index]);
                }
            }
            else{
                if(this->indexCheck(v2, v1Query, queryGraph) &&
                    this->indexCheck(v1, v2Query, queryGraph)
                ) {
                    queryCandidate.push_back(querySeries[index]);
                }
            }
        }
    }
    int i = 0;
    for(auto item : queryCandidate){
        auto  _edge = this->queryVec[item.first].GetEdge(item.second);
        uint _edgeV1label = _edge.GetV1Label();
        uint _edgeV2label = _edge.GetV2Label();
        uint v1query = _edge.GetV1();
        uint v2query = _edge.GetV2();
        const auto & matchOrder = this->queryVec[item.first].GetMatchOrder(item.second);
        if(_edgeV1label != _edgeV2label){
            if(this->data_.GetVertexLabel(v1) != this->queryVec[item.first].GetVertexLabel(matchOrder[0])){
                std::swap(v1, v2);
            }
            this->matchVertex(v1);
            this->matchVertex(v2);
            this->queryVec[item.first].isolatedVertexTimes.clear();
            this->searchVertex(item.first, item.second, type);

            this->popVertex(v2);
            this->popVertex(v1);
        }
        else{
            for(int i = 0; i < 2; i++){
                if(this->indexCheck(v1, matchOrder[0], item.first)){
                    this->matchVertex(v1);
                    this->matchVertex(v2);
                    this->queryVec[item.first].isolatedVertexTimes.clear();
                    this->searchVertex(item.first, item.second, type);
                    this->popVertex(v2);
                    this->popVertex(v1);
                }
                std::swap(v1,v2);// round 2 need 
            }
        }
        if(reach_time_limit) return;
        i++;
      }
}

void csmpp::searchVertex(uint queryIndex, uint edgeIndex, searchType type){
    uint depth = this->match.size();
    const auto & matchOrder = this->queryVec[queryIndex].GetMatchOrder(edgeIndex);
    const uint queryVexter = matchOrder[depth];
    uint vertexLabel = this->queryVec[queryIndex].GetVertexLabel(queryVexter);
    const auto & desItem = this->queryVec[queryIndex].GetDescList(edgeIndex, depth);//<v1Index,v1Label,eLabel>
    vertexType currentSearchVertexType = this->queryVec[queryIndex].getVertexType(edgeIndex, depth);
    const auto & freezeIndex = this->queryVec[queryIndex].getUnfreezeList(edgeIndex, depth);
    LRAndIndexCheckType dicision = this->decisionMakeSystem.makeDecision(currentSearchVertexType, desItem.size(), freezeIndex.size());
    std::vector<uint>candidate;
    if(desItem.size() != 0){
        uint min_u_index = INT_MAX;
        uint min_u_neighborSize = INT_MAX;
        for(int i = 0; i < desItem.size(); i++){
            const auto & item = desItem[i];
            uint v = this->match[std::get<0>(item)];
            uint Size = this->data_.getIndexValue(v, vertexLabel);
            if(Size < min_u_neighborSize){
                min_u_index = i;
                min_u_neighborSize = Size;
            }
            else if(Size == min_u_neighborSize){
                uint preID = this->match[std::get<0>(desItem[min_u_index])];
                if(this->data_.GetNeighbors(preID).size() > this->data_.GetNeighbors(v).size()){
                    min_u_index = i;
                    min_u_neighborSize = Size;
                }
            }
        }
        uint min_u = this->match[std::get<0>(desItem[min_u_index])];
        uint min_u_elabel = std::get<2>(desItem[min_u_index]);//elabel
        const auto& MinNeighbor = this->data_.GetNeighbors(min_u);
        const auto& q_nbr_labels = this->data_.GetNeighborLabels(min_u);
        for(int i = 0; i < MinNeighbor.size(); i++){
            const uint v = MinNeighbor[i];
            if(
                this->data_.GetVertexLabel(v) != vertexLabel ||
                q_nbr_labels[i] != min_u_elabel
            )continue;
            if(!homomorphism_ && this->visited_[v] == true) continue;
            bool joinable = true;
            for(int k = 0; k < desItem.size(); k++){
                if(k == min_u_index) continue;
                const uint data_V = this->match[std::get<0>(desItem[k])];
                const uint elabel = std::get<2>(desItem[k]);
                const auto & dataVNeighbor = this->data_.GetNeighbors(data_V);
                auto it = std::lower_bound(dataVNeighbor.begin(), dataVNeighbor.end(), v);
                uint dis = std::distance(dataVNeighbor.begin(), it);
                if(it == dataVNeighbor.end() || 
                this->data_.GetNeighborLabels(data_V)[dis] != elabel||
                *it != v
                ){
                    joinable = false;
                    break;
                }
            }
            if(!joinable)continue;
            candidate.push_back(v);
        }
        if(candidate.size() == 0){
            return;
        }
        if(dicision == Part1DoLR){
            if(this->decisionMakeSystem.predict(candidate, queryIndex, edgeIndex, depth, this->data_)){
                dicision = Part1JustCheck;
            }
        }
        if(dicision == Part1JustCheck){
            std::vector<uint> candidateCopy;
            candidateCopy.reserve(candidate.size());
            std::swap(candidate, candidateCopy);
            for(int i = 0; i < candidateCopy.size(); ++i){
                if(this->indexCheck(candidateCopy[i], queryVexter, queryIndex)){
                    candidate.push_back(candidateCopy[i]);
                }
            }
        }
        if(candidate.size() == 0){
            return;
        }
    }
    if(!freezeIndex.empty()){
        this->headRecord[depth].clear();
        this->intersectionResult[depth].clear();
        this->combineStack[depth].clear();  
        std::vector<std::vector<uint>> needToCombine;
        needToCombine.reserve(freezeIndex.size());
        for(int i = 0; i < freezeIndex.size(); i++){
            if(!(this->queryVec[queryIndex].getDestListSize(edgeIndex, freezeIndex[i]) != 0 && this->queryVec[queryIndex].getUnfreezeListSize(edgeIndex, freezeIndex[i]) != 0)){
                std::vector<uint> freezeCandidateCopy;
                freezeCandidateCopy.reserve(this->matchCandidate[freezeIndex[i]].size());
                std::swap(freezeCandidateCopy, this->matchCandidate[freezeIndex[i]]);
                for(int k = 0; k < freezeCandidateCopy.size(); ++k){
                    if(this->indexCheck(freezeCandidateCopy[k], matchOrder[freezeIndex[i]], queryIndex)){
                        this->matchCandidate[freezeIndex[i]].push_back(freezeCandidateCopy[k]);
                    }
                }
                if(this->matchCandidate[freezeIndex[i]].size() == 0){
                    return;//all kill
                }
            }
            needToCombine.push_back(this->matchCandidate[freezeIndex[i]]);
            this->headRecord[depth].push_back(this->matchCandidate[freezeIndex[i]].size());
        }
        this->stackSize[depth] = needToCombine.size();
        this->stackHead[depth] = 0;
        this->combineStack[depth].resize(this->stackSize[depth]);
        this->intersectionResult[depth].resize(this->stackSize[depth]);
        this->type[depth] = runningStack;
        while(this->headRecord[depth][0] >= 0){
            while(this->stackHead[depth] < this->stackSize[depth]){
                if(this->headRecord[depth][this->stackHead[depth]] == 0){
                    bool pending = this->headChange(needToCombine, depth);
                    if(pending == false){
                        return;
                    }
                }
                int replaceIndex = this->stackHead[depth];// need to add
                uint currentVertex = needToCombine[replaceIndex][this->headRecord[depth][replaceIndex] - 1];
                if(vertexPushCheck(currentVertex, vertexLabel, candidate, depth, queryIndex, queryVexter) == false){
                    this->headRecord[depth][replaceIndex]--;
                }
                else{
                    this->combinePushBack(currentVertex, replaceIndex, depth);
                }
            }
            this->queryVec[queryIndex].setVertexStatus(edgeIndex, freezeIndex, freeVertex);
            this->setUnVisitedPatch(this->combineStack[depth]);
            this->setMatchVertex(freezeIndex, this->combineStack[depth]);
            if(dicision == Part2DoLR){
                if(this->decisionMakeSystem.predict(this->getItersectionTop(depth), queryIndex, edgeIndex, depth, this->data_)){
                    dicision == Part2JustCheck;
                }
            }
            if(dicision == Part2JustCheck){
                std::vector<uint> Part2Candidate;
                Part2Candidate.reserve(this->getItersectionTop(depth).size());
                for(auto dataV : this->getItersectionTop(depth)){
                    if(this->indexCheck(dataV, queryVexter, queryIndex)){
                        Part2Candidate.push_back(dataV);
                    }
                }
                if(Part2Candidate.size() != 0){
                    if(currentSearchVertexType == freeVertex){
                        for(auto dataV : Part2Candidate){
                            this->matchVertex(dataV);
                            this->searchVertex(queryIndex, edgeIndex, type);
                            this->popVertex(dataV);
                        }
                    }
                    else{         
                        this->matchVertex(Part2Candidate);
                        if(currentSearchVertexType == isolatedVertex){
                            this->queryVec[queryIndex].isolatedVertexTimesAdd(Part2Candidate);
                        }
                        if(depth == this->queryVec[queryIndex].NumVertices() - 1){
                            this->addMatchResult(queryIndex, edgeIndex, type);
                        }
                        else{
                            this->searchVertex(queryIndex, edgeIndex, type);
                        }
                        if(currentSearchVertexType == isolatedVertex){
                            this->queryVec[queryIndex].isolatedVertexTimesMinus(Part2Candidate);
                        }
                        this->popVertex();
                    }
                }
            }
            else{
                if(currentSearchVertexType == freeVertex){
                    for(auto dataV : this->getItersectionTop(depth)){
                        this->matchVertex(dataV);
                        this->searchVertex(queryIndex, edgeIndex, type);
                        this->popVertex(dataV);
                    }
                }
                else{         
                    this->matchVertex(this->getItersectionTop(depth));
                    if(currentSearchVertexType == isolatedVertex){
                        this->queryVec[queryIndex].isolatedVertexTimesAdd(this->getItersectionTop(depth));
                    }
                    if(depth == this->queryVec[queryIndex].NumVertices() - 1){
                        this->addMatchResult(queryIndex, edgeIndex, type);
                    }
                    else{
                        this->searchVertex(queryIndex, edgeIndex, type);
                    }
                    if(currentSearchVertexType == isolatedVertex){
                        this->queryVec[queryIndex].isolatedVertexTimesMinus(this->getItersectionTop(depth));
                    }
                    this->popVertex();
                }
            }
            this->unsetMatchVertex(freezeIndex);
            this->setVisitedPatch(this->combineStack[depth]);
            this->queryVec[queryIndex].setVertexStatus(edgeIndex, freezeIndex, freezeVertex);
            this->combineStackPopTail(depth);
        }
    }
    else{
        if(currentSearchVertexType == freeVertex){
            for(auto dataV : candidate){
                this->matchVertex(dataV);
                this->searchVertex(queryIndex, edgeIndex, type);
                this->popVertex(dataV);
            }
        }
        else{
            this->matchVertex(candidate);
            if(currentSearchVertexType == isolatedVertex){
                this->queryVec[queryIndex].isolatedVertexTimesAdd(candidate);
            }
            if(depth == this->queryVec[queryIndex].NumVertices() - 1){
                this->addMatchResult(queryIndex, edgeIndex, type);
            }
            else{
                this->searchVertex(queryIndex, edgeIndex, type);
            }
            if(currentSearchVertexType == isolatedVertex){
                this->queryVec[queryIndex].isolatedVertexTimesMinus(candidate);
            }
            this->popVertex();
        }
    }
}

void csmpp::setMatchVertex(const std::vector<uint> & matchingIndex, const std::vector<int>& vertexs){
    for(int i = 0; i < matchingIndex.size(); i++){
        this->match[matchingIndex[i]] = vertexs[i];
        this->visited_[vertexs[i]] = true;
    }
}

void csmpp::unsetMatchVertex(const std::vector<uint> & matchingIndex){
    for(int i = 0; i < matchingIndex.size(); i++){
        this->visited_[this->match[matchingIndex[i]]] = false;
        this->match[matchingIndex[i]] = -1;
    }
}

void csmpp::addMatchResult(uint queryIndex, uint edgeIndex, searchType type){
    const auto & isolateVertexIndex = this->queryVec[queryIndex].getIsolatedVertexIndex(edgeIndex);
    if(print_enumeration_results_){
        std::vector<std::vector<uint>> needToCombine;
        needToCombine.reserve(isolateVertexIndex.size());
        for(int i = 0; i < isolateVertexIndex.size(); i++){
            needToCombine.push_back(this->matchCandidate[isolateVertexIndex[i]]);
        }
        uint depth = this->queryVec[queryIndex].NumVertices();
        this->headRecord[depth].clear();
        this->intersectionResult[depth].clear();
        this->combineStack[depth].clear();
        for(int i = 0; i < needToCombine.size(); i++){
            headRecord[depth].push_back(needToCombine[i].size());
        }
        this->stackSize[depth] = needToCombine.size();
        this->stackHead[depth] = 0;
        this->combineStack[depth].resize(this->stackSize[depth]);
        this->intersectionResult[depth].resize(this->stackSize[depth]);
        this->type[depth] = finalStack;
        while(this->headRecord[depth][0] >= 0){
            while(this->stackHead[depth] < this->stackSize[depth]){
                if(this->headRecord[depth][this->stackHead[depth]] == 0){
                    bool pending = this->headChange(needToCombine, depth);
                    if(pending == false){
                        return;
                    }
                }
                int replaceIndex = this->stackHead[depth];// need to add
                uint currentVertex = needToCombine[replaceIndex][this->headRecord[depth][replaceIndex] - 1];
                if(this->visited_[currentVertex] == true){
                    this->headRecord[depth][replaceIndex]--;
                }
                else{
                    this->combinePushBack(currentVertex, replaceIndex, depth);
                }
            }
            if(type == pos){
                num_positive_results_ ++;
            }
            else{
                if(type == neg){
                    num_negative_results_++;
                }
                else{
                    this->num_initial_results_++;
                }
            }
            this->setUnVisitedPatch(this->combineStack[depth]);
            this->setMatchVertex(isolateVertexIndex, this->combineStack[depth]);
            this->printMatch();
            this->unsetMatchVertex(isolateVertexIndex);
            this->setVisitedPatch(this->combineStack[depth]);
            this->combineStackPopTail(depth);
        }
    }
    else{
        auto & isolatedVertexMap = this->queryVec[queryIndex].isolatedVertexTimes;
        std::vector<std::vector<int>> needToCombineV1;
        std::vector<uint> NoOverLeafWeight;
        bool allSame = true;
        for(int i = 0; i < isolateVertexIndex.size(); i++){
            const auto & I_isolateVertexCandidate = this->matchCandidate[isolateVertexIndex[i]];
            std::vector<int> I_needToCombine;
            I_needToCombine.reserve(I_isolateVertexCandidate.size());
            for(int k = 0; k < I_isolateVertexCandidate.size(); k++){
                int vertex = I_isolateVertexCandidate[k];
                if(this->visited_[vertex] == true){
                    continue;
                }
                I_needToCombine.push_back(vertex);
            }
            if(I_needToCombine.empty()){
                return;
            }
            needToCombineV1.push_back(I_needToCombine);
            if(allSame && i >= 1){
                if(needToCombineV1[i - 1].size() != needToCombineV1[i].size()){
                    allSame = false;
                }
            }
        }
        if(needToCombineV1.size() == 1){
            if(type == pos){
                num_positive_results_ += needToCombineV1[0].size();
            }
            else{
                if(type == neg){
                    num_negative_results_ += needToCombineV1[0].size();
                }
                else{
                    num_initial_results_ += needToCombineV1[0].size();
                }
            }
            
            return;
        }
        if(allSame == true){
            const auto & firstItem = needToCombineV1[0];
            for(int k = 1; k < needToCombineV1.size(); k++){
                const auto & kSelf = needToCombineV1[k];
                if(firstItem != kSelf){
                    allSame = false;
                }
            }
            if(allSame == true){
                uint Matchresult = 1;
                for(int i = 0; i < needToCombineV1.size(); i++){
                    Matchresult *= (firstItem.size() - i);
                }
                if(type == pos){
                    num_positive_results_ += Matchresult;
                }
                else{
                    if(type == neg){
                        num_negative_results_ += Matchresult;
                    }
                    else{
                        num_initial_results_ += Matchresult;
                    }
                }
                
                return;
            }
        }
        std::vector<std::vector<int>> needToCombine;
        for(int i = 0; i < needToCombineV1.size(); i++){
            const auto & I_isolateVertexCandidate = needToCombineV1[i];
            std::vector<int> I_needToCombine;
            I_needToCombine.reserve(I_isolateVertexCandidate.size());
            int I_NoOverLeafWeight = 0;
            for(int k = 0; k < I_isolateVertexCandidate.size(); k++){
                int vertex = I_isolateVertexCandidate[k];
                if(isolatedVertexMap[vertex] > 1){
                    I_needToCombine.push_back(vertex);
                }
                else{
                    I_NoOverLeafWeight++;
                }
            }
            NoOverLeafWeight.push_back(I_NoOverLeafWeight);
            if(I_NoOverLeafWeight > 0){
                I_needToCombine.push_back(-1);
            }
            needToCombine.push_back(I_needToCombine);
        }
        uint depth = this->queryVec[queryIndex].NumVertices();
        this->headRecord[depth].clear();
        this->intersectionResult[depth].clear();
        this->combineStack[depth].clear();
        for(int i = 0; i < needToCombine.size(); i++){
            headRecord[depth].push_back(needToCombine[i].size());
        }
        this->stackSize[depth] = needToCombine.size();
        this->stackHead[depth] = 0;
        this->combineStack[depth].resize(this->stackSize[depth]);
        this->intersectionResult[depth].resize(this->stackSize[depth]);
        this->type[depth] = finalStack;
        uint totalMatch = 1;
        int currentVertex;
        int replaceIndex;
        while(this->headRecord[depth][0] >= 0){
            while(this->stackHead[depth] < this->stackSize[depth]){
                if(this->headRecord[depth][this->stackHead[depth]] == 0){
                    bool pending = this->headChange(needToCombine, depth, totalMatch, NoOverLeafWeight);
                    if(pending == false){
                        return;
                    }
                }
                replaceIndex = this->stackHead[depth];// need to add
                currentVertex = needToCombine[replaceIndex][this->headRecord[depth][replaceIndex] - 1];
                if(currentVertex != -1 && this->visited_[currentVertex] == true){
                    this->headRecord[depth][replaceIndex]--;
                }
                else{
                    this->combinePushBack(currentVertex, replaceIndex, depth, totalMatch, NoOverLeafWeight);
                }
            }
            if(type == pos){
                num_positive_results_ += totalMatch;
            }
            else{
                if(type == neg){
                    num_negative_results_ += totalMatch;
                }
                else{
                    num_initial_results_ += totalMatch;
                }
                
            }
            
            this->combineStackPopTail(depth, totalMatch, NoOverLeafWeight);
        }
    }
}

void csmpp::printMatch(){
    std::vector<int>matchCopy(match);
    std::cout << "vertex list" << std::endl;
    for(int i = 0; i < this->match.size(); i++){
        std::cout << this->match[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "sort vertex list" << std::endl;
    std::sort(matchCopy.begin(), matchCopy.end());
    for(int i = 0; i < matchCopy.size(); i++){
        std::cout << matchCopy[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "vertex Labels " << std::endl;
    for(int i = 0; i < this->match.size(); i++){
        std::cout << this->data_.GetVertexLabel(match[i])<< " ";
    }
    std::cout << std::endl;
}

void csmpp::matchVertex(uint vertex){
    this->match.push_back(vertex);
    this->matchCandidate.push_back({});
    this->visited_[vertex] = true;
}

void csmpp::matchVertex(const std::vector<uint> & Candidate){
    this->match.push_back(-1);
    this->matchCandidate.push_back(Candidate);
}

void csmpp::popVertex(uint vertex){
    this->match.pop_back();
    this->matchCandidate.pop_back();
    this->visited_[vertex] = false;
}

void csmpp::popVertex(){
    this->match.pop_back();
    this->matchCandidate.pop_back();
}

bool csmpp::vertexPushCheck(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection, int depth, uint queryIndex, uint queryVertex){
    if(this->visited_[vertex] == true){
        return false;
    }
    if(this->type[depth] == runningStack){
        if(this->getIntersection(vertex, queryVertexLabel, needToIntersection, depth, queryIndex, queryVertex) == false){
            return false;
        }
    }
    return true;
}

bool csmpp::getIntersection(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection, int depth, uint queryIndex, uint queryVertex){
    std::vector<uint>result;
    result.reserve(100000);
    const auto & vNeighbors = this->data_.GetNeighbors(vertex);
    if(!needToIntersection.empty()){
        std::set_intersection(vNeighbors.begin(), vNeighbors.end(), needToIntersection.begin(), needToIntersection.end(), std::back_inserter(result));
        if(result.empty()){
            return false;
        }
    }
    else{
        for(int i = 0; i < vNeighbors.size(); i++){
            if(this->visited_[vNeighbors[i]] == false && this->data_.GetVertexLabel(vNeighbors[i]) == queryVertexLabel){
                result.push_back(vNeighbors[i]);
            }
        }
    }
    if(this->stackHead[depth] == 0){
        this->intersectionResult[depth][this->stackHead[depth]] = result;
        return true;
    }
    std::vector<uint> finalResult;
    finalResult.reserve(100000);
    const std::vector<uint> & preVec = this->getItersectionTop(depth);
    std::set_intersection(preVec.begin(), preVec.end(), result.begin(), result.end(), std::back_inserter(finalResult));
    if(finalResult.empty()){
        return false;
    }
    this->intersectionResult[depth][this->stackHead[depth]] = finalResult;
    return true;
}

const std::vector<uint> & csmpp::getItersectionTop(int depth) const{
    return this->intersectionResult[depth][this->stackHead[depth] - 1];
}

void csmpp::combineStackPopTail(int depth){
    uint vertex = this->combineStack[depth][this->stackHead[depth] - 1];
    this->stackHead[depth]--;
    this->visited_[vertex] = false;
}

void csmpp::combineStackPopTail(int depth, uint & totalMatch, const std::vector<uint>& NoOverLeftWeight){
    uint vertex = this->combineStack[depth][this->stackHead[depth] - 1];
    
    if(vertex != -1){
        this->visited_[vertex] = false;
    }
    else{
        totalMatch /= NoOverLeftWeight[this->stackHead[depth] - 1];
    }
    this->stackHead[depth]--;
}

void csmpp::combinePushBack(uint vertex, uint vertexInWhichCandidateNeighbor, int depth){
    this->visited_[vertex] = true;
    this->combineStack[depth][this->stackHead[depth]] = vertex;
    this->stackHead[depth]++;
    this->headRecord[depth][vertexInWhichCandidateNeighbor]--;
}

void csmpp::combinePushBack(int vertex, uint vertexInWhichCandidateNeighbor, int depth, uint & totalMatch, const std::vector<uint>& NoOverLeftWeight){
    if(vertex != -1){
        this->visited_[vertex] = true;
    }
    else{
        totalMatch *= NoOverLeftWeight[vertexInWhichCandidateNeighbor];
    }
    this->combineStack[depth][this->stackHead[depth]] = vertex;
    this->stackHead[depth]++;
    this->headRecord[depth][vertexInWhichCandidateNeighbor]--;
}

bool csmpp::headChange(const std::vector<std::vector<uint>> & needToCombine, int depth){
    int iHead = this->stackHead[depth];
    while(this->headRecord[depth][iHead] <= 0){
        if(iHead <= 0){
            return false;
        }
        this->headRecord[depth][iHead] = needToCombine[iHead].size();
        iHead--;
        this->combineStackPopTail(depth);
    }
    return true;
}

bool csmpp::headChange(const std::vector<std::vector<int>> & needToCombine, int depth, uint & totalMatch, const std::vector<uint>& NoOverLeftWeight){
    int iHead = this->stackHead[depth];
    while(this->headRecord[depth][iHead] <= 0){
        if(iHead <= 0){
            return false;
        }
        this->headRecord[depth][iHead] = needToCombine[iHead].size();
        iHead--;
        this->combineStackPopTail(depth, totalMatch, NoOverLeftWeight);
    }
    return true;
}

void csmpp::setVisitedPatch(const std::vector<int> & vertex){
    for(int i = 0; i < vertex.size(); i++){
        this->visited_[vertex[i]] = true;
    }
}

void csmpp::setUnVisitedPatch(const std::vector<int> & vertex){
    for(int i = 0; i < vertex.size(); i++){
        this->visited_[vertex[i]] = false;
    }
}

void csmpp::searchVertexSampleCollect(uint queryIndex, uint edgeIndex, bool recordButton, int collectDepth){
    uint depth = this->match.size();
    vertexType currentSearchVertexType = this->queryVec[queryIndex].getVertexType(edgeIndex, depth);
    const auto & matchOrder = this->queryVec[queryIndex].GetMatchOrder(edgeIndex);
    const uint queryVexter = matchOrder[depth];
    uint vertexLabel = this->queryVec[queryIndex].GetVertexLabel(queryVexter);
    const auto & desItem = this->queryVec[queryIndex].GetDescList(edgeIndex, depth);//<v1Index,v1Label,eLabel>
    const auto & freezeIndex = this->queryVec[queryIndex].getUnfreezeList(edgeIndex, depth);
    LRAndIndexCheckType decision = this->decisionMakeSystem.makeDecision(currentSearchVertexType, this->queryVec[queryIndex].getDestListSize(edgeIndex, depth), this->queryVec[queryIndex].getUnfreezeListSize(edgeIndex, depth));
    if(collectDepth == depth){
        recordButton = true;
    }
    std::vector<uint>candidate;
    if(desItem.size() != 0){
        uint min_u_index = INT_MAX;
        uint min_u_neighborSize = INT_MAX;
        for(int i = 0; i < desItem.size(); i++){
            const auto & item = desItem[i];
            uint v = this->match[std::get<0>(item)];
            uint Size = this->data_.getIndexValue(v, vertexLabel);
            if(Size < min_u_neighborSize){
                min_u_index = i;
                min_u_neighborSize = Size;
            }
            else if(Size == min_u_neighborSize){
                uint preID = this->match[std::get<0>(desItem[min_u_index])];
                if(this->data_.GetNeighbors(preID).size() > this->data_.GetNeighbors(v).size()){
                    min_u_index = i;
                    min_u_neighborSize = Size;
                }
            }
        }
        uint min_u = this->match[std::get<0>(desItem[min_u_index])];
        uint min_u_elabel = std::get<2>(desItem[min_u_index]);//elabel
        const auto& MinNeighbor = this->data_.GetNeighbors(min_u);
        const auto& q_nbr_labels = this->data_.GetNeighborLabels(min_u);
        for(int i = 0; i < MinNeighbor.size(); i++){
            const uint v = MinNeighbor[i];
            if(
                this->data_.GetVertexLabel(v) != vertexLabel ||
                q_nbr_labels[i] != min_u_elabel
            )continue;
            if(!homomorphism_ && this->visited_[v] == true) continue;
            bool joinable = true;
            for(int k = 0; k < desItem.size(); k++){
                if(k == min_u_index) continue;
                const uint data_V = this->match[std::get<0>(desItem[k])];
                const uint elabel = std::get<2>(desItem[k]);
                const auto & dataVNeighbor = this->data_.GetNeighbors(data_V);
                auto it = std::lower_bound(dataVNeighbor.begin(), dataVNeighbor.end(), v);
                uint dis = std::distance(dataVNeighbor.begin(), it);
                if(it == dataVNeighbor.end() || 
                this->data_.GetNeighborLabels(data_V)[dis] != elabel||
                *it != v
                ){
                    joinable = false;
                    break;
                }
            }
            if(!joinable)continue;
            if(decision == Part1JustCheck){ 
                if(this->indexCheckForSampleCollect(v, queryVexter, queryIndex) == false)continue;
            }
            candidate.push_back(v);
        }
        if(candidate.size() == 0){
            return;
        }
    }
    if(recordButton == true && decision == Part1DoLR){
        recordButton = false;
        sampleUnit Unit;
        {
            Unit.descNeighborMinus1 = desItem.size() - 1;
            Unit.minSize = candidate.size();
            double DegreeCount = 0;
            for(int i = 0; i < candidate.size(); i++){
                DegreeCount += this->data_.GetDegree(candidate[i]);
            }
            Unit.avgDegreeSize = DegreeCount / candidate.size();
            Unit.motifSize = this->queryVec[queryIndex].NumELabels() + this->queryVec[queryIndex].NumVLabels();
            Unit.rank = this->data_.getFrequencyRank(true, vertexLabel);//true for vertex; false for edge
            uint queryGraphSize = this->queryVec[queryIndex].NumVertices();
            Unit.sizeMinusDepth = queryGraphSize - (depth + 1);
            Unit.sizeMinusDepthDivisionSize = (Unit.sizeMinusDepth + 0.0) / queryGraphSize;
            Unit.freezeVertexAfter = this->queryVec[queryIndex].getFreezeVertexNumAfter(edgeIndex, depth);
            Unit.isolatedVertexAfter = this->queryVec[queryIndex].getIsolatedVertexNumAfter(edgeIndex, depth);
        }
        
        Unit.nomotifTimer.StartTimer();
        if(!freezeIndex.empty()){
            this->headRecord[depth].clear();
            this->intersectionResult[depth].clear();
            this->combineStack[depth].clear();  
            std::vector<std::vector<uint>> needToCombine;
            needToCombine.reserve(freezeIndex.size());
            for(int i = 0; i < freezeIndex.size(); i++){
                needToCombine.push_back(this->matchCandidate[freezeIndex[i]]);
                this->headRecord[depth].push_back(this->matchCandidate[freezeIndex[i]].size());
            }
            this->stackSize[depth] = needToCombine.size();
            this->stackHead[depth] = 0;
            this->combineStack[depth].resize(this->stackSize[depth]);
            this->intersectionResult[depth].resize(this->stackSize[depth]);
            this->type[depth] = runningStack;
            while(this->headRecord[depth][0] >= 0){
                while(this->stackHead[depth] < this->stackSize[depth]){
                    if(this->headRecord[depth][this->stackHead[depth]] == 0){
                        bool pending = this->headChange(needToCombine, depth);
                        if(pending == false){
                            return;
                        }
                    }
                    int replaceIndex = this->stackHead[depth];// need to add
                    uint currentVertex = needToCombine[replaceIndex][this->headRecord[depth][replaceIndex] - 1];
                    if(vertexPushCheck(currentVertex, vertexLabel, candidate, depth, queryIndex, queryVexter) == false){
                        this->headRecord[depth][replaceIndex]--;
                    }
                    else{
                        this->combinePushBack(currentVertex, replaceIndex, depth);
                    }
                }
                this->queryVec[queryIndex].setVertexStatus(edgeIndex, freezeIndex, freeVertex);
                this->setUnVisitedPatch(this->combineStack[depth]);
                this->setMatchVertex(freezeIndex, this->combineStack[depth]);
                if(currentSearchVertexType == freeVertex){
                    for(auto dataV : this->getItersectionTop(depth)){
                        this->matchVertex(dataV);
                        //this->searchVertex(queryIndex, edgeIndex);
                        this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                        this->popVertex(dataV);
                    }
                }
                else{         
                    this->matchVertex(this->getItersectionTop(depth));
                    if(currentSearchVertexType == isolatedVertex){
                        this->queryVec[queryIndex].isolatedVertexTimesAdd(this->getItersectionTop(depth));
                    }
                    if(depth == this->queryVec[queryIndex].NumVertices() - 1){
                        this->addMatchResultSampleCollect(queryIndex, edgeIndex);
                    }
                    else{
                        this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                    }
                    if(currentSearchVertexType == isolatedVertex){
                        this->queryVec[queryIndex].isolatedVertexTimesMinus(this->getItersectionTop(depth));
                    }
                    this->popVertex();
                }
                this->unsetMatchVertex(freezeIndex);
                this->setVisitedPatch(this->combineStack[depth]);
                this->queryVec[queryIndex].setVertexStatus(edgeIndex, freezeIndex, freezeVertex);
                this->combineStackPopTail(depth);
            }
        }
        else{
            if(currentSearchVertexType == freeVertex){
                for(auto dataV : candidate){
                    this->matchVertex(dataV);
                    this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                    this->popVertex(dataV);
                }
            }
            else{
                this->matchVertex(candidate);
                if(currentSearchVertexType == isolatedVertex){
                    this->queryVec[queryIndex].isolatedVertexTimesAdd(candidate);
                }
                if(depth == this->queryVec[queryIndex].NumVertices() - 1){
                    this->addMatchResultSampleCollect(queryIndex, edgeIndex);
                }
                else{
                    this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                }
                if(currentSearchVertexType == isolatedVertex){
                    this->queryVec[queryIndex].isolatedVertexTimesMinus(candidate);
                }
                this->popVertex();
            }
        }
        
        Unit.nomotifTimer.StopTimer();
        std::vector<uint> candidateCopy(candidate);
        candidate.clear();
        candidate.reserve(candidateCopy.size());
        Unit.motifTimer.StartTimer();
        for(int i = 0; i < candidateCopy.size(); ++i){
            if(this->indexCheckForSampleCollect(candidateCopy[i], queryVexter, queryIndex)){
                candidate.push_back(candidateCopy[i]);
            }
        }
        if(!freezeIndex.empty()){
            this->headRecord[depth].clear();
            this->intersectionResult[depth].clear();
            this->combineStack[depth].clear();  
            std::vector<std::vector<uint>> needToCombine;
            needToCombine.reserve(freezeIndex.size());
            for(int i = 0; i < freezeIndex.size(); i++){
                needToCombine.push_back(this->matchCandidate[freezeIndex[i]]);
                this->headRecord[depth].push_back(this->matchCandidate[freezeIndex[i]].size());
            }
            this->stackSize[depth] = needToCombine.size();
            this->stackHead[depth] = 0;
            this->combineStack[depth].resize(this->stackSize[depth]);
            this->intersectionResult[depth].resize(this->stackSize[depth]);
            this->type[depth] = runningStack;
            while(this->headRecord[depth][0] >= 0){
                while(this->stackHead[depth] < this->stackSize[depth]){
                    if(this->headRecord[depth][this->stackHead[depth]] == 0){
                        bool pending = this->headChange(needToCombine, depth);
                        if(pending == false){
                            return;
                        }
                    }
                    int replaceIndex = this->stackHead[depth];// need to add
                    uint currentVertex = needToCombine[replaceIndex][this->headRecord[depth][replaceIndex] - 1];
                    if(vertexPushCheck(currentVertex, vertexLabel, candidate, depth, queryIndex, queryVexter) == false){
                        this->headRecord[depth][replaceIndex]--;
                    }
                    else{
                        this->combinePushBack(currentVertex, replaceIndex, depth);
                    }
                }
                this->queryVec[queryIndex].setVertexStatus(edgeIndex, freezeIndex, freeVertex);
                this->setUnVisitedPatch(this->combineStack[depth]);
                this->setMatchVertex(freezeIndex, this->combineStack[depth]);
                if(currentSearchVertexType == freeVertex){
                    for(auto dataV : this->getItersectionTop(depth)){
                        this->matchVertex(dataV);
                        this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                        this->popVertex(dataV);
                    }
                }
                else{         
                    this->matchVertex(this->getItersectionTop(depth));
                    if(currentSearchVertexType == isolatedVertex){
                        this->queryVec[queryIndex].isolatedVertexTimesAdd(this->getItersectionTop(depth));
                    }
                    if(depth == this->queryVec[queryIndex].NumVertices() - 1){
                        this->addMatchResultSampleCollect(queryIndex, edgeIndex);
                    }
                    else{
                        this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                    }
                    if(currentSearchVertexType == isolatedVertex){
                        this->queryVec[queryIndex].isolatedVertexTimesMinus(this->getItersectionTop(depth));
                    }
                    this->popVertex();
                }
                this->unsetMatchVertex(freezeIndex);
                this->setVisitedPatch(this->combineStack[depth]);
                this->queryVec[queryIndex].setVertexStatus(edgeIndex, freezeIndex, freezeVertex);
                this->combineStackPopTail(depth);
            }
        }
        else{
            if(currentSearchVertexType == freeVertex){
                for(auto dataV : candidate){
                    this->matchVertex(dataV);
                    this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                    this->popVertex(dataV);
                }
            }
            else{
                this->matchVertex(candidate);
                if(currentSearchVertexType == isolatedVertex){
                    this->queryVec[queryIndex].isolatedVertexTimesAdd(candidate);
                }
                if(depth == this->queryVec[queryIndex].NumVertices() - 1){
                    this->addMatchResultSampleCollect(queryIndex, edgeIndex);
                }
                else{
                    this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                }
                if(currentSearchVertexType == isolatedVertex){
                    this->queryVec[queryIndex].isolatedVertexTimesMinus(candidate);
                }
                this->popVertex();
            }
        }
        Unit.motifTimer.StopTimer();
        this->FileSystem.writeSampleToFile(this->Sample.sampleUnitToStr(Unit));
        this->sampleNum--;
    }
    else{
        if(!freezeIndex.empty()){
            this->headRecord[depth].clear();
            this->intersectionResult[depth].clear();
            this->combineStack[depth].clear();  
            std::vector<std::vector<uint>> needToCombine;
            needToCombine.reserve(freezeIndex.size());
            for(int i = 0; i < freezeIndex.size(); i++){
                needToCombine.push_back(this->matchCandidate[freezeIndex[i]]);
                this->headRecord[depth].push_back(this->matchCandidate[freezeIndex[i]].size());
            }
            this->stackSize[depth] = needToCombine.size();
            this->stackHead[depth] = 0;
            this->combineStack[depth].resize(this->stackSize[depth]);
            this->intersectionResult[depth].resize(this->stackSize[depth]);
            this->type[depth] = runningStack;
            while(this->headRecord[depth][0] >= 0){
                while(this->stackHead[depth] < this->stackSize[depth]){
                    if(this->headRecord[depth][this->stackHead[depth]] == 0){
                        bool pending = this->headChange(needToCombine, depth);
                        if(pending == false){
                            return;
                        }
                    }
                    int replaceIndex = this->stackHead[depth];
                    uint currentVertex = needToCombine[replaceIndex][this->headRecord[depth][replaceIndex] - 1];
                    if(vertexPushCheck(currentVertex, vertexLabel, candidate, depth, queryIndex, queryVexter) == false){
                        this->headRecord[depth][replaceIndex]--;
                    }
                    else{
                        this->combinePushBack(currentVertex, replaceIndex, depth);
                    }
                }
                this->queryVec[queryIndex].setVertexStatus(edgeIndex, freezeIndex, freeVertex);
                this->setUnVisitedPatch(this->combineStack[depth]);
                this->setMatchVertex(freezeIndex, this->combineStack[depth]);
                if(currentSearchVertexType == freeVertex){
                    for(auto dataV : this->getItersectionTop(depth)){
                        this->matchVertex(dataV);
                        this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                        this->popVertex(dataV);
                    }
                }
                else{         
                    this->matchVertex(this->getItersectionTop(depth));
                    if(currentSearchVertexType == isolatedVertex){
                        this->queryVec[queryIndex].isolatedVertexTimesAdd(this->getItersectionTop(depth));
                    }
                    if(depth == this->queryVec[queryIndex].NumVertices() - 1){
                        this->addMatchResultSampleCollect(queryIndex, edgeIndex);
                    }
                    else{
                        this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                    }
                    if(currentSearchVertexType == isolatedVertex){
                        this->queryVec[queryIndex].isolatedVertexTimesMinus(this->getItersectionTop(depth));
                    }
                    this->popVertex();
                }
                this->unsetMatchVertex(freezeIndex);
                this->setVisitedPatch(this->combineStack[depth]);
                this->queryVec[queryIndex].setVertexStatus(edgeIndex, freezeIndex, freezeVertex);
                this->combineStackPopTail(depth);
            } 
        }
        else{
            if(currentSearchVertexType == freeVertex){
                for(auto dataV : candidate){
                    this->matchVertex(dataV);
                    this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                    this->popVertex(dataV);
                }
            }
            else{
                this->matchVertex(candidate);
                if(currentSearchVertexType == isolatedVertex){
                    this->queryVec[queryIndex].isolatedVertexTimesAdd(candidate);
                }
                if(depth == this->queryVec[queryIndex].NumVertices() - 1){
                    this->addMatchResultSampleCollect(queryIndex, edgeIndex);
                }
                else{
                    this->searchVertexSampleCollect(queryIndex, edgeIndex, recordButton, collectDepth);
                }
                if(currentSearchVertexType == isolatedVertex){
                    this->queryVec[queryIndex].isolatedVertexTimesMinus(candidate);
                }
                this->popVertex();
            }
        }
    }
}

void csmpp::addMatchResultSampleCollect(uint queryIndex, uint edgeIndex){
    const auto & isolateVertexIndex = this->queryVec[queryIndex].getIsolatedVertexIndex(edgeIndex);
    if(print_enumeration_results_){
        std::vector<std::vector<uint>> needToCombine;
        needToCombine.reserve(isolateVertexIndex.size());
        for(int i = 0; i < isolateVertexIndex.size(); i++){
            needToCombine.push_back(this->matchCandidate[isolateVertexIndex[i]]);
        }
        uint depth = this->queryVec[queryIndex].NumVertices();
        this->headRecord[depth].clear();
        this->intersectionResult[depth].clear();
        this->combineStack[depth].clear();
        for(int i = 0; i < needToCombine.size(); i++){
            headRecord[depth].push_back(needToCombine[i].size());
        }
        this->stackSize[depth] = needToCombine.size();
        this->stackHead[depth] = 0;
        this->combineStack[depth].resize(this->stackSize[depth]);
        this->intersectionResult[depth].resize(this->stackSize[depth]);
        this->type[depth] = finalStack;
        while(this->headRecord[depth][0] >= 0){
            //1.push back stack
            while(this->stackHead[depth] < this->stackSize[depth]){
                if(this->headRecord[depth][this->stackHead[depth]] == 0){
                    bool pending = this->headChange(needToCombine, depth);
                    if(pending == false){
                        return;
                    }
                }
                int replaceIndex = this->stackHead[depth];// need to add
                uint currentVertex = needToCombine[replaceIndex][this->headRecord[depth][replaceIndex] - 1];
                if(this->visited_[currentVertex] == true){
                    this->headRecord[depth][replaceIndex]--;
                }
                else{
                    this->combinePushBack(currentVertex, replaceIndex, depth);
                }
            }
            this->combineStackPopTail(depth);
        }
    }
    else{
        auto & isolatedVertexMap = this->queryVec[queryIndex].isolatedVertexTimes;
        std::vector<std::vector<int>> needToCombineV1;
        std::vector<uint> NoOverLeafWeight;
        bool allSame = true;
        for(int i = 0; i < isolateVertexIndex.size(); i++){
            const auto & I_isolateVertexCandidate = this->matchCandidate[isolateVertexIndex[i]];
            std::vector<int> I_needToCombine;
            I_needToCombine.reserve(I_isolateVertexCandidate.size());
            for(int k = 0; k < I_isolateVertexCandidate.size(); k++){
                int vertex = I_isolateVertexCandidate[k];
                if(this->visited_[vertex] == true){
                    continue;
                }
                I_needToCombine.push_back(vertex);
            }
            if(I_needToCombine.empty()){
                return;
            }
            needToCombineV1.push_back(I_needToCombine);
            if(allSame && i >= 1){
                if(needToCombineV1[i - 1].size() != needToCombineV1[i].size()){
                    allSame = false;
                }
            }
        }
        if(needToCombineV1.size() == 1){
            return;
        }
        if(allSame == true){
            const auto & firstItem = needToCombineV1[0];
            for(int k = 1; k < needToCombineV1.size(); k++){
                const auto & kSelf = needToCombineV1[k];
                if(firstItem != kSelf){
                    allSame = false;
                }
            }
            if(allSame == true){
                uint Matchresult = 1;
                for(int i = 0; i < needToCombineV1.size(); i++){
                    Matchresult *= (firstItem.size() - i);
                }
                return;
            }
        }
        std::vector<std::vector<int>> needToCombine;
        for(int i = 0; i < needToCombineV1.size(); i++){
            const auto & I_isolateVertexCandidate = needToCombineV1[i];
            std::vector<int> I_needToCombine;
            I_needToCombine.reserve(I_isolateVertexCandidate.size());
            int I_NoOverLeafWeight = 0;
            for(int k = 0; k < I_isolateVertexCandidate.size(); k++){
                int vertex = I_isolateVertexCandidate[k];
                if(isolatedVertexMap[vertex] > 1){
                    I_needToCombine.push_back(vertex);
                }
                else{
                    I_NoOverLeafWeight++;
                }
            }
            NoOverLeafWeight.push_back(I_NoOverLeafWeight);
            if(I_NoOverLeafWeight > 0){
                I_needToCombine.push_back(-1);
            }
            needToCombine.push_back(I_needToCombine);
        }
        uint depth = this->queryVec[queryIndex].NumVertices();
        this->headRecord[depth].clear();
        this->intersectionResult[depth].clear();
        this->combineStack[depth].clear();
        for(int i = 0; i < needToCombine.size(); i++){
            headRecord[depth].push_back(needToCombine[i].size());
        }
        this->stackSize[depth] = needToCombine.size();
        this->stackHead[depth] = 0;
        this->combineStack[depth].resize(this->stackSize[depth]);
        this->intersectionResult[depth].resize(this->stackSize[depth]);
        this->type[depth] = finalStack;
        //combine
        uint totalMatch = 1;
        int currentVertex;
        int replaceIndex;
        while(this->headRecord[depth][0] >= 0){
            //1.push back stack
            while(this->stackHead[depth] < this->stackSize[depth]){
                if(this->headRecord[depth][this->stackHead[depth]] == 0){
                    bool pending = this->headChange(needToCombine, depth, totalMatch, NoOverLeafWeight);
                    if(pending == false){
                        return;
                    }
                }
                replaceIndex = this->stackHead[depth];// need to add
                currentVertex = needToCombine[replaceIndex][this->headRecord[depth][replaceIndex] - 1];
                if(currentVertex != -1 && this->visited_[currentVertex] == true){
                    this->headRecord[depth][replaceIndex]--;
                }
                else{
                    this->combinePushBack(currentVertex, replaceIndex, depth, totalMatch, NoOverLeafWeight);
                }
            }
            this->combineStackPopTail(depth, totalMatch, NoOverLeafWeight);
        }
    }
}

void csmpp::searchInitSampleCollect(uint v1, uint v2, uint label, bool op){
    srand((unsigned int)time(NULL));
    const auto &  querySeries = this->UpdateEdgeFineQuery(v1, v2, label, pos);
    std::vector<std::pair<uint, uint>> queryCandidate;
    for(int index = 0; index < querySeries.size(); index++){
        uint queryGraph = querySeries[index].first;
        uint queryEdge = querySeries[index].second;
        Edge _queryEdge = this->queryVec[queryGraph].GetEdge(queryEdge);
        uint v1label = this->data_.GetVertexLabel(v1);
        uint v2label = this->data_.GetVertexLabel(v2);
        uint v1Query = _queryEdge.GetV1();
        uint v2Query = _queryEdge.GetV2();
        if(v1label == v2label){
            if((this->indexCheckForSampleCollect(v1, v1Query, queryGraph) && this->indexCheckForSampleCollect(v2, v2Query, queryGraph))||
            (this->indexCheckForSampleCollect(v1, v2Query, queryGraph) && this->indexCheckForSampleCollect(v2, v1Query, queryGraph))){
                queryCandidate.push_back(querySeries[index]);
            }
        }
        else{
            if(v1label == _queryEdge.GetV1Label()){
                if(this->indexCheckForSampleCollect(v1, v1Query, queryGraph) &&
                    this->indexCheckForSampleCollect(v2, v2Query, queryGraph)
                ) {
                    queryCandidate.push_back(querySeries[index]);
                }
            }
            else{
                if(this->indexCheckForSampleCollect(v2, v1Query, queryGraph) &&
                    this->indexCheckForSampleCollect(v1, v2Query, queryGraph)
                ) {
                    queryCandidate.push_back(querySeries[index]);
                }
            }
        }
    }
    int i = 0;
    for(auto item : queryCandidate){
        if(randomUtil::randomProbability(0.5)){
            continue;
        }
        auto  _edge = this->queryVec[item.first].GetEdge(item.second);
        uint _edgeV1label = _edge.GetV1Label();
        uint _edgeV2label = _edge.GetV2Label();
        const auto & matchOrder = this->queryVec[item.first].GetMatchOrder(item.second);
        int collectDepth = randomUtil::getRandomIntBetweenAandB(2, this->queryVec[item.first].NumVertices() - 1);
        if(_edgeV1label != _edgeV2label){
            if(this->data_.GetVertexLabel(v1) != this->queryVec[item.first].GetVertexLabel(matchOrder[0])){
                std::swap(v1, v2);
            }
            this->matchVertex(v1);
            this->matchVertex(v2);
            this->queryVec[item.first].isolatedVertexTimes.clear();

            this->searchVertexSampleCollect(item.first, item.second, false, collectDepth);

            this->popVertex(v2);
            this->popVertex(v1);
        }
        else{
            for(int i = 0; i < 2; i++){

                if(this->indexCheckForSampleCollect(v1, matchOrder[0], item.first)){
                    this->matchVertex(v1);
                    this->matchVertex(v2);
                    this->queryVec[item.first].isolatedVertexTimes.clear();
                    this->searchVertexSampleCollect(item.first, item.second, false, collectDepth);
                    collectDepth = randomUtil::getRandomIntBetweenAandB(2, this->queryVec[item.first].NumVertices() - 1);
                    this->popVertex(v2);
                    this->popVertex(v1);
                }
                std::swap(v1,v2);
            }
        }
        if(reach_time_limit) return;
    }
}

bool csmpp::indexCheckForSampleCollect(uint data_v, uint query_v, uint queryID){
    const uint queryGraphVLabelMaxNum = this->queryVec[queryID].NumVLabels();
    const uint queryGraphELabelMaxNum = this->queryVec[queryID].NumELabels();
    const uint dataGraphVLabelMaxNum = this->data_.NumVLabels();
    const uint dataGraphELabelMaxNum = this->data_.NumELabels();
    const auto & dataGraph = this->data_.index[data_v];
    const auto & queryGraph = this->queryVec[queryID].index[query_v];
    for(int i = 0; i < queryGraphVLabelMaxNum; i++){
        if(i < dataGraphVLabelMaxNum){
            if(dataGraph[i] < queryGraph[i]){
                return false;
            }
        }
        else{
            if(queryGraph[i] > 0){
                return false;
            }
        }
    }
    for(int i = 0; i < queryGraphELabelMaxNum; i++){
        if(i < dataGraphELabelMaxNum){
            if(dataGraph[dataGraphVLabelMaxNum + i] < queryGraph[queryGraphVLabelMaxNum + i]){
                return false;
            }
        }
        else{
            if(queryGraph[queryGraphVLabelMaxNum + i] > 0){
                return false;
            }
        }
    }
    return true;
}