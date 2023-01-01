#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <tuple>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../utils/types.h"
#include "../utils/utils.h"
#include "graph.h"

void Graph::AddVertex(uint id, uint label)
{
    if (id >= vlabels_.size())
    {
        vlabels_.resize(id + 1, NOT_EXIST);
        vlabels_[id] = label;
        neighbors_.resize(id + 1);
        elabels_.resize(id + 1);
    }
    else if (vlabels_[id] == NOT_EXIST)
    {
        vlabels_[id] = label;
    }
    
    vlabel_count_ = std::max(vlabel_count_, label + 1);
    // print graph
    /*std::cout << "labels: ";
    for (uint i = 0; i < vlabels_.size(); i++)
    {
        std::cout << i << ":" << vlabels_[i] << " (";
        for (uint j = 0; j < neighbors_[i].size(); j++)
        {
            std::cout << neighbors_[i][j] << ":" << elabels_[i][j] << " ";
        }
        std::cout << ")" << std::endl;
    }*/
}

void Graph::RemoveVertex(uint id)
{
    vlabels_[id] = NOT_EXIST;
    neighbors_[id].clear();
    elabels_[id].clear();
}

void Graph::AddEdge(uint v1, uint v2, uint label)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower != neighbors_[v1].end() && *lower == v2) return;;
    size_t dis = std::distance(neighbors_[v1].begin(), lower);
    neighbors_[v1].insert(lower, v2);
    elabels_[v1].insert(elabels_[v1].begin() + dis, label);
    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    dis = std::distance(neighbors_[v2].begin(), lower);
    neighbors_[v2].insert(lower, v1);
    elabels_[v2].insert(elabels_[v2].begin() + dis, label);
    if(v1 > v2){
        std::swap(v1, v2);
    }
    Edge edge(v1, v2, this->GetVertexLabel(v1), this->GetVertexLabel(v2), label, this->edge.size());
    this->edge.push_back(edge);
    edge_count_++;
    elabel_count_ = std::max(elabel_count_, label + 1);
    // print graph
    /*std::cout << "labels: ";
    for (uint i = 0; i < vlabels_.size(); i++)
    {
        std::cout << i << ":" << vlabels_[i] << " (";
        for (uint j = 0; j < neighbors_[i].size(); j++)
        {
            std::cout << neighbors_[i][j] << ":" << elabels_[i][j] << " ";
        }
        std::cout << ")" << std::endl;
    }*/
}

void Graph::RemoveEdge(uint v1, uint v2)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower == neighbors_[v1].end() || *lower != v2)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v1].erase(lower);
    elabels_[v1].erase(elabels_[v1].begin() + std::distance(neighbors_[v1].begin(), lower));
    
    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    if (lower == neighbors_[v2].end() || *lower != v1)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v2].erase(lower);
    elabels_[v2].erase(elabels_[v2].begin() + std::distance(neighbors_[v2].begin(), lower));
    // if(v1 > v2){
    //     std::swap(v1, v2);
    // }
    // int i;
    // for(i = 0; i < this->edge.size(); i++){
    //     if(this->edge[i].GetV1() == v1 && this->edge[i].GetV2() == v2){
    //         std::cout << "find index = " << i << std::endl;
    //     }
    // }
    // uint v1label = this->GetVertexLabel(v1);
    // uint v2label = this->GetVertexLabel(v2);
    // uint elabel = std::get<2>(this->GetEdgeLabel(v1, v2));
    // auto it = find(this->edge.begin(), this->edge.end(), Edge(v1, v2, v1label, v2label, elabel, 0));
    // assert(it != this->edge.end());
    // this->edge.erase(it);
}

uint Graph::GetVertexLabel(uint u) const
{
    return vlabels_[u];
}

const std::vector<uint>& Graph::GetNeighbors(uint v) const
{
    return neighbors_[v];
}

const std::vector<uint>& Graph::GetNeighborLabels(uint v) const
{
    return elabels_[v];
}

std::tuple<uint, uint, uint> Graph::GetEdgeLabel(uint v1, uint v2) const
{
    uint v1_label, v2_label, e_label;
    v1_label = GetVertexLabel(v1);
    v2_label = GetVertexLabel(v2);

    const std::vector<uint> *nbrs;
    const std::vector<uint> *elabel;
    uint other;
    if (GetDegree(v1) < GetDegree(v2))
    {
        nbrs = &GetNeighbors(v1);
        elabel = &elabels_[v1];
        other = v2;
    }
    else
    {
        nbrs = &GetNeighbors(v2);
        elabel = &elabels_[v2];
        other = v1;
    }
    
    long start = 0, end = nbrs->size() - 1, mid;
    while (start <= end)
    {
        mid = (start + end) / 2;
        if (nbrs->at(mid) < other)
        {
            start = mid + 1;
        }
        else if (nbrs->at(mid) > other)
        {
            end = mid - 1;
        }
        else
        {
            e_label = elabel->at(mid);
            return {v1_label, v2_label, e_label};
        }
    }
    return {v1_label, v2_label, -1};
}

uint Graph::GetDegree(uint v) const
{
    return neighbors_[v].size();
}

uint Graph::GetDiameter() const
{
    uint diameter = 0;
    for (uint i = 0u; i < NumVertices(); i++)
    if (GetVertexLabel(i) != NOT_EXIST)
    {
        std::queue<uint> bfs_queue;
        std::vector<bool> visited(NumVertices(), false);
        uint level = UINT_MAX;
        bfs_queue.push(i);
        visited[i] = true;
        while (!bfs_queue.empty())
        {
            level++;
            uint size = bfs_queue.size();
            for (uint j = 0u; j < size; j++)
            {
                uint front = bfs_queue.front();
                bfs_queue.pop();

                const auto& nbrs = GetNeighbors(front);
                for (const uint nbr: nbrs)
                {
                    if (!visited[nbr])
                    {
                        bfs_queue.push(nbr);
                        visited[nbr] = true;
                    }
                }
            }
        }
        if (level > diameter) diameter = level;
    }
    return diameter;
}

void Graph::LoadFromFile(const std::string &path)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    char type;
    while (ifs >> type)
    {
        if (type == 't')
        {
            char temp1;
            uint temp2;
            ifs >> temp1 >> temp2;
        }
        else if (type == 'v')
        {
            uint vertex_id, label;
            ifs >> vertex_id >> label;
            AddVertex(vertex_id, label);
        }
        else
        {
            uint from_id, to_id, label;
            ifs >> from_id >> to_id >> label;
            AddEdge(from_id, to_id, label);
        }
    }
    ifs.close();
}

void Graph::LoadUpdateStream(const std::string &path)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    std::string type;
    while (ifs >> type)
    {
        if (type == "v" || type == "-v")
        {
            uint vertex_id, label;
            ifs >> vertex_id >> label;
            updates_.emplace('v', type == "v", vertex_id, 0u, label);
        }
        else
        {
            uint from_id, to_id, label;
            ifs >> from_id >> to_id >> label;
            updates_.emplace('e', type == "e", from_id, to_id, label);
        }
    }
    ifs.close();
    std::cout << "# update size = " << this->updates_.size() << std::endl;
}

void Graph::PrintMetaData() const
{
    std::cout << "# vertices = " << NumVertices() << 
        "  # edges = " << NumEdges() << 
        "  vlabel size = " << this->NumVLabels() << 
        "  elabel size = " << this->NumELabels() << std::endl;
}

/**
 * @description: 获取match order
 * @param {uint} index
 * @return {*}
 */
const std::vector<uint> & Graph::GetMatchOrder(uint index) const{
    return this->matchOrder[index];
}

/**
 * @description: 建立match order的初始化部分
 * @return {*}
 */
void Graph::MatchOrderInit(){
    uint numVertices = this->NumVertices();
    for(uint i = 0; i < this->edge.size(); i++){
        auto edge = this->edge[i];
        std::vector<bool>visit(numVertices, false);
        //std::cout << "match order begin set " << std::endl;
        std::vector<uint> order{};
        EdgeMakeOrder(i, edge.GetV1(), edge.GetV2(), visit, order);
        //std::cout << "order make" << std::endl;
        this->matchOrder.push_back(order);
        /*
        for(int i = 0; i < order.size(); i++){
            std::cout << order[i] << " ";
        }
        std::cout << std::endl;
        */
    }
    /*
    std::cout << "match order set" << std::endl;
    for(int i = 0; i < this->matchOrder.size(); i++){
        for(int j = 0; j < this->matchOrder[i].size(); j++){
            std::cout << this->matchOrder[i][j] << "(";
            std::vector<std::tuple<uint,uint,uint>> des = this->descList[i][j];
            for(auto item : des){
                std::cout << this->matchOrder[i][std::get<0>(item)] << " ";
            }
            std::cout << ")->";
        }
        std::cout << "|" <<  std::endl;
    }
    std::cout << std::endl;
    */
}

/**
 * @description: 根据已成的match order 获取对应的位置的DescList
 * @param {uint} edgeIndex
 * @param {vector<uint>} &
 * @param {uint} visitedVertex
 * @return <orderIndex, vertexlabel, elabel>
 */
void Graph::DescListInit(uint edgeIndex, std::vector<uint> & order, uint visitedVertex){
    //std::cout << "in DescList Init" << std::endl;
    std::vector<std::tuple<uint,uint,uint>> neighborBefore{};
    const std::vector<uint> & VertexNeighbor = this->GetNeighbors(visitedVertex);
    for(int i = 0; i < order.size(); i++){
        if(find(VertexNeighbor.begin(), VertexNeighbor.end(), order[i]) != VertexNeighbor.end()){
            uint vertexLabel = this->GetVertexLabel(order[i]);
            uint eLabel = std::get<2>(this->GetEdgeLabel(visitedVertex, order[i]));
            neighborBefore.push_back(std::tuple(i, vertexLabel, eLabel));
        }
    }
    if(this->descList.size() == edgeIndex){
        std::vector<std::vector<std::tuple<uint,uint,uint>>> edgeIndexMatchDesList {};
        edgeIndexMatchDesList.push_back(neighborBefore);
        this->descList.push_back(edgeIndexMatchDesList);
    }
    else{
        this->descList[edgeIndex].push_back(neighborBefore);
    }
}

/**
 * @description: 根据触发边获取match order
 * @param {uint} edgeIndex
 * @param {uint} v1
 * @param {uint} v2
 * @param {vector<bool>} &
 * @param {vector<uint>} &
 * @return {*}
 */
void Graph::EdgeMakeOrder(uint edgeIndex, uint v1, uint v2, std::vector<bool> & visit, std::vector<uint> & order){
    
    // std::cout << "length = " << order.size() << std::endl;
    // for(int i = 0; i < order.size(); i++){
    //     std::cout << order[i] << " ";
    // }
    // std::cout << std::endl;
    // std::cout << std::endl;
    uint vMax = this->NumVertices();
    if(order.size() == vMax){
        return;
    } 
    if(order.size() < 2){
        if(this->GetDegree(v1) > this->GetDegree(v2)){
            std::swap(v1, v2);
        }
        this->DescListInit(edgeIndex, order, v1);
        order.push_back(v1);
        visit[v1] = true;

        this->DescListInit(edgeIndex, order, v2);
        order.push_back(v2);
        visit[v2] = true;
        
        EdgeMakeOrder(edgeIndex, v1, v2, visit, order);
    }
    else{
        //round 1 check neighbor has been visited
        std::map<uint,uint>compareUse;
        /*
        {
            std::cout << " visit vector" << std::endl;
            for(auto item : visit){
                std::cout << item << " ";
            }
            std::cout << std::endl;
        }
        */

        for(uint i = 0; i < vMax; i++){
            if(visit[i] == false){
                std::vector<uint>neighbor = this->GetNeighbors(i);
                for(auto item : order){
                    if(find(neighbor.begin(), neighbor.end(), item) != neighbor.end()){
                        compareUse[i]++;
                    }
                }
            }
        }
        /*
        {
            for(auto iter : compareUse){
                std::cout << iter.first << "->" << iter.second << std::endl;
            }
        }
        */
        std::vector<std::pair<uint,uint> > arr;
        for (std::map<uint, uint>::iterator it=compareUse.begin();it!=compareUse.end();++it)
        {
            arr.push_back(std::make_pair(it->first,it->second));
        }
        sort(arr.begin(),arr.end(),cmp);
        std::vector<uint>Candidate;
        uint max = 0;
        for(auto item : arr){
            if(item.second >= max){
                max = item.second;
                Candidate.push_back(item.first);
            }
            else break;
        }
        if(Candidate.size() == 1){
            this->DescListInit(edgeIndex, order, Candidate[0]);
            visit[Candidate[0]] = true;
            order.push_back(Candidate[0]);
            EdgeMakeOrder(edgeIndex, v1, v2, visit, order);
            return;
        }
        //std::cout << "1" << std::endl;
        
        // {
        //     for(auto iter : Candidate){
        //         std::cout << iter << " " ;
        //     }
        //     std::cout << std::endl;
        // }

        //std::cout << "round 2 finish" << std::endl;
        
        //round 2 check degree
        //std::cout << "before 2 candidate size = " << Candidate.size() << std::endl;
        compareUse.clear();
        arr.clear();
        max = 0;
        for(uint i = 0; i < Candidate.size(); i++){
            compareUse[Candidate[i]] = this->GetDegree(Candidate[i]);
        }
        for (std::map<uint, uint>::iterator it=compareUse.begin();it!=compareUse.end();++it)
        {
            arr.push_back(std::make_pair(it->first,it->second));
        }
        sort(arr.begin(),arr.end(),cmp);
        Candidate.clear();
        max = 0;
        for(auto item : arr){
            if(item.second >= max){
                max = item.second;
                Candidate.push_back(item.first);
            }
            else break;
        }
        if(Candidate.size() == 1){
            this->DescListInit(edgeIndex, order, Candidate[0]);
            visit[Candidate[0]] = true;
            order.push_back(Candidate[0]);
            EdgeMakeOrder(edgeIndex, v1, v2, visit, order);
            return;
        }
        //std::cout << "2" << std::endl;
        //round 3 check id
        //std::cout << "candidate size = " << Candidate.size() << std::endl;
        sort(Candidate.begin(), Candidate.end());
        this->DescListInit(edgeIndex, order, Candidate[0]);
        visit[Candidate[0]] = true;
        order.push_back(Candidate[0]);
        EdgeMakeOrder(edgeIndex, v1, v2, visit, order);
        //std::cout << "3" << std::endl;
        return;
    }
    
}

/**
 * @description: 根据edge获取对应的match order
 * @param {Edge} edge
 * @return {*}
 */
const std::vector<uint> & Graph::EdgeGetOrder(Edge edge) const{
    return this->matchOrder[edge.GetIndex()];
}

/**
 * @description: 根据边的index获取对应的match order
 * @param {uint} index
 * @return {*}
 */
const std::vector<uint> & Graph::EdgeGetOrder(uint index) const{
    return this->matchOrder[index];
}

/**
 * @description: 标签分布初始化
 * @return {*}
 */
void Graph::indexInit(){
    uint vLabelMax = this->NumVLabels();
    uint eLabelMax = this->NumELabels();
    uint vMax = this->NumVertices();
    this->index.resize(vMax, nullptr);
    this->eLabelFrequery = new uint [eLabelMax]();
    this->vLabelFrequery = new uint [vLabelMax]();
    this->eLabelIndex2Rank = new uint [eLabelMax]();
    this->vLabelIndex2Rank = new uint [vLabelMax]();
    this->eLabelRank2Index = new uint [eLabelMax]();
    this->vLabelRank2Index = new uint [vLabelMax]();
    //this->data_label.resize(vMax);
    for(uint k = 0; k < vMax; k++){
        //std::vector<uint>vLabelCount(vLabelMax, 0);
        //std::vector<uint>eLabelCount(eLabelMax, 0);
        std::vector<uint>iNeighbor = this->GetNeighbors(k);
        std::vector<uint>iNeighborEdgeLabel = this->GetNeighborLabels(k);
        int * indexDistribution = new int [vLabelMax + eLabelMax]();
        // for(int j = 0; j < vLabelMax + eLabelMax; j++){
        //     indexDistribution[j] = 0;
        // }//初始化
        //std::cout << "length: " << vLabelMax + eLabelMax << std::endl;
        for(int i = 0; i < iNeighbor.size(); i++){
            /*
            std::cout << "vLabel index: " << this->GetVertexLabel(iNeighbor[i]) << std::endl;
            std::cout << "vLabel list before:  " <<  indexDistribution[this->GetVertexLabel(iNeighbor[i])] << std::endl;
            */
            indexDistribution[this->GetVertexLabel(iNeighbor[i])]++;
            this->vLabelFrequery[this->GetVertexLabel(iNeighbor[i])]++;
            /*
            std::cout << "vLabel list: " <<  indexDistribution[this->GetVertexLabel(iNeighbor[i])] << std::endl;

            std::cout << "eLabel index: " << vLabelMax + iNeighborEdgeLabel[i] << std::endl;
            std::cout << "eLabel list before: " <<  indexDistribution[vLabelMax + iNeighborEdgeLabel[i]] << std::endl;
            */
            indexDistribution[vLabelMax + iNeighborEdgeLabel[i]]++;
            //this->eLabelFrequery[iNeighborEdgeLabel[i]]++;
            /*
            std::cout << "eLabel list: " <<  indexDistribution[vLabelMax + iNeighborEdgeLabel[i]] << std::endl;
            */
        }
        this->index[k] = indexDistribution;
        //std::cout << "index's size = " << this->index.size() << std::endl;
        /*
        {
            std::cout << "vMax" << this->NumVertices() << std::endl;
            std::cout << "vLabelMax" << this->NumVLabels() << std::endl;
            std::cout << "eLabelMax" << this->NumELabels() << std::endl;
            for(int k = 0; k < this->NumVertices(); k++){
                std::cout << "vertex id : " << k << std::endl;
                std::cout << "vLabel: " ;
                for(int i = 0; i < vLabelMax + eLabelMax; i++){
                    std::cout << this->index[k][i] << " ";
                    if(i == vLabelMax - 1){
                        std::cout << std::endl;
                        std::cout << "eLabel: ";
                    }
                }
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
        */
    }
    this->frequencyRank(true, 0, 0, 0);
}

/**
 * @description: 频率排序
 * @return {*}
 */
void Graph::frequencyRank(bool init, int vIndex, int eIndex, bool add){
    uint vLabelMax = this->NumVLabels();
    uint eLabelMax = this->NumELabels();
    // std::cout << "vLabelFrequency->";
    // for(int i = 0; i < vLabelMax; i++){
    //     std::cout << "[" << i << "->" << this->vLabelFrequery[i] << "]" << ",";
    // }
    // std::cout << std::endl;
    //vlabel
    uint * vLabelFrequeryCopy = new uint [vLabelMax];
    memcpy(vLabelFrequeryCopy, this->vLabelFrequery, vLabelMax * sizeof(uint));
    std::sort(vLabelFrequeryCopy, vLabelFrequeryCopy + vLabelMax, std::greater<uint>());
    // std::cout << "vLabelFrequeryCopy->";
    // for(int i = 0; i < vLabelMax; i++){
    //     std::cout << "[" << i << "->" << vLabelFrequeryCopy[i] << "]" << ",";
    // }
    // std::cout << std::endl;
    for(int i = 0; i < vLabelMax; i++){
        uint findKey = this->vLabelFrequery[i];
        uint pos = 0;
        while(pos < vLabelMax){
            if(vLabelFrequeryCopy[pos] != findKey){
                pos++;
            }
            else{
                break;
            }
        }
        this->vLabelIndex2Rank[i] = pos + 1;
    }
    // {
    //     for(int i = 0; i < vLabelMax; ++i){
    //         std::cout << this->vLabelIndex2Rank[i] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    /*
    //elabel
    uint * eLabelFrequeryCopy = new uint [eLabelMax];
    memcpy(eLabelFrequeryCopy, this->eLabelFrequery, eLabelMax * sizeof(uint));
    std::sort(eLabelFrequeryCopy, eLabelFrequeryCopy + eLabelMax, std::greater<uint>());
    for(int i = 0; i < eLabelMax; i++){
        auto it = std::lower_bound(eLabelFrequeryCopy, eLabelFrequeryCopy + eLabelMax, this->eLabelFrequery[i]);
        uint dis = std::distance(eLabelFrequeryCopy, it);
        this->eLabelIndex2Rank[i] = dis + 1;
    }
    */
}

/**
 * @description: 获取频率rank
 * @param {bool} eOrv
 * @param {int} index
 * @return {*}
 */
const uint Graph::getFrequencyRank(bool eOrv, int index) const{
    if(eOrv == true){
        return this->vLabelIndex2Rank[index];
    }
    return this->eLabelIndex2Rank[index];
}

/**
 * @description: 更新标签分布
 * @param {uint} v1
 * @param {uint} v2
 * @param {uint} label
 * @param {bool} op
 * @return {*}
 */
void Graph::indexUpdate(uint v1, uint v2, uint label, bool op){
    uint v1Label = this->GetVertexLabel(v1);
    uint v2Label = this->GetVertexLabel(v2);
    uint Label = this->NumVLabels() + label;
    if(op == true){
        this->index[v1][v2Label]++;
        this->index[v1][Label]++;
        this->index[v2][v1Label]++;
        this->index[v2][Label]++;
    }
    else{
        this->index[v1][v2Label]--;
        this->index[v1][Label]--;
        this->index[v2][v1Label]--;
        this->index[v2][Label]--;
    }
}


/**
 * @description: 根据index获取边元素
 * @param {uint} index
 * @return {*}
 */
Edge Graph::GetEdge(uint index){
    return this->edge[index];
}

/**
 * @description: 获取所有的边元素
 * @return {*}
 */
std::vector<Edge> Graph::GetEdge(){
    return this->edge;
}

/**
 * @description: 增加match mapping 元素 并检测是否已经完全匹配完成
 * @param {uint} vertex
 * @param {uint} beginEdge
 * @return {*}
 */
bool Graph::mappingAdd(uint vertex, uint beginEdge){
    this->mapping.push_back(vertex);
    std::cout << "# Mapping Add size is : " << this->mapping.size() << std::endl;
    if(this->mapping.size() == this->NumVertices()){
        std::cout << "match success" << std::endl;
        for(uint i = 0; i < this->NumVertices(); i++){
            std::cout << "query vertex: " << this->GetMatchOrder(beginEdge)[i] << "--> data graph: " << this->mapping[i] << std::endl;
        }
        return true;
    }
    return false;
}

/**
 * @description: 删除match mapping 的最后一位
 * @return {*}
 */
void Graph::mappingPopTail(){
    this->mapping.pop_back();
}

/**
 * @description: 获取对应触发边的对应的match order中的第index个元素的前置位邻居
 * @param {uint} edgeIndex
 * @param {uint} depth
 * @return {*}
 */
std::vector<std::tuple<uint, uint, uint>> Graph::GetDescList(uint edgeIndex, uint depth){
    //std::cout << "descList size = " << this->descList[edgeIndex].size() << " depth: " << depth << std::endl;
    return this->descList[edgeIndex][depth];
}



/**
 * @description: 批量读取查询图
 * @param {string} &prefixPath
 * @param {uint} prefixBegin
 * @param {uint} prefixEnd
 * @param {vector<Graph>} &
 * @return {*}
 */
void Graph::LoadFromFile(const std::string &prefixPath, std::vector<Graph> & queryGraph){
    std::vector<std::string>fileNames;
    struct stat s;
    stat(prefixPath.c_str(), &s);
    DIR* open_dir = opendir(prefixPath.c_str());
    if (NULL == open_dir) {
        std::exit(EXIT_FAILURE);
    }
    dirent* p = nullptr;
    while( (p = readdir(open_dir)) != nullptr) {
        struct stat st;
        if (p->d_name[0] != '.') {
            std::string name = prefixPath + std::string("/") + std::string(p->d_name);
            stat(name.c_str(), &st);
            fileNames.push_back(name);
        }
    }
    closedir(open_dir);
    /*
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(fileNames.begin(), fileNames.end(), std::default_random_engine(seed));
    
    int i = 0;
    while(i < 100){
        Graph queryGraphItem;
        queryGraphItem.LoadFromFile(fileNames[i]);
        queryGraph.push_back(queryGraphItem);
        i++;
    }
    */
    
    for(int i = 0; i < fileNames.size(); i++){
        Graph queryGraphItem;
        std::cout << "index " << i <<" -> name " << fileNames[i] << std::endl;
        queryGraphItem.LoadFromFile(fileNames[i]);
        queryGraph.push_back(queryGraphItem);
    }
    
    std::cout << "Load query graph number: " << queryGraph.size() << std::endl;
}

void Graph::indexPrint(){
    int vLabelMax = this->NumVLabels();
    int eLabelMax = this->NumELabels();
    int vMAx = this->NumVertices();
    for(int i = 0; i < vMAx; i++){
        std::cout << "vLabel: ";
        for(int k = 0; k < vLabelMax; k++){
            std::cout << this->index[i][k] << " ";
        }
        std::cout << std::endl;
        std::cout << "eLabel: ";
        for(int k = 0; k < eLabelMax; k++){
            std::cout << this->index[i][k + vLabelMax] << " ";
        }
        std::cout << std::endl;
        std::cout << std::endl;

    }
}

void Graph::MatchOrderAllPrint(){
    for(int i = 0; i < this->NumEdges(); i++){
        std::cout << "egde index: " << i << std::endl;
        // i is edge index
        const auto & matchOrder = this->matchOrder[i];
        const auto & vertexTypes = this->matchVertexTypes[i];
        const auto & desitemList = this->descList[i];
        const auto & unfreezeRecord_ = this->unfreezeRecord[i];
        //2.cout a match order
        for(int k = 0; k < matchOrder.size(); k++){
            std::cout << matchOrder[k] << " ";
        }
        std::cout << std::endl;
        for(int k = 0; k < matchOrder.size(); k++){
            std::cout << "vertex : " << matchOrder[k] << std::endl;
            std::cout << "vertex index: " << k << std::endl;
            std::cout << "vertex type: " << vertexTypes[k] << std::endl;
            const auto & unfreezerecord_ = unfreezeRecord_[k];
            std::cout << "unfreeze list : ";
            for(auto vertex : unfreezerecord_){
                std::cout << vertex << " ";
            }
            std::cout << std::endl;
            const auto & desitemList_ = desitemList[k];
            std::cout << "desItem list : " << std::endl;
            for(auto item : desitemList_){
                std::cout << "index: " << std::get<0>(item) << " item order index: " << this->matchOrder[i][std::get<0>(item)] << " item label: " << std::get<1>(item) << " elabel: " << std::get<2>(item) << std::endl ;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }
}

void Graph::DescListPrint(uint edgeIndex){
    std::vector<std::vector<std::tuple<uint,uint,uint>>> temp = this->descList[edgeIndex];
    std::cout << "desList size = " << temp.size() << std::endl;
    auto matchOrder = this->GetMatchOrder(edgeIndex);
    for(int index = 0; index < this->NumVertices(); index++){
        std::vector<std::tuple<uint,uint,uint>> temp_ = temp[index];
        std::cout << "match order index: " << index << " item is " << matchOrder[index] << " descList size is " << temp_.size() << std::endl;
        for(auto item : temp_){
            std::cout << "item in order : " << this->matchOrder[edgeIndex][std::get<0>(item)] << " item label:" << std::get<1>(item) << " elabel: " << std::get<2>(item) << std::endl ;
        }
        std::cout << std::endl;
    }
}

uint Graph::GetMappingSize(){
    return this->mapping.size();
}

void Graph::PrintGraphNature() const{
    double avgDegree = 0;
    uint maxDegree = 0;
    uint DegreeBiggerAvg = 0;
    const uint vMax = this->NumVertices();
    for(int i = 0; i < vMax; i++){
        uint tempDegree = this->GetDegree(i);
        avgDegree  += tempDegree;
        maxDegree = maxDegree > tempDegree? maxDegree : tempDegree;
    }
    avgDegree = avgDegree / vMax;
    for(int i = 0; i < vMax; i++){
        if(this->GetDegree(i) > avgDegree){
            DegreeBiggerAvg++;
        }
    }
    std::cout << "# avgDegree = " << avgDegree << std::endl;
    std::cout << "# maxDegree = " << maxDegree << std::endl;
    std::cout << "# DegreeBiggerAvg++ = " << DegreeBiggerAvg << std::endl;
}

void Graph::indexAvg() const{
    uint vLabelMax = this->NumVLabels();
    uint eLabelMax = this->NumELabels();
    const uint vMax = this->NumVertices();
    double * avg = new double [vLabelMax + eLabelMax]();
    for(int i =  0; i < vMax; i++){
        for(int k = 0; k < vLabelMax + eLabelMax; k++){
            avg[k] = avg[k] + this->index[i][k] ;
        }
    }
    std::cout << "avg index dirtribution: ";
    for(int k = 0; k < vLabelMax + eLabelMax; k++){
        avg[k] = avg[k] / vMax;
        std::cout << avg[k] << " ";
    }
    std::cout << std::endl;
}

const uint Graph::getIndexValue(uint v, uint pos) const{
    return this->index[v][pos];
}

void Graph::matchOrderTypeSet(){
    // set vertexTypes
    for(int edgeIndex = 0; edgeIndex < this->NumEdges(); edgeIndex++){
        const auto & matchOrder = this->matchOrder[edgeIndex];
        const auto & desitemList = this->descList[edgeIndex];
        //0. set unfreeze record;
        std::vector<std::vector<uint>> currentMatchOrderUnfreezeRecord;
        currentMatchOrderUnfreezeRecord.resize(matchOrder.size(), {});
        //1. get vertexTypes and get unfreeze record
        std::vector<vertexType> currentMatchOrderTypes;
        for(int i = 0; i < matchOrder.size(); i++){
            vertexType iType = isolatedVertex;
            if(i == 0 || i == 1){
                iType = freeVertex;
            }
            else{
                for(int k = i + 1; k < matchOrder.size(); k++){
                    const auto & desitemListTuples = desitemList[k];
                    for(int j = 0; j < desitemListTuples.size(); j++){
                        if(std::get<0>(desitemListTuples[j]) == i){
                            if(k == i + 1){
                                iType = freeVertex;
                            }
                            if(k != i + 1){
                                currentMatchOrderUnfreezeRecord[k].push_back(i); // not the after vertex of i , should record it can unfreeze which vertices
                                iType = freezeVertex;
                            }
                            break;
                        }
                    }
                    if(iType == freeVertex || iType == freezeVertex){
                        break;
                    }
                }
            }
            currentMatchOrderTypes.push_back(iType);
        }
        this->matchVertexTypes.push_back(currentMatchOrderTypes);
        this->unfreezeRecord.push_back(currentMatchOrderUnfreezeRecord);
    }
    //delete freeze vertex from unfreeze vertex in desItemList
    for(int edgeIndex = 0; edgeIndex < this->NumEdges(); edgeIndex++){// for every match order
        const auto & unfreezeRecord_ = this->unfreezeRecord[edgeIndex];// get match order's unfreezeRecord
        for(int i = 0; i < unfreezeRecord_.size(); i++){
            if(unfreezeRecord_[i].size() != 0){//now unfreezeRecord_[i] is a unfreeze vertex
                auto & desitemItemOfi = this->descList[edgeIndex][i];// get unfreezeRecord_[i]'s desitem
                for(int k = 0; k < unfreezeRecord_[i].size(); k++){// loop of unfreezeRecord_[i]'s unfreeze record
                    auto iter = desitemItemOfi.begin();
                    while(iter != desitemItemOfi.end()){
                        if(std::get<0>(*iter) == unfreezeRecord_[i][k]){
                            desitemItemOfi.erase(iter);//delete freeze vertex from unfreeze vertex 's desitemList
                            break;
                        }
                        iter++;
                    }
                }
            }
        }
    }
    //record isolate vertex
    for(int i = 0; i < this->NumEdges(); i++){
        const auto & vertexTypes = this->matchVertexTypes[i];
        std::vector<uint>isolateRecord;
        for(int k = 0; k < vertexTypes.size(); k++){
            if(vertexTypes[k] == isolatedVertex){
                isolateRecord.push_back(k);
            }
        }
        this->matchIsolatedOrder.push_back(isolateRecord);
    }
    //record freezeVertex + isolatedVertex Num After
    for(int i = 0; i < this->NumEdges(); i++){
        const auto & vertexTypes = this->matchVertexTypes[i];
        std::vector<uint> freezeVertexRecorder;
        std::vector<uint> isolatedVertexRecorder;
        freezeVertexRecorder.resize(vertexTypes.size());
        isolatedVertexRecorder.resize(vertexTypes.size());
        uint freezeCount = 0;
        uint isolateCount = 0;
        for(int k = vertexTypes.size() - 1; k >= 0; k--){
            freezeVertexRecorder[k] = freezeCount;
            isolatedVertexRecorder[k] = isolateCount;
            if(vertexTypes[k] == freezeVertex){
                freezeCount++;
            }
            if(vertexTypes[k] == isolatedVertex){
                isolateCount++;
            }
        }
        this->freezeVertexNumAfter.push_back(freezeVertexRecorder);
        this->isolatedVertexNumAfter.push_back(isolatedVertexRecorder);
    }
}

const vertexType Graph::getVertexType(uint edgeIndex, uint pos) const{
    return this->matchVertexTypes[edgeIndex][pos];
}

const std::vector<uint> & Graph::getUnfreezeList(uint edgeIndex, uint pos)const {
    return this->unfreezeRecord[edgeIndex][pos];
}

void Graph::setVertexFree(uint edgeIndex, uint pos){
    this->matchVertexTypes[edgeIndex][pos] = freeVertex;
}

void Graph::setVertexStatus(uint edgeIndex, const std::vector<uint> & vertexs, vertexType type){
    for(int i = 0; i < vertexs.size(); i++){
        this->matchVertexTypes[edgeIndex][vertexs[i]] = type;
    }
}

const std::vector<uint> & Graph::getIsolatedVertexIndex(uint edgeIndex) const{
   return this->matchIsolatedOrder[edgeIndex];
}

const std::vector<vertexType> & Graph::getVertexType(uint edgeIndex) const{
    return this->matchVertexTypes[edgeIndex];
}

void Graph::isolatedVertexTimesAdd(const std::vector<uint> & candidateVertexs){
    for(int i = 0; i < candidateVertexs.size(); i++){
        this->isolatedVertexTimes[candidateVertexs[i]]++;
    }
}

void Graph::isolatedVertexTimesMinus(const std::vector<uint> & candidateVertexs){
    for(int i = 0; i < candidateVertexs.size(); i++){
        this->isolatedVertexTimes[candidateVertexs[i]]--;
    }
}

const uint Graph::getFreezeVertexNumAfter(uint edgeIndex, uint pos) const{
    return this->freezeVertexNumAfter[edgeIndex][pos];
}

const uint Graph::getIsolatedVertexNumAfter(uint edgeIndex, uint pos) const{
    return this->isolatedVertexNumAfter[edgeIndex][pos];
}

const uint Graph::getDestListSize(uint egdeIndex, uint depth) const{
    return this->descList[egdeIndex][depth].size();
}

const uint Graph::getUnfreezeListSize(uint edgeIndex, uint depth) const{
    return this->unfreezeRecord[edgeIndex][depth].size();
}
