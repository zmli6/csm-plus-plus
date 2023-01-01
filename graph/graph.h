#ifndef GRAPH_GRAPH
#define GRAPH_GRAPH

#include <queue>
#include <tuple>
#include <vector>
#include <map>
#include "../utils/types.h"
#include "../utils/utils.h"
#define LABLE_TYPE std::vector<std::pair<std::vector<uint>, std::vector<uint>>> 

class Edge{
private:
    uint v1;
    uint v2;
    uint v1Label;
    uint v2Label;
    uint eLabel;
    uint index;
    bool exist;
public:
    Edge(uint v1, uint v2, uint v1Label, uint v2Label, uint eLabel, uint index):v1(v1),v2(v2),v1Label(v1Label),v2Label(v2Label),eLabel(eLabel),index(index){
        this->exist = true;
    }
    Edge(uint v1Label, uint v2Label, uint eLabel):v1Label(v1Label),v2Label(v2Label),eLabel(eLabel){}
    bool GetExist(){return this->exist;}
    void EdgeDelete(){this->exist = false;}
    const uint GetV1() const{return this->v1;}
    const uint GetV2() const{return this->v2;}
    const uint GetV1Label() const{return this->v1Label;}
    const uint GetV2Label() const{return this->v2Label;}
    const uint GeteLabel() const{return this->eLabel;}
    uint GetIndex(){return this->index;}
    bool operator == (const Edge& edge)const {
        if((edge.v1Label == this->v1Label && edge.v2Label == this->v2Label && edge.eLabel == this->eLabel) ||
        (edge.v2Label == this->v1Label && edge.v1Label == this->v2Label && edge.eLabel == this->eLabel)){
            // std::cout << "v1Label : " << edge.v1Label << " v2Label: " << edge.v2Label << " eLabel: " << edge.eLabel << std::endl;
            // std::cout << "v1Label : " << this->v1Label << " v2Label: " << this->v2Label << " eLabel: " << this->eLabel << std::endl;
            // std::cout << std::endl;
            return true;
        }
        return false;
    }
    bool operator < (const Edge & edge) const{
        if(this->v1Label < edge.v1Label){
            return true;
        } else if (this->v1Label > edge.v1Label) {
            return false;
        }

        if(this->v2Label < edge.v2Label){
            return true;
        } else if (this->v2Label > edge.v2Label) {
            return false;
        }


        if(this->eLabel < edge.eLabel){
            return true;
        } else if (this->eLabel > edge.eLabel) {
            return false;
        }

        return false;
    }
};
class Graph
{
protected:
    uint edge_count_;
    uint vlabel_count_;
    uint elabel_count_;
    std::vector<std::vector<uint>> neighbors_;
    std::vector<std::vector<uint>> elabels_;
    std::vector<std::vector<uint>> matchOrder;
    std::vector<std::vector<uint>> matchIsolatedOrder;
    std::vector<std::vector<vertexType>> matchVertexTypes;
    std::vector<std::vector<std::vector<uint>>> unfreezeRecord;
    std::vector<std::vector<std::vector<std::tuple<uint,uint,uint>>>> descList;//tuple <orderIndex, neighborLabel, elabel>
    std::vector<Edge> edge;
    std::vector<uint> mapping;
    std::vector<std::vector<uint>> freezeVertexNumAfter;
    std::vector<std::vector<uint>> isolatedVertexNumAfter;

    uint * eLabelFrequery;
    uint * vLabelFrequery;

    uint * eLabelIndex2Rank;
    uint * vLabelIndex2Rank;

    uint * eLabelRank2Index;
    uint * vLabelRank2Index;
    
public:
    std::queue<InsertUnit> updates_;
    std::vector<uint> vlabels_;// vertex's label
    bool matchOrderBuild;// build index flag
    LABLE_TYPE data_label; //data graph label check
    std::vector<int*>index;
    //make queryGraph use
    std::vector<uint>vertexLabel;
    std::vector<uint>vertexDataGraphID;
    std::map<uint, uint> isolatedVertexTimes;

public:
    Graph()
    : edge_count_(0)
    , vlabel_count_(0)
    , elabel_count_(0)
    , neighbors_{}
    , elabels_{}
    , updates_{}
    , vlabels_{}
    , matchOrderBuild(false)
    {}
    Graph(bool matchOrderBuild):matchOrderBuild(matchOrderBuild){}

    virtual uint NumVertices() const { return vlabels_.size(); }
    virtual uint NumEdges() const { return edge_count_; }
    uint NumVLabels() const { return vlabel_count_; }
    uint NumELabels() const { return elabel_count_; }
    uint GetDiameter() const;

    void AddVertex(uint id, uint label);
    void RemoveVertex(uint id);
    void AddEdge(uint v1, uint v2, uint label);
    void RemoveEdge(uint v1, uint v2);

    uint GetVertexLabel(uint u) const;
    const std::vector<uint>& GetNeighbors(uint v) const;
    const std::vector<uint>& GetNeighborLabels(uint v) const;
    uint GetDegree(uint v) const;
    std::tuple<uint, uint, uint> GetEdgeLabel(uint v1, uint v2) const;

    void LoadFromFile(const std::string &path);
    void LoadFromFile(const std::string &prefixPath, std::vector<Graph> & queryGraph);
    void LoadUpdateStream(const std::string &path);
    void PrintMetaData() const;
    void PrintGraphNature() const;
    void indexAvg() const;

    const std::vector<uint> &  GetMatchOrder(uint index) const;
    void DescListInit(uint edgeIndex, std::vector<uint> & order, uint visitedVertex);
    std::vector<std::tuple<uint, uint, uint>> GetDescList(uint edgeIndex, uint depth);
    void MatchOrderInit();
    void EdgeMakeOrder(uint edgeIndex, uint v1, uint v2, std::vector<bool> & visit, std::vector<uint> & order);
    const std::vector<uint> & EdgeGetOrder(Edge edge) const;
    const std::vector<uint> & EdgeGetOrder(uint index) const;
    void matchOrderTypeSet();
    const vertexType getVertexType(uint edgeIndex, uint pos) const; 
    const uint getFreezeVertexNumAfter(uint edgeIndex, uint pos) const ;
    const uint getIsolatedVertexNumAfter(uint edgeIndex, uint pos) const;
    const std::vector<vertexType> & getVertexType(uint edgeIndex) const;
    void setVertexFree(uint edgeIndex, uint pos);
    const std::vector<uint> & getUnfreezeList(uint edgeIndex, uint pos) const;
    void setVertexStatus(uint edgeIndex, const std::vector<uint> & vertexs, vertexType type);
    const std::vector<uint>& getIsolatedVertexIndex(uint edgeIndex) const;
    const uint getDestListSize(uint egdeIndex, uint depth) const;
    const uint getUnfreezeListSize(uint egdeIndex, uint depth) const;

    void indexUpdate(uint v1, uint v2, uint label, bool op);//true insert; false delete
    void indexInit();
    void indexPrint();
    void MatchOrderAllPrint();
    void DescListPrint(uint edgeIndex);
    void frequencyRank(bool init, int vIndex, int eIndex, bool add);
    const uint getFrequencyRank(bool eOrv, int index) const;
    const uint getIndexValue(uint v, uint pos) const;

    Edge GetEdge(uint index);
    std::vector<Edge> GetEdge();

    bool mappingAdd(uint vertex, uint beginEdge);
    void mappingPopTail();
    uint GetMappingSize();

    void isolatedVertexTimesAdd(const std::vector<uint> & candidateVertexs);
    void isolatedVertexTimesMinus(const std::vector<uint> & candidateVertexs);


    static bool cmp(const std::pair<int,int> &p1,const std::pair<int,int> &p2){
        return p1.second > p2.second;
    }

    void addVertexAndEdge(uint id,const Graph & dataGraph_){
        this->vertexDataGraphID.push_back(id);
        this->vertexLabel.push_back(dataGraph_.GetVertexLabel(id));
        uint id_ = this->vertexDataGraphID.size() - 1; 
        this->AddVertex(id_, dataGraph_.GetVertexLabel(id));
        std::vector<uint> idNeighbor = dataGraph_.GetNeighbors(id);
        for(uint i = 0; i < id_; i++){
            uint ids = this->vertexDataGraphID[i];
            auto it = std::lower_bound(idNeighbor.begin(), idNeighbor.end(), ids);
            if(*it == ids){
                uint dis = std::distance(idNeighbor.begin(), it);
                this->AddEdge(i, id_, dataGraph_.GetNeighborLabels(id)[dis]);
            }
        }
    }
    std::string toStr(){
        std::string result = "";
        for(int i = 0; i < this->vertexLabel.size(); i++){
            result += "v " + std::to_string(i) + " " + std::to_string(this->vertexLabel[i]) + "\n";
        }
        for(int i = 0; i < this->edge.size(); i++){
            result += "e " + std::to_string(this->edge[i].GetV1()) + " " + std::to_string(this->edge[i].GetV2()) + " " + std::to_string(this->edge[i].GeteLabel()) + "\n";
        }
        return result;
    }
};

#endif //GRAPH_GRAPH
