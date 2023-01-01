#ifndef MATCHING_HNU
#define MATCHING_HNU

#include <vector>
#include <map>

#include "utils/types.h"
#include "utils/globals.h"
#include "graph/graph.h"
#include "matching/matching.h"
#include "DecisionMakingSystem/DecisionMakingSystem.h"
enum searchType{
    init, pos, neg
};
class csmpp : public matching
{
private:
    
    std::vector<std::pair<Edge, std::vector<std::pair<uint,uint>>>> updateEdgeFindQueryEqual;// record the edge update 
    std::map<Edge,std::vector<std::pair<uint,uint>>> updateEdgeFindQuery;
    std::map<Edge,std::vector<std::pair<uint,uint>>> initEdgeFindQuery;
   
    std::vector<int>match;
    std::vector<std::vector<uint>> matchCandidate;
    bool print_data;
    bool motif;
    uint edgeNum;
    uint graphNum;
    uint updateTime;
    sample Sample;
    fileSystem FileSystem;
    DMS decisionMakeSystem;

    //LR part
    std::string sample_ = "data/sample.txt";
    std::string sampleCollectPath = "";
    int predict = 0;
    std::queue<InsertUnit> updates_;
    bool collectBotton;
    
    bool LRButton;
    int sampleNum;

    //combine part
    std::vector<std::vector<std::vector<uint>>> intersectionResult;
    std::vector<std::vector<int>> combineStack;
    std::vector<std::vector<uint>> headRecord;
    std::vector<uint> stackSize;
    std::vector<unFreezeStackType> type;
    std::vector<int> stackHead;
    

public:
    csmpp(Graph& data_graph, Graph& query_graph, std::vector<Graph> & multiQueryGraph, bool print_data, uint max_num_results,
            bool print_prep, bool print_enum, bool homo, std::string sampleCollectPath, bool LR, uint sampleNum, std::string SRP);
    ~csmpp() override{
        // std::cout << "predict count = " << this->decisionMakeSystem.predictStatusCount << std::endl;
        // std::cout << "predict time = " << this->decisionMakeSystem.timer.CountTime() / 1000.0 << " us" << std::endl;
        // std::cout << "predict1 = " << this->decisionMakeSystem.predict1 << std::endl;
        // std::cout << "predict2 = " << this->decisionMakeSystem.predict2 << std::endl;
    };

    void Preprocessing() override;
    void InitialMatching() override;

    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;
    void TimePrint();
 
    

private:
    //LR part
    
    
    void setLRSampleUpdate();

    bool indexCheck(uint data_v, uint query_v, uint queryID);

    bool indexCheckForSampleCollect(uint data_v, uint query_v, uint queryID);

    void EdgeFindQueryInit();
    void initEdgeFindQueryInit();
    std::vector<std::pair<uint, uint>> UpdateEdgeFineQuery(uint v1, uint v2, uint edgeLabel, searchType type);
    void setMatchVertex(const std::vector<uint> & matchingIndex, const std::vector<int> & vertexs);
    void unsetMatchVertex(const std::vector<uint> & matchingIndex);
    void addMatchResult(uint queryIndex, uint edgeIndex, searchType type);
    void searchInit(uint v1, uint v2, uint label, searchType type);
    void searchVertex(uint queryIndex, uint edgeIndex, searchType type);
    void matchVertex(uint vertex);
    void matchVertex(const std::vector<uint> & Candidate);
    void popVertex(uint vertex);
    void popVertex();
    
    void printMatch();

    void searchInitSampleCollect(uint v1, uint v2, uint label, bool op);
    void searchVertexSampleCollect(uint queryIndex, uint edgeIndex, bool collectButton, int collectDepth);
    void addMatchResultSampleCollect(uint queryIndex, uint edgeIndex);

    bool vertexPushCheck(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection, int depth, uint queryIndex, uint queryVertex);
    bool getIntersection(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection, int depth, uint queryIndex, uint queryVertex);
    const std::vector<uint> & getItersectionTop(int depth) const;
    void combineStackPopTail(int depth);
    void combineStackPopTail(int depth, uint & totalMatch, const std::vector<uint>& NoOverLeftWeight);
    void combinePushBack(uint vertex, uint vertexInWhichCandidateNeighbor, int depth);
    void combinePushBack(int vertex, uint vertexInWhichCandidateNeighbor, int depth, uint & totalMatch, const std::vector<uint>& NoOverLeftWeight);
    bool headChange(const std::vector<std::vector<uint>> & needToCombine, int depth);
    bool headChange(const std::vector<std::vector<int>> & needToCombine, int depth, uint & totalMatch, const std::vector<uint>& NoOverLeftWeight);

    void setVisitedPatch(const std::vector<int> & vertex);
    void setUnVisitedPatch(const std::vector<int> & vertex);

    static bool cmp(const std::pair<int,int> &p1,const std::pair<int,int> &p2){
        return p1.second < p2.second;
    }
};

#endif 
