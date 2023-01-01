/*
 * @Author: zmli6 liziming033@gmail.com
 * @Date: 2022-11-17 08:17:04
 * @LastEditors: zmli6 liziming033@gmail.com
 * @LastEditTime: 2022-11-18 00:46:01
 * @FilePath: /CSM_matchOrderStatusAddLR/DecisionMakingSystem/DecisionMakingSystem.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef DecisionMakingSystem_DecisionMakingSystem
#define DecisionMakingSystem_DecisionMakingSystem
#include "model/logistic_regression.h"
#include "utils/types.h"
#include "utils/globals.h"
#include "graph/graph.h"
class DMS{
private:
    //use for lr
    std::string sample_;
    std::string LRParasGetFilePath;
    ANN::LogisticRegression LR;
    std::vector<std::vector<std::vector<float *>>> params;
    uint featureSize;
    
public:
    void init(std::vector<Graph> & queryGraphs, Graph & dataGraph, uint featrueSize, std::string SampleFilePath);
    LRAndIndexCheckType makeDecision(const vertexType & currentType, const uint & destListSize, const uint & unfreezeListSize);
    //LR Part
    bool predict(const std::vector<uint> & candidateS, uint queryIndex, uint edgeIndex, uint depth, Graph & dataGraph);
    void SetPredict(uint PosSampleCount, uint NegSampleCount);
    uint predictStatus = -1;
    // uint predictStatusCount = 0;
    // uint predict1 = 0;
    // uint predict2 = 0;
    // Timer timer;
private:
    int LRinit(uint PosSampleCount, uint NegSampleCount);
    int findTableSizeof2(const int target);
};
#endif