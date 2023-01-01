/*
 * @Author: zmli6 liziming033@gmail.com
 * @Date: 2022-11-15 21:19:54
 * @LastEditors: zmli6 liziming033@gmail.com
 * @LastEditTime: 2022-11-18 01:55:22
 * @FilePath: /CSM_matchOrderStatusAddLR/DecisionMakingSystem/DecisionMakingSystem.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%
 */
#include "DecisionMakingSystem.h"

void DMS::init(std::vector<Graph> & queryGraphs, Graph & dataGraph, uint featureSize, std::string SampleFilePath){
    this->featureSize = featureSize;
    for(int i = 0; i < queryGraphs.size(); ++i){
        Graph & queryGraph = queryGraphs[i];
        std::vector<std::vector<float *>> OneGraphParams;
        for(int j = 0; j < queryGraph.NumEdges(); ++j){
            std::vector<float *> OneMatchOrderParams;
            const auto & matchOrder = queryGraph.GetMatchOrder(j);
            for(int k = 0; k < matchOrder.size(); ++k){
                float * param;
                param = new float[this->featureSize];
                param[0] = dataGraph.getFrequencyRank(true, queryGraph.GetVertexLabel(matchOrder[k]));//Label rank
                param[1] = matchOrder.size() - k;//size - depth
                param[2] = param[1] / matchOrder.size();//(size - depth)/size
                const auto & destLIst = queryGraph.GetDescList(j, k);
                param[3] = destLIst.size() - 1;//desList.size() - 1
                param[4] = queryGraph.getFreezeVertexNumAfter(j, k);//freeze vertex number after current vertex
                param[5] = queryGraph.getIsolatedVertexNumAfter(j, k);//isolated vertex number after current vertex
                //param[6] candidate list
                //param[7] avagere degree
                OneMatchOrderParams.push_back(param);
            }
            OneGraphParams.push_back(OneMatchOrderParams);
        }
        this->params.push_back(OneGraphParams);
    }

    this->sample_ = SampleFilePath;
}

LRAndIndexCheckType DMS::makeDecision(const vertexType & currentType, const uint & destListSize, const uint & unfreezeListSize){
    if(currentType == isolatedVertex){
        return ioslatedVertexNothing;
    }
    if(unfreezeListSize == 0){
        if(currentType == freeVertex){
            if(this->predictStatus == 0){
                return Part1Nothing;
            }
            if(this->predictStatus == 1){
                return Part1JustCheck;
            }
            //return Part1JustCheck;
            return Part1DoLR;
        }
        return Part1Nothing;
    }
    else{/*need to unfreeze*/
        if(destListSize != 0){
            return Part1JustCheck;
        }
        else{
            if(currentType == freeVertex){
                if(this->predictStatus == 0){
                    return Part2Nothing;
                }
                if(this->predictStatus == 1){
                    return Part2JustCheck;
                }
                //return Part2JustCheck;
                return Part2DoLR;
            }
            return Part2Nothing;
        }
    }
}

int DMS::findTableSizeof2(const int target){
    int temp = target -1;
    temp |= temp >> 1;
    temp |= temp >> 2;
    temp |= temp >> 4;
    temp |= temp >> 8;
    temp |= temp >> 16;
    return (temp < 0) ? 1 : temp + 1;
}

int DMS::LRinit(uint PosSampleCount, uint NegSampleCount){
    std::cout << "PosSampleCount: " << PosSampleCount << std::endl;
    std::cout << "NegSampleCount: " << NegSampleCount << std::endl;
    if(NegSampleCount + PosSampleCount < 0){
        return 1;
    }
    // if(PosSampleCount > NegSampleCount * 10){
    //     return 1;
    // }
    // else if(NegSampleCount > PosSampleCount * 10){
    //     return 0;
    // }
    else{
        std::cout << "can make LR" << std::endl;
        int rowCount = 0;
        std::ifstream sampleFP(this->sample_, std::ios::in);
        std::string rowData;
        std::vector<std::string>rowDatas;
        while(std::getline(sampleFP, rowData)){
            rowDatas.push_back(rowData);
            rowCount++;
        }
        sampleFP.close();
        double Proportion = 0.7;
        int trainCount = (int)floor(rowCount * Proportion);
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(rowDatas.begin(), rowDatas.end(), std::default_random_engine(seed));
    
        //1. make for train set
        std::cout << "make for train set" << std::endl;
        auto data = std::make_unique<ANN::Database>();
        std::string lineData;
        std::stringstream _ss;
        int item = 0;
        while(item < trainCount){
            _ss << rowDatas[item];
            int label;
            _ss >> label;
            data->labels.push_back(label);

            double postime, negtime;
            _ss >> postime;
            data->posTime.push_back(postime);
            _ss >> negtime;
            data->negTime.push_back(negtime);
            double timeDiff;
            _ss >> timeDiff;
            data->timeDiff.push_back(timeDiff);
            std::vector<float> Parameters;
            float parameter;
            int time = featureSize;
            while(time--){
                _ss >> parameter;
                Parameters.push_back(parameter);
            }
            data->samples.push_back(Parameters);
            _ss.clear();
            item++;
        }
        //2. make for predict set
        _ss.clear();
        std::cout << "make for predict set" << std::endl;
        auto predict_ = std::make_unique<ANN::Database>();
        while(item < rowCount){
            _ss << rowDatas[item];
            int label;
            _ss >> label;
            predict_->labels.push_back(label);
            double postime, negtime;
            _ss >> postime;
            predict_->posTime.push_back(postime);
            _ss >> negtime;
            predict_->negTime.push_back(negtime);
            double timeDiff;
            _ss >> timeDiff;
            predict_->timeDiff.push_back(timeDiff);
            std::vector<float> Parameters;
            float parameter;
            int time = featureSize;
            while(time--){
                _ss >> parameter;
                Parameters.push_back(parameter);
            }
            predict_->samples.push_back(Parameters);
            _ss.clear();
            item++;
        }
        std::cout << "train size is " << data->samples.size() << std::endl;
        std::cout << "predict size is " << predict_->samples.size() << std::endl;
        //3. set for LR
        int ret = this->LR.init(std::move(data), std::move(predict_), featureSize, .0001, 100, ANN::Optimization::MBGD, findTableSizeof2(data->samples.size()/20));
        if (ret != 0) {
            fprintf(stderr, "logistic regression init fail: %d\n", ret);
            return -1;
        }
        const std::string model{ "data/logistic_regression.model" };
        ret = this->LR.train(model);
        if (ret != 0) {
            fprintf(stderr, "logistic regression train fail: %d\n", ret);
            return -1;
        }
        this->LR.predict();
        return 2;
    }
}

void DMS::SetPredict(uint PosSampleCount, uint NegSampleCount){
    this->predictStatus = this->LRinit(PosSampleCount, NegSampleCount);
    if(this->predictStatus == 1 || this->predictStatus == 0){
        return;
    } 
    this->LR.load_model("data/logistic_regression.model");
    return;
}

bool DMS::predict(const std::vector<uint> & candidateS, uint queryIndex, uint edgeIndex, uint depth, Graph & dataGraph){
    // this->predictStatusCount++;
    // this->timer.StartTimer();
    this->params[queryIndex][edgeIndex][depth][this->featureSize - 2] = candidateS.size();
    double degree = 0;
    for(int i = 0; i < candidateS.size(); ++i){
        degree += dataGraph.GetDegree(candidateS[i]);
    }
    this->params[queryIndex][edgeIndex][depth][this->featureSize - 1] = degree / candidateS.size();
    if(this->LR.predict(this->params[queryIndex][edgeIndex][depth], this->featureSize) > 0.5){
        // this->timer.StopTimer();
        // this->predict1++;
        return true;
    }
    // this->timer.StopTimer();
    // this->predict2++;
    return false;
}