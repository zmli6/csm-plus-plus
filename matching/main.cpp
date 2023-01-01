#include <chrono>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include "utils/CLI11.hpp"
#include "utils/globals.h"
#include "utils/types.h"

#include "graph/graph.h"
#include "matching/matching.h"
#include "matching/csmpp.h"
void splitFile(std::string& OriginFile, const std::string & NewFile1, const std::string & NewFile2, uint splitNum){
    std::ifstream OriginFileFp(OriginFile, std::ios::in);
    std::ofstream NewFile1Fp(NewFile1, std::ios::out);
    std::ofstream NewFile2Fp(NewFile2, std::ios::out);
    int Count = 0;
    std::string rowData;
    std::vector<std::string>rowDatas;
    while(getline(OriginFileFp, rowData)){
        rowDatas.push_back(rowData);
    }
    OriginFileFp.close();
    int Ndata = rowDatas.size() - splitNum;
    for(int i = 0; i < rowDatas.size(); i++){
        if(i < Ndata){
            NewFile1Fp << rowDatas[i] << std::endl;
        }
        else{
            NewFile2Fp << rowDatas[i] << std::endl;
        }
    }
    NewFile1Fp.close();
    NewFile2Fp.close();
}
int main(int argc, char *argv[])
{
    CLI::App app{"App description"};
    std::string NewAlgorithm = "csmpp";
    uint queryGraphIndex = 0;
    std::string query_path = "", initial_path = "", stream_path = "", FolderPath = "", NewInitail_path = "", sample_path = "", sampleRecordPath = "";
    uint max_num_results = UINT_MAX, time_limit = UINT_MAX, initial_time_limit = UINT_MAX;
    bool print_prep = false, print_enum = false, homo = false, report_initial = true, print_hnu = false, LR = true;
    uint sampleNumUpdateUse = 1; int sampleCollectNum = 1;
    std::vector<std::vector<std::vector<uint>>> orders;

    app.add_option("-q,--query", query_path, "query graph path/Prefixed path")->required();
    app.add_option("-d,--data", initial_path, "initial data graph path")->required();
    app.add_option("-u,--update", stream_path, "data graph update stream path")->required();
    app.add_option("--time-limit", time_limit, "time limit for the incremental matching (second)");
    app.add_option("--print-prep", print_prep, "print processing results or not");
    app.add_option("--print-enum", print_enum, "print enumeration results or not");
    app.add_option("--report-initial", report_initial, "report the result of initial matching or not");
    app.add_option("--initial-time-limit", initial_time_limit, "time limit for the initial matching (second)");
    app.add_option("--I2sn,--Insert2sampleNum", sampleNumUpdateUse, "The initail data use to collect LR sample");
    app.add_option("--SCN,--SampleCollectNum", sampleCollectNum, "LR sample collect number");
    app.add_option("--SRP", sampleRecordPath, "LR sample collect file ");
    CLI11_PARSE(app, argc, argv);
    
    std::chrono::high_resolution_clock::time_point start, lstart;

    start = Get_Time();
    std::cout << "----------- Loading graphs ------------" << std::endl;
    //use for LR
    int index = initial_path.find_last_of("/");
    FolderPath = initial_path.substr(0, index);
    NewInitail_path = FolderPath + "/NdataInit.graph";
    sample_path = FolderPath + "/NdataSample.graph";
    splitFile(initial_path, NewInitail_path, sample_path, sampleNumUpdateUse);
    initial_path = NewInitail_path;
    std::cout << "file split success" << std::endl;
    
    

    Graph data_graph {};
    data_graph.LoadFromFile(initial_path);
    data_graph.PrintMetaData();
    Print_Time("Load Graphs: ", start);

    Graph query_graph {};
    std::vector<Graph> multiQueryGraph;
    query_graph.LoadFromFile(query_path, multiQueryGraph);


    
    std::cout << "------------ Preprocessing ------------" << std::endl;
    matching *mm = nullptr;
    csmpp *csmpp_ = nullptr;

    start = Get_Time();
    mm = csmpp_ = new csmpp(data_graph, query_graph, multiQueryGraph, print_hnu, max_num_results, print_prep, print_enum, homo, sample_path, LR, sampleCollectNum, sampleRecordPath);
    
    mm->Preprocessing();
    Print_Time("Preprocessing: ", start);
    
    if (report_initial)
    {
        std::cout << "----------- Initial Matching ----------" << std::endl;
        
        start = Get_Time();
        auto InitialFun = [&mm]()
        {
            mm->InitialMatching();
        };
        execute_with_time_limit(InitialFun, initial_time_limit, reach_time_limit);
        Print_Time("Initial Matching: ", start);
        
        size_t num_results = 0ul;
        mm->GetNumInitialResults(num_results);
        std::cout << num_results << " matches. \n";

        if (reach_time_limit) return 0;
    }
    
    std::cout << "--------- Incremental Matching --------" << std::endl;
    data_graph.LoadUpdateStream(stream_path);
    size_t num_v_updates = 0ul, num_e_updates = 0ul;
    
    auto IncrementalFun = [&data_graph, &mm, &num_v_updates, &num_e_updates, &start]()
    {
        int i = 0;
        size_t temp = 0;
        while (!data_graph.updates_.empty())
        {
            i++;
            const InsertUnit & insert = data_graph.updates_.front();
            if (insert.type == 'v' && insert.is_add)
            {
                mm->AddVertex(insert.id1, insert.label);
                num_v_updates ++;
            }
            else if (insert.type == 'v' && !insert.is_add)
            {
                mm->RemoveVertex(insert.id1);
                num_v_updates ++;
            }
            else if (insert.type == 'e' && insert.is_add)
            {
                mm->AddEdge(insert.id1, insert.id2, insert.label);
                num_e_updates ++;
            }
            else if (insert.type == 'e' && !insert.is_add)
            {
                mm->RemoveEdge(insert.id1, insert.id2);
                num_e_updates ++;
            }
            data_graph.updates_.pop();
            if (reach_time_limit) break;
        }
    };
    start = Get_Time();

    execute_with_time_limit(IncrementalFun, time_limit, reach_time_limit);
    Print_Time("Incremental Matching: ", start);

    

    std::cout << num_v_updates << " vertex updates.\n";
    std::cout << num_e_updates << " edge updates.\n";

    size_t positive_num_results = 0ul, negative_num_results = 0ul;
    mm->GetNumPositiveResults(positive_num_results);
    mm->GetNumNegativeResults(negative_num_results);
    std::cout << positive_num_results << " positive matches.\n";
    std::cout << negative_num_results << " negative matches.\n";

    mm->PrintCounter();
    mm->TimePrint(true);
    size_t num_edges = 0u, num_vertices = 0ul;
    mm->GetMemoryCost(num_edges, num_vertices);
    std::cout << "\n# edges in index in total: " << num_edges;
    std::cout << "\n# vertices in index in total: " << num_vertices;

    std::cout << "\nPeak Virtual Memory: " << mem::getValue() / 1024.0 << " MB";
    std::cout << "\n\n----------------- End -----------------" << std::endl;
    //delete incIsoMatch;
    
    delete mm;
    
    
}
