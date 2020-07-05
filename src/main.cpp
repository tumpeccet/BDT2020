//
// Created by nera on 17. 03. 2020..
//

#include "main.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <bits/unique_ptr.h>
#include <bioparser/bioparser.hpp>
#include "spoa/spoa.hpp"
#include "../alignment/alignment.h"
#include "../KMeans/k_mean.h"
#include <fstream>

using namespace std;


// define a class for sequences in FASTA format
class FASTAformat {
public:
    std::string name;
    std::uint32_t name_length;
    std::string sequence;
    std::uint32_t sequence_length;
    FASTAformat(
            const char* name, std::uint32_t name_length,
            const char* sequence, std::uint32_t sequence_length):
        name{std::string(name, name_length)},
        sequence{std::string(sequence, sequence_length)} {}

};


// define a class for sequences in FASTQ format
class FASTQformat {
public:
    std::string name;
    std::string sequence;
    std::string quality;
    FASTQformat(
            const char* name, std::uint32_t name_length,
            const char* sequence, std::uint32_t sequence_length,
            const char* quality, std::uint32_t quality_length):
            name{std::string(name, name_length)},
            sequence{std::string(sequence, sequence_length)},
            quality{std::string(quality, quality_length)} {}

};


int main(int argc, char* argv[]){
    std::vector<std::unique_ptr<FASTQformat>> fastq_objects;
    std::vector<std::string> filtered_objects;
    bool invalidArgument = false;
    switch (argc) {
        case 1: {
            fprintf(stderr, "No arguments passed to program!\n");
            break;
        }

        case 9: {
            ofstream evaluationFile;
            cout << "Evaluation" << endl;
            int eval_param = stoi(argv[1]);
            string filePath1 = argv[2], filePath2 = argv[3];
            int match = stoi(argv[4]), mismatch = stoi(argv[5]), gap = stoi(argv[6]);
            int length = stoi(argv[8]);
            int alignment_type = stoi(argv[7]);
            alignment::AlignmentType alignment_al;
            switch (alignment_type) {
                case 1: {
                    alignment_al = alignment::AlignmentType::kGlobal;
                    break;
                }
                case 2: {
                    alignment_al = alignment::AlignmentType::kSemiGlobal;
                    break;
                }
                default: {
                    fprintf(stderr, "Invalid alignment type argument!\n");
                    cout << "Please choose one of the following:\n" <<
                         "   1 for global alignment\n" <<
                         "   2 for semi-global alignment\n" << endl;
                    invalidArgument = true;
                }
            }

            if (invalidArgument) break;

            cout << filePath1 << endl << filePath2 << endl;

            switch (eval_param) {
                case 1: {
                    //aleli medusobno
                    std::vector<std::unique_ptr<FASTAformat>> fasta_objects_1, fasta_objects_2;
                    auto fasta_parser_1 = bioparser::createParser<bioparser::FastaParser, FASTAformat>(filePath1);
                    fasta_parser_1->parse(fasta_objects_1, -1);
                    auto fasta_parser_2 = bioparser::createParser<bioparser::FastaParser, FASTAformat>(filePath2);
                    fasta_parser_2->parse(fasta_objects_2, -1);
                    evaluationFile.open("../evaluation.txt", ios::out | ios::app);
                    evaluationFile << filePath1 << endl << filePath2 << endl;

                    for (int i = 0; i < fasta_objects_1.size() - 1; i++) {
                        for (int j = 1; j < fasta_objects_1.size(); j++) {
                            if (i == j) continue;
                            auto score = alignment::Align(fasta_objects_1[i]->sequence.c_str(), fasta_objects_1[i]->sequence.length(),
                                                          fasta_objects_1[j]->sequence.c_str(), fasta_objects_1[j]->sequence.length(), match, mismatch, gap,
                                                          alignment_al);
                            int smaller = fasta_objects_1[i]->sequence.length();
                            if (fasta_objects_1[j]->sequence.length() < smaller)
                                smaller = fasta_objects_1[j]->sequence.length();
                            //std::cout << smaller << ',' << score << std::endl;
                            evaluationFile << "Distance between " <<
                                      fasta_objects_1[i]->name << " and " << fasta_objects_1[j]->name << " is " << smaller - score << std::endl;
                            std::cout << "Distance between " <<
                                      fasta_objects_1[i]->name << " and " << fasta_objects_1[j]->name << " is " << smaller - score << std::endl;
                        }
                    }

                    for (int i = 0; i < fasta_objects_1.size(); i++) {
                        for (int j = 0; j < fasta_objects_2.size(); j++) {
                            auto score = alignment::Align(fasta_objects_1[i]->sequence.c_str(), fasta_objects_1[i]->sequence.length(),
                                                          fasta_objects_2[j]->sequence.c_str(), fasta_objects_2[j]->sequence.length(), match, mismatch, gap,
                                                          alignment_al);
                            int smaller = fasta_objects_1[i]->sequence.length();
                            if (fasta_objects_2[j]->sequence.length() < smaller)
                                smaller = fasta_objects_2[j]->sequence.length();

                            evaluationFile << "Distance between " <<
                                           fasta_objects_1[i]->name << " and " << fasta_objects_2[j]->name << " is " << smaller - score << std::endl;
                            std::cout << "Distance between " <<
                                      fasta_objects_1[i]->name << " and " << fasta_objects_2[j]->name << " is " << smaller - score << std::endl;

                        }
                    }

                    for (int i = 0; i < fasta_objects_2.size() - 1; i++) {
                        for (int j = 1; j < fasta_objects_2.size(); j++) {
                            if (i == j) continue;
                            auto score = alignment::Align(fasta_objects_2[i]->sequence.c_str(), fasta_objects_2[i]->sequence.length(),
                                                          fasta_objects_2[j]->sequence.c_str(), fasta_objects_2[j]->sequence.length(), match, mismatch, gap,
                                                          alignment_al);
                            int smaller = fasta_objects_2[i]->sequence.length();
                            if (fasta_objects_2[j]->sequence.length() < smaller)
                                smaller = fasta_objects_2[j]->sequence.length();

                            evaluationFile << "Distance between " <<
                                           fasta_objects_2[i]->name << " and " << fasta_objects_2[j]->name << " is " << smaller - score << std::endl;
                            std::cout << "Distance between " <<
                                      fasta_objects_2[i]->name << " and " << fasta_objects_2[j]->name << " is " << smaller - score << std::endl;
                        }
                    }

                    evaluationFile << "------------------------------------------------------------------" << endl;
                    evaluationFile.close();
                    break;
                }

                case 2: {
                    //unutar grupe
                    std::vector<std::unique_ptr<FASTAformat>> fasta_objects;
                    //std::vector<std::unique_ptr<FASTQformat>> fastq_objects;
                    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAformat>(filePath1);
                    fasta_parser->parse(fasta_objects, -1);
                    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FASTQformat>(filePath2);
                    fastq_parser->parse(fastq_objects, -1);
                    //filtering by length
                    //std::vector<std::string> filtered_objects;
                    std::map<int, int> lengths;
                    for (const auto& it : fastq_objects) {
                        if (lengths.count(it->sequence.length()) > 0) {
                            lengths[it->sequence.length()]++;
                        } else {
                            lengths[it->sequence.length()] = 1;
                        }
                    }

                    auto max_length = std::max_element(lengths.begin(), lengths.end(),
                                                       [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
                                                           return p1.second < p2.second;
                                                       });
                    for (const auto& it : fastq_objects) {
                        if (it->sequence.length() <= max_length->first + length &&
                            it->sequence.length() >= max_length->first - length) filtered_objects.push_back(it->sequence);
                    }
                    cout << "Length of " << max_length->first << " nucleotide bases is most common (" << max_length->second << " appearances)" << endl;
                    cout << "Filtering by length of most common +/-"
                            << length << endl;


                    for(const auto& elem : fasta_objects){
                        std::cout << elem->name << std::endl;
                    }
                    map<string, map<int, int>> expectedDifsMap;
                    for (int j = 0; j < fasta_objects.size(); j++) {
                        map<int, int> newMap;
                        expectedDifsMap[fasta_objects[j]->sequence] = newMap;
                    }
                    for (int i = 0; i < filtered_objects.size(); i++) {
                        for (int j = 0; j < fasta_objects.size(); j++) {
                            //poravnaj i spremi u mapu mape
                            auto score = alignment::Align(filtered_objects[i].c_str(), filtered_objects[i].length(),
                                                          fasta_objects[j]->sequence.c_str(), fasta_objects[j]->sequence.length(), match, mismatch, gap,
                                                          alignment_al);
                            auto min = filtered_objects[i].length();
                            if (fasta_objects[j]->sequence.length() < filtered_objects[i].length()) min = fasta_objects[j]->sequence.length();
                            if (expectedDifsMap[fasta_objects[j]->sequence].find(min - score) == expectedDifsMap[fasta_objects[j]->sequence].end())
                                expectedDifsMap[fasta_objects[j]->sequence][min - score] = 1;
                            else
                                expectedDifsMap[fasta_objects[j]->sequence][min - score] += 1;
                        }
                    }
                    map<string, map<int, int>>::iterator it2 = expectedDifsMap.begin();
                    evaluationFile.open("../evaluation.txt", ios::out | ios::app);
                    evaluationFile << filePath1 << endl << filePath2 << endl;
                    evaluationFile << "Length +/- " << to_string(length) << " most common " << endl;
                    while (it2 != expectedDifsMap.end()) {
                        cout << it2->first << endl;
                        evaluationFile << "Alel (" << it2->first.size()
                                       << ", " << ")" << '\n';
                        evaluationFile << it2->first << '\n' << '\n';
                        map<int, int>::iterator it_inner = it2->second.begin();
                        evaluationFile << "numOfReps,diff" << endl;
                        cout << "diff,numOfReps" << endl;
                        while (it_inner != it2->second.end()) {
                            cout << it_inner->first << ',' << it_inner->second << endl;
                            evaluationFile << it_inner->first << ' ' << it_inner->second << '\n';
                            it_inner++;
                        }
                        it2++;
                    }
                    evaluationFile << "--------------------------------------------" << '\n';
                    evaluationFile.close();

                    break;
                }

                default: {
                    fprintf(stderr, "Invalid argument for evaluation type!\n");
                    cout << "Please choose one of the following:" << endl <<
                            "   1 for allele evaluation" << endl <<
                            "   2 for group evaluation" << endl;
                    invalidArgument = true;
                    break;
                }
            }

            if (invalidArgument) break;

            break;

        }

        case 10: {
            //ovo mogu bit 1. i 2.
            ofstream firstAlgorithmFile, secondAlgorithmFile,
                    thirdAlgorithmFile, kMeansAlgorithmFile, evaluationFile, msaFile;
            //cout << "Dva" << endl;
            string filePath = argv[1];
            int match = stoi(argv[2]), mismatch = stoi(argv[3]), gap = stoi(argv[4]);
            int length = stoi(argv[5]);
            int alignment_type = stoi(argv[6]);
            auto algorithm = stoi(argv[7]);
            int cluster_size = stoi(argv[8]);
            int addition_score = stoi(argv[9]);
            if (algorithm != 1) {
                if (algorithm != 2) {
                    invalidArgument = true;
                    fprintf(stderr, "Invalid type of algorithm for number of hyperparameters!\n");
                    cout << "Please choose the FIRST or SECOND algorithm (1/2) or USE DIFFERENT PARAMETERS and choose one of the following:\n" <<
                         "   3 for THIRD\n   4 for K-MEANS" << endl;
                }
            }
            if(invalidArgument) break;
            spoa::AlignmentType alignment_spoa;
            alignment::AlignmentType alignment_al;
            cout << filePath << endl;
            auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FASTQformat>(filePath);
            fastq_parser->parse(fastq_objects, -1);
            //filtering by length
            std::map<int, int> lengths;

            for (const auto& it : fastq_objects) {
                if (lengths.count(it->sequence.length()) > 0) {
                    lengths[it->sequence.length()]++;
                } else {
                    lengths[it->sequence.length()] = 1;
                }
            }

            auto max_length = std::max_element(lengths.begin(), lengths.end(),
                                               [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
                                                   return p1.second < p2.second;
                                               });
            for (const auto& it : fastq_objects) {
                if (it->sequence.length() <= max_length->first + length &&
                    it->sequence.length() >= max_length->first - length) filtered_objects.push_back(it->sequence);
            }
            cout << "Length of " << max_length->first << " nucleotide bases is most common (" << max_length->second << " appearances)" << endl;
            cout << "Filtering by length of most common +/-" << length << endl;

            switch (alignment_type) {
                case 1: {
                    alignment_spoa = spoa::AlignmentType::kNW;
                    alignment_al = alignment::AlignmentType::kGlobal;
                    break;
                }
                case 2: {
                    alignment_spoa = spoa::AlignmentType::kOV;
                    alignment_al = alignment::AlignmentType::kSemiGlobal;
                    break;
                }
                default: {
                    fprintf(stderr, "Invalid alignment type argument!\n");
                    cout << "Please choose one of the following:\n" <<
                            "   1 for global alignment\n" <<
                            "   2 for semi-global alignment\n" << endl;
                    invalidArgument = true;
                }
            }

            if (invalidArgument) break;

            auto alignment_engine = spoa::createAlignmentEngine(alignment_spoa,
                                                                1, -1, -1, -1);

            auto graph = spoa::createGraph();

            for (const auto& it: filtered_objects) {
                auto alignment = alignment_engine->align(it, graph);
                graph->add_alignment(alignment, it);
            }

            std::string consensus = graph->generate_consensus();
            vector<string> msa;
            graph->generate_multiple_sequence_alignment(msa);

            std::cout <<"Consensus (" << consensus.size() << ")" << std::endl;
            std::cout << consensus.c_str() << std::endl;

            switch (algorithm) {
                case 1: {
                    cout << "First algorithm:" << endl;
                    if (alignment_type == 1) {
                        cout << "global alignment, ";
                    } else {
                        cout << "semi-global alignment, ";
                    }
                    cout << "M = " << match << ", MM = " << mismatch << ", G = " << gap << endl;
                    // PRVI ALGORITAMMMMMM
                    vector<vector<string>> clusters;
                    vector<vector<string>> clusters_orig;
                    vector<string> vec;
                    vector<string> vec2;
                    string first_allel = msa.front();
                    vec.push_back(first_allel);
                    vec2.push_back(filtered_objects.front());
                    clusters.push_back(vec);
                    clusters_orig.push_back(vec2);

                    bool matched = false, this_cluster = true;

                    for (int i = 1; i < filtered_objects.size(); i++) {
                        auto new_obj = filtered_objects[i];
                        auto new_cons = msa[i];
                        for (int j = 0; j < clusters.size(); j++) {
                            auto cluster_vector = clusters[j];
                            for (auto it : cluster_vector) {
                                int distance = 0;
                                for (int i = 0; i < new_cons.size(); i++) {
                                    if (new_cons[i] != it[i]) distance++;
                                }
                                if (distance > addition_score) {
                                    this_cluster = false;
                                    break;
                                }
                            }
                            if (!this_cluster) {
                                this_cluster = true;
                                continue;
                            } else {
                                clusters[j].push_back(new_cons);
                                clusters_orig[j].push_back(new_obj);
                                //tu napraviti konsezusnu sekvencu
                                matched = true;
                                break;
                            }
                        }
                        if (!matched) {
                            vector<string> vec;
                            vector<string> vec2;
                            vec.push_back(new_cons);
                            vec2.push_back(new_obj);
                            clusters.push_back(vec);
                            clusters_orig.push_back(vec2);
                        }
                        matched = false;
                    }

                    cout << clusters_orig.size() << endl;
                    sort(clusters_orig.begin(), clusters_orig.end(),
                              [](const vector<string>& v1, const vector<string>& v2){
                                  return v1.size() > v2.size();
                              });
                    int the_size2 = 2, i2 = 0;
                    auto len2 = clusters_orig[i2].size();
                    firstAlgorithmFile.open("../first_algorithm.txt", ios::out | ios::app);
                    firstAlgorithmFile << "File " << filePath << endl;
                    firstAlgorithmFile << "Minimal size of clusters: " << cluster_size << '\n';
                    firstAlgorithmFile << "Parameters: " << '\n';
                    firstAlgorithmFile << "match = " << match << ", mismatch = "
                                       << mismatch << ", gap = " << gap <<
                                       ", score for making a new cluster = " << addition_score << '\n';
                    firstAlgorithmFile << "Number of clusters: " << clusters.size() << '\n';
                    firstAlgorithmFile << '\n';
                    cout << "First Algorithm Results: " << endl;
                    while(len2 >= cluster_size) {
                        auto alignment_engine_2 = spoa::createAlignmentEngine(alignment_spoa,
                                                                              match, mismatch, gap, gap);

                        auto graph_2 = spoa::createGraph();

                        for (const auto& it: clusters_orig[i2]) {
                            auto alignment_2 = alignment_engine_2->align(it, graph_2);
                            graph_2->add_alignment(alignment_2, it);
                        }

                        string consensus_2 = graph_2->generate_consensus();
                        cout << "Size of cluster: " << clusters_orig[i2].size() << endl;
                        cout <<"Allele (" << consensus_2.size() << ")" << endl;
                        cout << consensus_2.c_str() << endl;
                        cout << endl;
                        firstAlgorithmFile << "Size of cluster: " << clusters_orig[i2].size() << '\n';
                        firstAlgorithmFile << "Allele (" << consensus_2.size() << ")" << '\n';
                        firstAlgorithmFile << consensus_2.c_str() << '\n' << '\n';

                        i2++;
                        len2 = clusters_orig[i2].size();

                    }
                    firstAlgorithmFile << "-----------------------------------" << '\n';
                    firstAlgorithmFile.close();
                    break;

                }
                case 2: {
                    cout << "Second algorithm:" << endl;
                    if (alignment_type == 1) {
                        cout << "global alignment, ";
                    } else {
                        cout << "semi-global alignment, ";
                    }
                    cout << "M = " << match << ", MM = " << mismatch << ", G = " << gap << endl;
                    //DRUGI ALGORITAM
                    vector<vector<string>> clusters;
                    map<string, vector<string>> cluster_consensus_map;
                    int count = 1;
                    //int OK_score = 15;
                    vector<string> vec;
                    string first_allel = filtered_objects.back();
                    vec.push_back(first_allel);
                    clusters.push_back(vec);
                    cluster_consensus_map[first_allel] = vec;
                    bool matched = false, this_cluster = true;
                    string newKey;
                    cout << "------------------------" << endl;
                    vector<string> consensuses;
                    for (int i = 0; i < filtered_objects.size() - 1; i++) {
                        //std::cout << "tu smo " << std::endl;
                        auto new_obj = filtered_objects[i];
                        auto msa_obj = msa[i];
                        map<string, vector<string>>::iterator it = cluster_consensus_map.begin();
                        do {
                            auto cluster_consensus = it->first;
                            auto cluster_vector = it->second;
                            //std::cout << it->first << std::endl;
                            //for (auto it : cluster_vector) {
                            auto score = alignment::Align(cluster_consensus.c_str(), cluster_consensus.length(), new_obj.c_str(), new_obj.length(), match, mismatch, gap,
                                        alignment_al);
                            //std::cout << score << std::endl;
                            auto smaller = cluster_consensus.size();
                            if (new_obj.size() < smaller) smaller = new_obj.size();
                            auto theScore = smaller - score;
                            //cout << theScore << endl;
                            //cout << score << ", smaller = " << smaller<< endl;
                            if (theScore > addition_score) {
                                this_cluster = false;
                //                if (it == cluster_consensus_map.end()) break;
                //                it++;
                //                continue;
                            }
                            //}
                            if (this_cluster == false) {
                                if (it == cluster_consensus_map.end()) {
                                    this_cluster = true;
                                    break;
                                }
                                it++;
                            } else {
                                it->second.push_back(new_obj);
                                //tu napraviti konsezusnu sekvencu
                                matched = true;
                                auto alignment_engine_3 = spoa::createAlignmentEngine(alignment_spoa,
                                                                                      match, mismatch, gap, gap);

                                auto graph_3 = spoa::createGraph();

                                for (const auto& alel : it->second) {
                                    auto alignment_3 = alignment_engine_3->align(alel, graph_3);
                                    graph_3->add_alignment(alignment_3, alel);
                                }
                                newKey = graph_3->generate_consensus();
                                break;
                            }
                            this_cluster = true;
                        } while (it != cluster_consensus_map.end());
                        if (matched == false) {
                            std::vector<std::string> vec;
                            vec.push_back(new_obj);
                            cluster_consensus_map[new_obj] = vec;
                            //std::cout << new_obj << std::endl;
                            //clusters.push_back(vec);
                        } else {
                            std::vector<std::string> new_cluster_vector;
                            for (int k = 0; k < it->second.size(); k++)
                                new_cluster_vector.push_back(it->second[k]);
                            //std::copy(it->second.begin(), it->second.end(), std::back_inserter(new_cluster_vector));
                            //std::cout << newKey << " zamjenjuje " << it->first << std::endl;
                            //std::cout << it->second.size() << std::endl;
                            cluster_consensus_map.erase(it->first);
                            cluster_consensus_map[newKey] = new_cluster_vector;
                            //std::cout << newKey << ", duljina:  " << cluster_consensus_map[newKey].size() << std::endl;
                        }
                        matched = false;
                    }

                    //for (auto it: clusters) {
                    //    std::cout << it[0] << std::endl;
                    //}
                    std::cout << cluster_consensus_map.size() << std::endl;

                    int i = 0;
                    secondAlgorithmFile.open("../second_algorithm.txt", ios::out | ios::app);
                    secondAlgorithmFile << "File " << filePath << endl;
                    secondAlgorithmFile << "Minimal size of clusters: " << cluster_size << '\n';
                    secondAlgorithmFile << "Parameters: " << '\n';
                    secondAlgorithmFile << "match = " << match << ", mismatch = "
                    << mismatch << ", gap = " << gap <<
                    ", score for making a new cluster = " << addition_score << '\n';
                    //UPISI TIP PORAVNANJA
                    if (alignment_type == 1) {
                        secondAlgorithmFile << "Global alignment" << '\n';
                    } else {
                        secondAlgorithmFile << "Semi-global alignment" << '\n';
                    }

                    secondAlgorithmFile << "Number of clusters: " << cluster_consensus_map.size() << '\n' << '\n';

                    std::map<std::string, std::vector<std::string>>::iterator it = cluster_consensus_map.begin();
                    auto len = it->second.size();
                    do {
                        //std::cout << it->second.size() << std::endl;
                        if (it->second.size() < cluster_size) {
                            if (it != cluster_consensus_map.end()){
                                it++;
                                continue;
                            } else break;
                        }
                        //std::cout << '!' << std::endl;
                        cout << "Size of cluster: " << it->second.size() << endl;
                        cout <<"Consensus (" << it->first.size() << ")" << endl;
                        cout << it->first << endl;
                        consensuses.push_back(it->first);
                        secondAlgorithmFile << "Size of cluster: " << it->second.size() << '\n';
                        secondAlgorithmFile << "Consensus (" << it->first.size() <<
                        ")" << '\n';
                        secondAlgorithmFile << it->first << '\n' << '\n';
                        it++;
                        len = it->second.size();
                    } while(it != cluster_consensus_map.end());
                    secondAlgorithmFile << "----------------------------------" << '\n';
                    secondAlgorithmFile.close();

                    break;
                }
                default: {
                    fprintf(stderr, "Invalid algorithm choosen for number of hyperparameters (10)!\n");
                    cout << "Please choose one of the following:" << '\n' <<
                         "   1 for first\n   2 for second\n  OR choose 3 for third\n   4 for kmeans\n WITH different parameters" << endl;
                    invalidArgument = true;
                    break;
                }
            }
            if (invalidArgument) break;

            break;
        }
        case 12: {
            //treci algoritam
            ofstream thirdAlgorithmFile;
            string filePath = argv[1];
            int match = stoi(argv[2]), mismatch = stoi(argv[3]), gap = stoi(argv[4]);
            int length = stoi(argv[5]);
            int alignment_type = stoi(argv[6]);
            auto algorithm = stoi(argv[7]);
            int cluster_size = stoi(argv[8]);
            int newCentScore = stoi(argv[9]);
            int minDist = stoi(argv[10]);
            int mergingScore = stoi(argv[11]);
            if (algorithm != 3) {
                invalidArgument = true;
                fprintf(stderr, "Invalid type of algorithm for number of hyperparameters!\n");
                cout << "Please choose the THIRD algorithm (3) or USE DIFFERENT PARAMETERS and choose one of the following:\n" <<
                     "   1 for FIRST\n   2 for SECOND\n   4 for K-MEANS" << endl;
            }
            if (invalidArgument) break;
            spoa::AlignmentType alignment_spoa;
            alignment::AlignmentType alignment_al;
            cout << filePath << endl;
            auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FASTQformat>(filePath);
            fastq_parser->parse(fastq_objects, -1);
            //filtering by length
            std::map<int, int> lengths;
            for (const auto& it : fastq_objects) {
                if (lengths.count(it->sequence.length()) > 0) {
                    lengths[it->sequence.length()]++;
                } else {
                    lengths[it->sequence.length()] = 1;
                }
            }

            auto max_length = std::max_element(lengths.begin(), lengths.end(),
                                               [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
                                                   return p1.second < p2.second;
                                               });
            for (const auto& it : fastq_objects) {
                if (it->sequence.length() <= max_length->first + length &&
                    it->sequence.length() >= max_length->first - length) filtered_objects.push_back(it->sequence);
            }
            cout << "Length of " << max_length->first << " nucleotide bases is most common (" << max_length->second << " appearances)" << endl;

            cout << "Filtering by length of most common +/-" << length << endl;
            switch (alignment_type) {
                case 1: {
                    alignment_spoa = spoa::AlignmentType::kNW;
                    alignment_al = alignment::AlignmentType::kGlobal;
                    break;
                }
                case 2: {
                    alignment_spoa = spoa::AlignmentType::kOV;
                    alignment_al = alignment::AlignmentType::kSemiGlobal;
                    break;
                }
                default: {
                    fprintf(stderr, "Invalid alignment type argument!\n");
                    cout << "Please choose one of the following:\n" <<
                         "   1 for global alignment\n" <<
                         "   2 for semi-global alignment\n" << endl;
                    invalidArgument = true;
                }
            }

            if (invalidArgument) break;

            auto alignment_engine = spoa::createAlignmentEngine(alignment_spoa,
                                                                match, mismatch, gap, gap);

            auto graph = spoa::createGraph();

            for (const auto& it: filtered_objects) {
                auto alignment = alignment_engine->align(it, graph);
                graph->add_alignment(alignment, it);
            }

            std::string consensus = graph->generate_consensus();
            vector<string> msa;
            graph->generate_multiple_sequence_alignment(msa);

            switch (algorithm) {
                case 3: {
                    ofstream thirdAlgorithmFile;
                    thirdAlgorithmFile.open("../third_algorithm.txt", ios::out | ios::app);
                    cout << "Third algorithm:" << endl;
                    if (alignment_type == 1) {
                        cout << "global alignment, ";
                    } else {
                        cout << "semi-global alignment, ";
                    }
                    cout << "M = " << match << ", MM = " << mismatch << ", G = " << gap << endl;

                    map<int, string> centroids;
                    int cluster_index = 0;
                    centroids[cluster_index] = *msa.begin();
                    vector<string>::iterator it_s = msa.begin();
                    while (it_s != msa.end()) {
                        bool cluster_not_found = true;
                        string msa_current;
                        if (it_s->c_str() != nullptr) msa_current = (*it_s);
                        else continue;
                        map<int, string>::iterator it_c = centroids.begin();
                        while (it_c != centroids.end()) {
                            string centroid;
                            if (!it_c->second.empty()) centroid = it_c->second;
                            else continue;
                            int distance = 0;
                            for (int i = 0; i < msa_current.size(); i++) {
                                if (msa_current[i] != centroid[i]) distance++;
                            }
                            if (distance < newCentScore) {
                                //cout << distance << endl;
                                cluster_not_found = false;
                            }
                            it_c++;
                        }

                        if (cluster_not_found) {
                            cluster_index += 1;
                            centroids[cluster_index] = msa_current;
                        }
                        it_s++;
                    }

                    map<int, string>::iterator it_c = centroids.begin();
                    map<int, vector<string>> clusters;
                    while (it_c != centroids.end()) {
                        int index = it_c->first;

                        string centroid = it_c->second;
                        vector<string>::iterator it_m = msa.begin();
                        while(it_m != msa.end()) {
                            string msa_current = (*it_m);
                            int distance = 0;
                            for (int i = 0; i < msa_current.size(); i++) {
                                if (msa_current[i] != centroid[i]) distance++;
                            }
                            if (distance < minDist) {
                                clusters[index].push_back(msa_current);
                            }
                            it_m++;
                        }
                        it_c++;
                    }

                    vector<string> consensuses;
                    //cout << "ima " << clusters.size() << " klastera" << endl;
                    map<int, vector<string>>::iterator it = clusters.begin();
                    while (it != clusters.end()) {
                        vector<string> cluster_sequences = it->second;
                        auto alignment_engine_c = spoa::createAlignmentEngine(alignment_spoa, match, mismatch, gap, gap);
                        auto graph_c = spoa::createGraph();

                        for (const auto &it_seq : cluster_sequences)
                        {
                            //cout << it_c << endl;
                            auto alignment_c = alignment_engine_c->align(it_seq, graph_c);
                            graph_c->add_alignment(alignment_c, it_seq);
                        }

                        string consensus = graph_c->generate_consensus();
                        consensuses.push_back(consensus);

                        consensus.erase(remove(consensus.begin(), consensus.end(), '-'), consensus.end());
                        it++;
                    }

                    map<string, vector<string>> cluster_map;
                    for (int i = 0; i < consensuses.size(); ++i) {
                        cluster_map[consensuses[i]] = clusters.find(i)->second;
                    }
                    map<int, vector<string>> merged_centroids;
                    int index = 0;
                    merged_centroids[index].push_back(*(consensuses.begin()));

                    //za svaki konsenzus osim prvog koji je vec dodan
                    vector<string>::iterator it_m = consensuses.begin() + 1;
                    while (it_m != consensuses.end()) {
                        bool cluster_found = false;
                        string current = (*it_m);
                        map<int, vector<string>>::iterator it2 = merged_centroids.begin();
                        while (it2 != merged_centroids.end()) {
                            int max_k = 0;
                            vector<string>::iterator it3 = it2->second.begin();
                            while (it3 != it2->second.end()) {
                                string seq = (*it3).substr(0);
                                seq.erase(remove(seq.begin(), seq.end(), '-'), seq.end());
                                string curSeq = current.substr(0);
                                curSeq.erase(remove(curSeq.begin(), curSeq.end(), '-'), curSeq.end());

                                auto smaller = seq.length();
                                if (curSeq.length() < smaller) smaller = curSeq.length();
                                auto dist = alignment::Align(seq.c_str(), seq.length(), curSeq.c_str(), curSeq.length(), match, mismatch, gap,
                                                             alignment_al);

                                if (smaller - dist > max_k) {
                                    max_k = smaller - dist;
                                }
                                it3++;
                            }
                            if (max_k < mergingScore) {
                                //cout << "manje" << endl;
                                //cout << current << endl;
                                (it2->second).push_back(current);
                                cluster_found = true;
                                break;
                            }
                            it2++;
                        }
                        if (!(cluster_found)) {
                            //cout << "nije manje" << endl;
                            index += 1;
                            merged_centroids[index].push_back(current);
                        }
                        it_m++;
                    }

                    //grupe
                    map<int, vector<string>> merged_clusters;
                    for (int i = 0; i < merged_centroids.size(); i++) {
                        //u cents je cluster tj sve sekvence
                        vector<string> cents = merged_centroids.find(i)->second;
                        for (int j = 0; j < cents.size(); j++) {
                            vector<string> cluster = cluster_map.find(cents[j])->second;

                            if (j == 0) {
                                for(int k = 0; k < cluster.size(); k++) {
                                    merged_clusters[i].push_back(cluster[k]);
                                }
                            } else {
                                for(int k = 0; k < cluster.size(); k++) {
                                    if (find(merged_clusters[i].begin(),
                                             merged_clusters[i].end(), cluster[k]) ==
                                        merged_clusters[i].end())
                                        merged_clusters[i].push_back(cluster[k]);
                                }
                            }

                        }
                    }
                    cout << "ima " << merged_clusters.size() << " klastera" << endl;
                    map<int, vector<string>>::iterator it_merged = merged_clusters.begin();
                    while (it_merged != merged_clusters.end()) {
                        //cout << "jedan" << endl;
                        if (it_merged->second.size() < cluster_size) {
                            //cout << "brisem " << it->second.size() << endl;
                            merged_clusters.erase(it_merged++);
                        } else {
                            //cout << it->second.size() << " ne brisem" << endl;
                            it_merged++;
                        }
                    }

                    vector<string> alelles;
                    vector<int> sizes;
                    it_merged = merged_clusters.begin();
                    while (it_merged != merged_clusters.end()) {
                        auto alignment_engine = spoa::createAlignmentEngine(alignment_spoa, match, mismatch, gap, gap);

                        auto graph = spoa::createGraph();

                        vector<string> cluster_sequences = it_merged->second;
                        //cout << "klaster broj " << it->first << " koji ima " <<
                        //     it->second.size() << " elemenata" << endl;
                        for (const auto &it2 : cluster_sequences) {
                            auto alignment = alignment_engine->align(it2, graph);
                            graph->add_alignment(alignment, it2);
                        }
                        string consensus = graph->generate_consensus();
                        alelles.push_back(consensus);
                        sizes.push_back(it_merged->second.size());
                        it_merged++;
                    }
                    for (auto it = alelles.begin(); it != alelles.end(); it++) {
                        (*it).erase(remove((*it).begin(), (*it).end(), '-'), (*it).end());
                    }
                    thirdAlgorithmFile << "File " << filePath << '\n';
                    thirdAlgorithmFile << "Length most common +/- " << length << '\n';
                    thirdAlgorithmFile << "Alignment type: " << alignment_al <<
                                       "match = " << match << ", mismatch = " << mismatch << ", gap = " << gap << '\n';

                    thirdAlgorithmFile << "Score for new centroid: " << newCentScore << '\n' <<
                                       "Score for query addition to cluster: " << minDist - 1 << " (<" << minDist << ")" << '\n' <<
                                       "Score for merging clusters: " << mergingScore - 1 << " (<" << mergingScore << ")" << '\n' <<
                                       "Minimal size of clusters: " << cluster_size << '\n';
                    thirdAlgorithmFile << '\n' << "Number of clusters: " << alelles.size() << '\n' << "-------" << '\n';
                    for (int i = 0; i < alelles.size(); i++) {
                        cout << endl << endl;
                        thirdAlgorithmFile << "Size of cluster: " << sizes[i] << '\n';
                        thirdAlgorithmFile << "Allele (" << alelles[i].size() << ")" << '\n';
                        thirdAlgorithmFile << alelles[i].c_str() << '\n' << '\n';
//                        if (all_alleles.find(alelles[i]) == all_alleles.end()) {
//                            all_alleles.insert(alelles[i]);
//                            alleles_num[alelles[i]] = 1;
//                            numOfDistinct++;
//                        } else {
//                            alleles_num[alelles[i]]++;
//                        }
//                        numOfAlleles++;
                        cout <<"Size of cluster: " << sizes[i] << endl;
                        cout <<"Allele (" << alelles[i].size() << ")";
                        cout << alelles[i].c_str() << endl;
                    }
                    thirdAlgorithmFile << "----------------------------------------------------------------------------" << '\n';
//                    vector<string>().swap(alelles);
//                    vector<int>().swap(sizes);
//                    std::vector<std::unique_ptr<FASTQformat>>().swap(fastq_objects_all);
//                    vector<string>().swap(filtered_objects_all);
//                    vector<string>().swap(consensuses_shorter);
//                    map<int, string>().swap(centroids_shorter);
//                    map<int, vector<string>>().swap(clusters_shorter);
//                    map<int, vector<string>>().swap(merged_clusters);
//                    map<int, vector<string>>().swap(merged_centroids);

//                    thirdAlgorithmFile << '\n' << '\n';
//                    thirdAlgorithmFile << "Number of alleles: " << numOfAlleles << '\n';
//                    thirdAlgorithmFile << "Number of distinct alleles: " << numOfDistinct << '\n';
//                    for (auto elem: all_alleles) {
//                        thirdAlgorithmFile << "Allele (" << elem.size() << ") found in " << alleles_num[elem] << " samples." << '\n';
//                        thirdAlgorithmFile << elem << '\n';
//                    }
                    thirdAlgorithmFile.close();

                    break;
                }
                default: {
                    invalidArgument = true;
                    fprintf(stderr, "Invalid type of algorithm for number of hyperparameters!\n");
                    cout << "Please choose the THIRD algorithm (3) or USE DIFFERENT PARAMETERS and choose one of the following:\n" <<
                    "   1 for FIRST\n   2 for SECOND\n   4 for K-MEANS" << endl;
                }
            }
            break;

        }

        case 13: {
            //k means
            ofstream thirdAlgorithmFile;
            cout << "13" << endl;
            string filePath = argv[1];
            int match = stoi(argv[2]), mismatch = stoi(argv[3]), gap = stoi(argv[4]);
            int length = stoi(argv[5]);
            int alignment_type = stoi(argv[6]);
            auto algorithm = stoi(argv[7]);
            int cluster_size = stoi(argv[8]);
            int initParam = stoi(argv[9]);
            int newClusterScore = stoi(argv[10]);
            int numOfIterations = stoi(argv[11]);
            int k = stoi(argv[12]);
            spoa::AlignmentType alignment_spoa;
            alignment::AlignmentType alignment_al;
            cout << filePath << endl;
            auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FASTQformat>(filePath);
            fastq_parser->parse(fastq_objects, -1);
            //filtering by length
            std::map<int, int> lengths;
            for (const auto& it : fastq_objects) {
                if (lengths.count(it->sequence.length()) > 0) {
                    lengths[it->sequence.length()]++;
                } else {
                    lengths[it->sequence.length()] = 1;
                }
            }

            auto max_length = std::max_element(lengths.begin(), lengths.end(),
                                               [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
                                                   return p1.second < p2.second;
                                               });
            for (const auto& it : fastq_objects) {
                if (it->sequence.length() <= max_length->first + length &&
                    it->sequence.length() >= max_length->first - length) filtered_objects.push_back(it->sequence);
            }
            cout << "Length of " << max_length->first << " nucleotide bases is most common (" << max_length->second << " appearances)" << endl;

            cout << "Filtering by length of most common +/-" << length << endl;
            switch (alignment_type) {
                case 1: {
                    alignment_spoa = spoa::AlignmentType::kNW;
                    alignment_al = alignment::AlignmentType::kGlobal;
                    break;
                }
                case 2: {
                    alignment_spoa = spoa::AlignmentType::kOV;
                    alignment_al = alignment::AlignmentType::kSemiGlobal;
                    break;
                }
                default: {
                    fprintf(stderr, "Invalid alignment type argument!\n");
                    cout << "Please choose one of the following:\n" <<
                         "   1 for global alignment\n" <<
                         "   2 for semi-global alignment\n" << endl;
                    invalidArgument = true;
                }
            }

            if (invalidArgument) break;

            switch (algorithm) {
                case 4: {
                    //k-means
                    ofstream kMeansAlgorithmFile;
                    kMeansAlgorithmFile.open("../kMeans.txt", ios::out | ios::app);
                    //int numOfClusters = 3, numOfIterations = 20;
                    //int match_k = 1, mismatch_k = -1, gap_k = -1;
                    //int newClusterScore = 30, clusterSize = 3;
                    switch (initParam) {
                        case 0: {
                            //fair
                            std::cout << "K-means algorithm FAIR: " << std::endl;
                            auto kMeansClusters = kMeans::generateClusters(filtered_objects, k, numOfIterations, match, mismatch, gap,
                                                                            kMeans::CentroidInit::kFair, newClusterScore, cluster_size);


                            kMeansAlgorithmFile << "K-means algorithm fair: " << '\n';
                            kMeansAlgorithmFile << "File " << filePath << '\n';
                            kMeansAlgorithmFile << "Score for new cluster: " << newClusterScore <<
                                                ", minimal size of clusters: " << cluster_size << '\n';
                            if (alignment_type == 2) kMeansAlgorithmFile << "Semi-global, ";
                            else kMeansAlgorithmFile << "Global, ";
                            kMeansAlgorithmFile << "match = " << match <<
                                                ", mismatch = " << mismatch << ", gap = " << gap << '\n';
                            for (auto it : kMeansClusters){
                                std::cout <<"Consensus (" << it.size() << ")" << std::endl;
                                std::cout << it << std::endl;
                                std::cout << std::endl;
                                kMeansAlgorithmFile << "Allele (" << it.size()
                                                    << ")" << '\n';
                                kMeansAlgorithmFile << it << '\n' << '\n';
                            }

                            break;
                        }

                        case 1: {
                            //worst case
                            std::cout << "K-means algorithm WORST CASE: " << std::endl;
                            auto kMeansClusters = kMeans::generateClusters(filtered_objects, k, numOfIterations, match, mismatch, gap,
                                                                            kMeans::CentroidInit::kWorst, newClusterScore, cluster_size);


                            kMeansAlgorithmFile << "K-means algorithm worst case: " << '\n';
                            kMeansAlgorithmFile << "File " << filePath << '\n';
                            kMeansAlgorithmFile << "Score for new cluster: " << newClusterScore <<
                            ", minimal size of clusters: " << cluster_size << '\n';
                            if (alignment_type == 2) kMeansAlgorithmFile << "Semi-global, ";
                            else kMeansAlgorithmFile << "Global, ";
                            kMeansAlgorithmFile << "match = " << match <<
                                                     ", mismatch = " << mismatch << ", gap = " << gap << '\n';
                            for (auto it : kMeansClusters){
                                std::cout <<"Consensus (" << it.size() << ")" << std::endl;
                                std::cout << it << std::endl;
                                std::cout << std::endl;
                                kMeansAlgorithmFile << "Allele (" << it.size()
                                                    << ")" << '\n';
                                kMeansAlgorithmFile << it << '\n' << '\n';
                            }

                            break;
                        }

                        default: {
                            fprintf(stderr, "Invalid number for initialization parameter!\n");
                            cout << "Please choose one of the following:" << endl <<
                            "   0 for FAIR initialization" << endl <<
                            "   1 for WORST CASE initialization" << endl;
                            invalidArgument = true;
                            kMeansAlgorithmFile.close();
                            break;
                        }
                    }

                    if (invalidArgument) break;

//
//                    kMeansAlgorithmFile << "results: " << '\n';
//                    for (int i = 0; i < fasta_objects.size(); i++) {
//                        auto maxScore = 0;
//                        string theSequence;
//                        for (auto seq : kMeansClusters3) {
//                            auto score = alignment::Align(fasta_objects[i]->sequence.c_str(), fasta_objects[i]->sequence.length(),
//                                                          seq.c_str(), seq.length(), 1, 0, -1,
//                                                          alignment::AlignmentType::kSemiGlobal);
//                            if (score > maxScore) {
//                                maxScore = score;
//                                theSequence = seq;
//                            }
//
//                        }
//                        kMeansAlgorithmFile << fasta_objects[i]->sequence << ' ' << fasta_objects[i]->sequence.length() << endl;
//                        kMeansAlgorithmFile << theSequence << ' ' << theSequence.length() << endl;
//                        if (theSequence.length() < fasta_objects[i]->sequence.length()) kMeansAlgorithmFile << theSequence.length() - maxScore << endl;
//                        else kMeansAlgorithmFile << fasta_objects[i]->sequence.length() - maxScore << endl;
//                    }
//                    kMeansAlgorithmFile << "--------------------------------------------" << '\n';
//

//                    std::vector<std::string> cons;
//                    cons.push_back("GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA");
//                    cons.push_back("GATCCTCTCTCTGCAGCACATTTCCTGCTGTATACTACGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTACGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGTCCGCCAAGTACTGGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACAGGTACTGCAGACACAACTACGGGGTTCTTGACAGTTTCGCTGGTGCAGCGGTCGAGGTGACGCGAAA");
//                    cons.push_back("GATCCTCTCTCTGCAGCACATTTCCTGGAGCATCATAAGTGCGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCAGTTCCTGCAGAGATACATCTATAACCGGGAAGAGTACGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGCCGTCCGCCAAGTACTATAACAGCCAGAAGGAGCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAGGTACTGCAGACACAACTACGGGGTCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA");
//                    auto kMeansClusters2 = kMeans::generateClusters(filtered_objects_longer, numOfClusters, numOfIterations, cons, match_k, mismatch_k, gap_k);
//                    std::cout << "K-means algorithm with given consensuses: " << std::endl;
//                    //kMeansAlgorithmFile.open("../kMeans.txt", ios::out | ios::app);
//                    kMeansAlgorithmFile << "File J29_B_CE_IonXpress_005" << '\n';
//                    kMeansAlgorithmFile << "with given consensuses" << '\n';
//                    kMeansAlgorithmFile << "Semi-global, m = " << match_k <<
//                    ", mm = " << mismatch_k << ", gap = " << gap_k << '\n';
//                    for (auto it : kMeansClusters2){
//                        std::cout <<"Consensus (" << it.size() << ")" << std::endl;
//                        std::cout << it << std::endl;
//                        std::cout << std::endl;
//                        kMeansAlgorithmFile << "Allele (" << it.size()
//                                            << ")" << '\n';
//                        kMeansAlgorithmFile << it << '\n' << '\n';
//                    }
//                    kMeansAlgorithmFile << "results: " << '\n';
//                    for (int i = 0; i < fasta_objects.size(); i++) {
//                        auto maxScore = 0;
//                        string theSequence;
//                        for (auto seq : kMeansClusters2) {
//                            auto score = alignment::Align(fasta_objects[i]->sequence.c_str(), fasta_objects[i]->sequence.length(),
//                                                          seq.c_str(), seq.length(), 1, 0, -1,
//                                                          alignment::AlignmentType::kSemiGlobal);
//                            if (score > maxScore) {
//                                maxScore = score;
//                                theSequence = seq;
//                            }
//
//                        }
//                        kMeansAlgorithmFile << fasta_objects[i]->sequence << ' ' << fasta_objects[i]->sequence.length() << endl;
//                        kMeansAlgorithmFile << theSequence << ' ' << theSequence.length() << endl;
//                        if (theSequence.length() < fasta_objects[i]->sequence.length()) kMeansAlgorithmFile << theSequence.length() - maxScore << endl;
//                        else kMeansAlgorithmFile << fasta_objects[i]->sequence.length() - maxScore << endl;
//                    }
//                    kMeansAlgorithmFile << "--------------------------------------------" << '\n';
//
//                    std::vector<std::string> cons_30;
//                    cons_30.push_back("GATCCTCTCTCTGCAGCACATTTCCTGGAGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGCGGTTCCTGGACAGATACTTCTATAACCGGGAAGAGTACGTGCGCTTCGACAGCGACTGGGGCGAGTTCCGGGCGGTGACCGAGCTGGGGCGGCCGTCCGCCAAGTACTGGAACAGCCAGAAGGATTTCATGGAGCAGAAGCGGGCCGAGGTGGACACGGTGTGCAGACACAACTACGGGGTTATTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA");
//                    cons_30.push_back("GATCCTCTCTCTGCAGCACATTTCCTGGAGCATCTTAAGGCCGAGTGTCATTTCTTCAACGGGACGGAGCGGATGCAGTTCCTGGCGAGATACTTCTATAACGGAGAAGAGTACGCGCGCTTCGACAGCGACGTGGGCGAGTTCCGGGCGGTGACCGAGCTGGGGCGGCCGGACGCCAAGTACTGGAACAGCCAGAAGGAGATCCTGGAGCAGCACGGGGCAGAGGTGGACAGGTACTGCAGACACAACTACGGGGTCGGTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA");
//                    cons_30.push_back("GATCCTCTCTCTGCAGCACATTTCCTGATGTATACTAAGAAAGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGCTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTTCCGGGCGGTGACCGAGCTGGGGCGGCCGGACGCCGAGGCTGGAACAGACAGAAGGAGCTCCTGGAGCAGAGGCGGGCCGCGGTGGACACGTACTGCAGACACAACTACGGGGTTATTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA");
//                    auto kMeansClusters = kMeans::generateClusters(filtered_objects_longer_30, numOfClusters, numOfIterations, cons_30, match_k, mismatch_k, gap_k);
//
//                    std::cout << "K-means algorithm: " << std::endl;
//                    kMeansAlgorithmFile << "File File J30_B_CE_IonXpress_006" << '\n';
//                    kMeansAlgorithmFile << "with given consensuses" << '\n';
//                    kMeansAlgorithmFile << "Semi-global, m = " << match_k <<
//                                        ", mm = " << mismatch_k << ", gap = " << gap_k << '\n';
//                    for (auto it : kMeansClusters){
//                        std::cout <<"Consensus (" << it.size() << ")" << std::endl;
//                        std::cout << it << std::endl;
//                        std::cout << std::endl;
//                        kMeansAlgorithmFile << "Allele (" << it.size()
//                        << ")" << '\n';
//                        kMeansAlgorithmFile << it << '\n' << '\n';
//                    }
//                    kMeansAlgorithmFile << "results: " << '\n';
//                    for (int i = 0; i < fasta_objects_30.size(); i++) {
//                        auto maxScore = 0;
//                        string theSequence;
//                        for (auto seq : kMeansClusters) {
//                            auto score = alignment::Align(fasta_objects_30[i]->sequence.c_str(), fasta_objects_30[i]->sequence.length(),
//                                                          seq.c_str(), seq.length(), 1, 0, -1,
//                                                          alignment::AlignmentType::kSemiGlobal);
//                            if (score > maxScore) {
//                                maxScore = score;
//                                theSequence = seq;
//                            }
//
//                        }
//                        kMeansAlgorithmFile << fasta_objects[i]->sequence << ' ' << fasta_objects[i]->sequence.length() << endl;
//                        kMeansAlgorithmFile << theSequence << ' ' << theSequence.length() << endl;
//                        if (theSequence.length() < fasta_objects[i]->sequence.length()) kMeansAlgorithmFile << theSequence.length() - maxScore << endl;
//                        else kMeansAlgorithmFile << fasta_objects[i]->sequence.length() - maxScore << endl;
//                    }
                    kMeansAlgorithmFile << "--------------------------------------------" << '\n';
                    kMeansAlgorithmFile.close();

                    break;
                }

                default: {
                    fprintf(stderr, "Invalid number of arguments for chosen algorithm!\n");
                    cout << "Please choose K-MEANS algorithm (4) or USE DIFFERENT PARAMETERS and choose one of the following:\n" <<
                         "   1 for FIRST\n   2 for SECOND\n   3 for THIRD" << endl;
                    invalidArgument = true;

                    break;
                }
            }
            break;
        }

        default: {
            fprintf(stderr, "Invalid number of arguments passed to program!\n");
            invalidArgument = true;
            break;
        }

        if (invalidArgument) break;

    }
    //imena fileova bez patha
//    vector<string> fileNames = {"J1_S_CE_IonXpress_018.fastq", "J2_S_CE_IonXpress_019.fastq",
//                                "J3_S_CE_IonXpress_020.fastq", "J4_GK_CE_IonXpress_051.fastq",
//                                "J4_S_CE_IonXpress_021.fastq", "J5_GK_CE_IonXpress_052.fastq",
//                                "J5_S_CE_IonXpress_022.fastq", "J6_GK_CE_IonXpress_001.fastq",
//                                "J6_S_CE_IonXpress_023.fastq", "J7_GK_CE_IonXpress_002.fastq",
//                                "J7_S_CE_IonXpress_024.fastq", "J8_GK_CE_IonXpress_003.fastq",
//                                "J8_S_CE_IonXpress_025.fastq", "J9_GK_CE_IonXpress_004.fastq",
//                                "J9_L_CE_IonXpress_029.fastq", "J9_S_CE_IonXpress_026.fastq",
//                                "J10_L_CE_IonXpress_030.fastq", "J10_S_CE_IonXpress_027.fastq",
//                                "J11_L_CE_IonXpress_031.fastq", "J11_S_CE_IonXpress_028.fastq",
//                                "J12_L_CE_IonXpress_032.fastq", "J13_L_CE_IonXpress_033.fastq",
//                                "J14_L_CE_IonXpress_034.fastq", "J15_L_CE_IonXpress_035.fastq",
//                                "J16_L_CE_IonXpress_036.fastq", "J17_L_CE_IonXpress_037.fastq",
//                                "J18_L_CE_IonXpress_038.fastq", "J19_L_CE_IonXpress_039.fastq",
//                                "J20_L_CE_IonXpress_040.fastq", "J21_L_CE_IonXpress_041.fastq",
//                                "J22_L_CE_IonXpress_042.fastq", "J23_L_CE_IonXpress_043.fastq",
//                                "J25_L_CE_IonXpress_045.fastq", "J27_L_CE_IonXpress_046.fastq",
//                                "J28_L_CE_IonXpress_015.fastq", "J29_B_CE_IonXpress_005.fastq",
//                                "J30_B_CE_IonXpress_006.fastq", "J31_B_CE_IonXpress_007.fastq",
//                                "J32_B_CE_IonXpress_008.fastq", "J33_B_CE_IonXpress_009.fastq",
//                                "J34_B_CE_IonXpress_010.fastq", "J35_B_CE_IonXpress_011.fastq",
//                                "J36_B_CE_IonXpress_012.fastq", "J37_B_CE_IonXpress_013.fastq",
//                                "J38_B_CE_IonXpress_014.fastq", "J40_B_CE_IonXpress_016.fastq",
//                                "J41_B_CE_IonXpress_017.fastq"};
//
//    set<string> all_alleles;
//    map<string, int> alleles_num;
//    int numOfAlleles = 0, numOfDistinct = 0;
//    thirdAlgorithmFile.open("../third_algorithm.txt", ios::out | ios::app);
//
//    for (int b = 0; b < fileNames.size(); b++){
//        string pathName = "../Bioinformatika-jeleni/fastq/";
//        pathName += fileNames[b];
//        std::vector<std::unique_ptr<FASTQformat>> fastq_objects_all;
//        auto fastq_parser_all = bioparser::createParser<bioparser::FastqParser, FASTQformat>(pathName);
//        fastq_parser_all->parse(fastq_objects_all, -1);
//        vector<string> filtered_objects_all;
//        for (const auto& it : fastq_objects_all) {
//            if (it->sequence.length() <= 296 + 5
//                && it->sequence.length() >= 296 - 5) filtered_objects_all.push_back(it->sequence);
//        }
//        if (filtered_objects_all.empty()) continue;
//        auto alignment_engine_all = spoa::createAlignmentEngine(spoa::AlignmentType::kNW,
//                                                                match, mismatch, gap, gap);
//
//        auto graph_all = spoa::createGraph();
//
//        for (const auto& it: filtered_objects_all) {
//            auto alignment_all = alignment_engine_all->align(it, graph_all);
//            graph_all->add_alignment(alignment_all, it);
//        }
//
//        std::string consensus_all = graph_all->generate_consensus();
//
//        vector<string> msa_all;
//        graph_all->generate_multiple_sequence_alignment(msa_all);
//
//        cout <<"File " << fileNames[b] << ", consensus +/-5 length (" << consensus_all.size() << ")" << endl;
//        cout << consensus_all.c_str() << endl;
//
//

    return 0;
}