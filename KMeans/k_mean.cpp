//
// Created by nera on 11. 05. 2020..
//

#include <limits>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <map>
#include <set>
#include <functional>
#include "k_mean.h"
#include "../spoa/include/spoa/spoa.hpp"
#include "../alignment/alignment.h"

//parameters random, fair, the worst
//if not fair, newClusterScore is ignored
namespace kMeans {

    std::vector<std::string> generateClusters(std::vector<std::string> sequences,
            int k, int numberOfIterations, int match, int mismatch, int gap,
            CentroidInit initParam, int newClusterScore, int clusterSize) {
        bool converged = true;
        int iterations = 0;
        std::vector<kMeans::Cluster> clusters;
        std::map<int, std::vector<std::string>> clusterConsensusesMap;

        if (initParam == CentroidInit::kRandom) {
            //assigning random k centroids
        } else if (initParam == CentroidInit::kWorst) {
            //determining the consensus sequence and
            //  choosing the sequences that are most different from the consensus
            //  for initial clusters' consensuses

            auto alignment_engine = spoa::createAlignmentEngine(spoa::AlignmentType::kOV,
                                                                match, mismatch, gap, gap);

            auto graph = spoa::createGraph();
            for (const auto& it: sequences) {
                auto alignment = alignment_engine->align(it, graph);
                graph->add_alignment(alignment, it);
            }

            std::string consensus = graph->generate_consensus();

            std::map<std::string, int> consensusAlignmentMap;

            for (const auto& it : sequences) {
                if (consensusAlignmentMap.count(it) == 0) {
                    auto score = alignment::Align(it.c_str(), it.length(),
                                                  consensus.c_str(), consensus.length(),
                                                  match, mismatch, gap, alignment::AlignmentType::kSemiGlobal);

                    consensusAlignmentMap[it] = score;
                } else {
                    continue;
                }
            }

            typedef std::function<bool(std::pair<std::string, int>, std::pair<std::string, int>)> Comparator;
            Comparator compFunctor =
                    [](std::pair<std::string, int> elem1 ,std::pair<std::string, int> elem2)
                    {
                        return elem1.second < elem2.second;
                    };

            std::set<std::pair<std::string, int>, Comparator> setOfCons(
                    consensusAlignmentMap.begin(), consensusAlignmentMap.end(), compFunctor);

            int i = 0;
            for (std::pair<std::string, int> el : setOfCons) {
                std::vector<std::string> cluster_vector;
                kMeans::Cluster newCluster;
                newCluster.centroid = el.first;
                newCluster.cluster_sequences = cluster_vector;
                clusters.push_back(newCluster);
                i++;
                if (i == k) break;
            }

            for (int i = 0; i < clusters.size(); i++) {
                std::vector<std::string> consensusesVector;
                clusterConsensusesMap[i] = consensusesVector;
            }
        } else if (initParam == CentroidInit::kFair) {
            //fairly choosing initial centroids
            //parameter k is ignored

            auto alignment_engine = spoa::createAlignmentEngine(spoa::AlignmentType::kOV,
                                                                match, mismatch, gap, gap);

            auto graph = spoa::createGraph();
            for (const auto& it: sequences) {
                auto alignment = alignment_engine->align(it, graph);
                graph->add_alignment(alignment, it);
            }

            std::string consensus = graph->generate_consensus();

            for (const auto& it: sequences) {
                bool newCent = true;
                for (auto it_c: clusters) {
                    auto score = alignment::Align(it.c_str(), it.size(),
                            it_c.centroid.c_str(), it_c.centroid.size(), match,
                            mismatch, gap, alignment::AlignmentType::kSemiGlobal);
                    auto smaller = it.size();
                    if (it_c.centroid.size() < smaller) smaller = it_c.centroid.size();
                    auto totalScore = smaller - score;
                    if (totalScore < newClusterScore) {
                        newCent = false;
                    }
                }
                if (newCent) {
                    std::vector<std::string> cluster_vector;
                    kMeans::Cluster newCluster;
                    newCluster.centroid = it;
                    newCluster.cluster_sequences = cluster_vector;
                    clusters.push_back(newCluster);
                }
            }

            for (int i = 0; i < clusters.size(); i++) {
                std::vector<std::string> consensusesVector;
                clusterConsensusesMap[i] = consensusesVector;
            }

        }

        //the process is ongoing until all of the clusters converge
        do {
            iterations++;
            std::cout << iterations << ". iteration" << std::endl;
            converged = true;
            //computing the distance between all centroids
            //for every query
            for (int i = 0; i < clusters.size(); i++ ) {
                clusters[i].cluster_sequences.clear();
            }
            for (int i = 0; i < sequences.size(); i++) {
                int min_diff = std::numeric_limits<int>::min();
                int belonging_cluster_index = 0;
                //for each cluster
                for (int j = 0; j < clusters.size(); j++) {
                    auto score = alignment::Align(sequences[i].c_str(), sequences[i].length(),
                                                  clusters[j].centroid.c_str(), clusters[j].centroid.length(),
                                                  match, mismatch, gap, alignment::AlignmentType::kSemiGlobal);
                    if (score > min_diff) {
                        min_diff = score;
                        belonging_cluster_index = j;
                    }
                }
                //adding the query to cluster
                clusters[belonging_cluster_index].cluster_sequences.push_back(sequences[i]);
            }
            //computing new clusters' centroids
            for (int i = 0; i < clusters.size(); i++) {

                if (clusters[i].cluster_sequences.size() != 0) {
                    auto alignment_engine = spoa::createAlignmentEngine(spoa::AlignmentType::kOV,
                                                                        match, mismatch, gap, gap);
                    auto graph = spoa::createGraph();

                    for (const auto &it: clusters[i].cluster_sequences) {
                        auto alignment = alignment_engine->align(it, graph);
                        graph->add_alignment(alignment, it);
                    }

                    std::string consensus = graph->generate_consensus();
                    if (clusters[i].centroid != consensus) converged = false;
                    clusters[i].centroid = consensus;

                    clusterConsensusesMap[i].push_back(consensus);

                } //else converged = false;
            }

            if (iterations == numberOfIterations) {
                break;
            }


        } while (converged == false);

        std::vector<std::string> centroids;
        for (auto each : clusters) {
            if (initParam == CentroidInit::kFair) {
                if (each.cluster_sequences.size() < clusterSize) {
                    continue;
                }
            }
            centroids.push_back(each.centroid + std::to_string(each.cluster_sequences.size()));
        }

        if (initParam == CentroidInit::kFair) {
            std::cout << "Number of clusters: " << centroids.size() << std::endl;
        }

        for (auto cluster : clusterConsensusesMap) {
            std::set<std::string> s(cluster.second.begin(),
                    cluster.second.end());
            std::map<std::string, int> centroidCountMap;
            for (auto centroid : cluster.second) {
                if (centroidCountMap.find(centroid) == centroidCountMap.end()) {
                    centroidCountMap[centroid] = 1;
                } else {
                    centroidCountMap[centroid]++;
                }
            }

        }

        std::cout << "Finished in " << iterations << " iterations" << std::endl;

        return centroids;
    }

    //given centroids
    std::vector<std::string> generateClusters(std::vector<std::string> sequences,
                                              int k,
                                              int numberOfIterations,
                                              std::vector<std::string> cents, int match, int mismatch, int gap) {
        bool converged = true;
        int iterations = 0;
        std::vector<kMeans::Cluster> clusters;
        std::map<int, std::vector<std::string>> clusterConsensusesMap;
        if (k > cents.size()) {
            //assigning random k centroids
            srand(time(NULL));
            for (int i = 0; i < (k - cents.size()); i++) {
                std::vector<std::string> cluster_vector;
                kMeans::Cluster newCluster;
                newCluster.centroid = sequences[rand() % sequences.size() - 1];
                newCluster.cluster_sequences = cluster_vector;
                clusters.push_back(newCluster);
            }
        }
        for (int i = 0; i < cents.size(); i++) {
            std::vector<std::string> cluster_vector;
            kMeans::Cluster newCluster;
            newCluster.centroid = cents[i];
            newCluster.cluster_sequences = cluster_vector;
            clusters.push_back(newCluster);
        }

        //the process is ongoing until all of the clusters converge
        do {
            iterations++;
            std::cout << iterations << ". iteration" << std::endl;
            converged = true;
            //computing the distance between all centroids
            //for every query
            for (int i = 0; i < clusters.size(); i++ ) {
                clusters[i].cluster_sequences.clear();
            }
            for (int i = 0; i < sequences.size(); i++) {
                int min_diff = std::numeric_limits<int>::min();
                int belonging_cluster_index = 0;
                //for each cluster
                for (int j = 0; j < clusters.size(); j++) {
                    auto score = alignment::Align(sequences[i].c_str(), sequences[i].length(),
                                                  clusters[j].centroid.c_str(), clusters[j].centroid.length(),
                                                  match, mismatch, gap, alignment::AlignmentType::kSemiGlobal);
                    if (score >= min_diff) {
                        if (score > min_diff) {
                            min_diff = score;
                            belonging_cluster_index = j;
                        } else {
                            if (clusters[belonging_cluster_index].cluster_sequences.size() <
                            clusters[j].cluster_sequences.size()) {
                                min_diff = score;
                                belonging_cluster_index = j;
                            }
                        }

                    }
                }
                //adding the query to cluster
                clusters[belonging_cluster_index].cluster_sequences.push_back(sequences[i]);
            }
            //computing new clusters' centroids
            for (int i = 0; i < clusters.size(); i++) {

                if (clusters[i].cluster_sequences.size() != 0) {
                    auto alignment_engine = spoa::createAlignmentEngine(spoa::AlignmentType::kOV,
                                                                        match, mismatch, gap, gap);
                    auto graph = spoa::createGraph();

                    for (const auto &it: clusters[i].cluster_sequences) {
                        auto alignment = alignment_engine->align(it, graph);
                        graph->add_alignment(alignment, it);
                    }

                    std::string consensus = graph->generate_consensus();
                    //std::cout << clusters[i].centroid << ' ' << consensus << std::endl;
                    if (clusters[i].centroid != consensus) converged = false;
                    clusters[i].centroid = consensus;

                    clusterConsensusesMap[i].push_back(consensus);

                } else converged = false;
            }

            if (iterations == numberOfIterations) {
                break;
            }


        } while (converged == false);

        std::vector<std::string> centroids;
        for (auto each : clusters) {
            centroids.push_back(each.centroid + std::to_string(each.cluster_sequences.size()));
        }

        for (auto cluster : clusterConsensusesMap) {
            std::cout << "For cluster number " << cluster.first + 1 << " these are the centroids: " << std::endl;
            std::set<std::string> s(cluster.second.begin(),
                                    cluster.second.end());
            std::map<std::string, int> centroidCountMap;
            for (auto centroid : cluster.second) {
                if (centroidCountMap.find(centroid) == centroidCountMap.end()) {
                    centroidCountMap[centroid] = 1;
                } else {
                    centroidCountMap[centroid]++;
                }
            }

            for(auto centroid : centroidCountMap) {
                std::cout << centroid.first << " with " << centroid.second << " appearances." << std::endl;

            }
        }

        std::cout << "Finished in " << iterations << " iterations" << std::endl;

        return centroids;
    }
}