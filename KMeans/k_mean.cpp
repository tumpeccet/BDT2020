//
// Created by nera on 11. 05. 2020..
//

#include <limits>
#include <iostream>
#include <map>
#include <set>
#include <functional>
#include "k_mean.h"
#include "../spoa/include/spoa/spoa.hpp"
#include "../alignment/alignment.h"


//ograniciti broj iteracija
//ako do kraja nije pronaden centroid vratiti centroid
//  svih pronadenih centroida od kad je pocelo konvergirati ili
// nesto bolje

namespace kMeans {
    std::vector<std::string> generateClusters(std::vector<std::string> sequences,
            int k, int numberOfIterations) {
        bool converged = true;
        int iterations = 0;
        //assigning random k centroids
        std::vector<kMeans::Cluster> clusters;
        for (int i = 0; i < k; i++) {
            std::vector<std::string> cluster_vector;
            //cluster_vector.push_back(sequences[k]);
            kMeans::Cluster newCluster;
            newCluster.centroid = sequences[i];
            newCluster.cluster_sequences = cluster_vector;
            clusters.push_back(newCluster);
        }

        //ili racunanje konsenzusne sekvence i uzimanje za pocetne
        //one s najvecom udaljenosti od konsenzusne
        auto alignment_engine = spoa::createAlignmentEngine(spoa::AlignmentType::kOV,
                                                            1, -1, -1, -1);

        auto graph = spoa::createGraph();
        for (const auto& it: sequences) {
            auto alignment = alignment_engine->align(it, graph);
            graph->add_alignment(alignment, it);
        }

        std::string consensus = graph->generate_consensus();
        std::map<std::string, int> consensusAlignmentMap;
        std::map<int, std::vector<std::string>> clusterConsensusesMap;
        for (const auto& it : sequences) {
            if (consensusAlignmentMap.count(it) == 0) {
                auto score = alignment::Align(it.c_str(), it.length(),
                                              consensus.c_str(), consensus.length(),
                                              1, -1, -1, alignment::AlignmentType::kSemiGlobal);

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

        clusters.clear();
        int i = 0;
        for (std::pair<std::string, int> el : setOfCons) {
            std::vector<std::string> cluster_vector;
            kMeans::Cluster newCluster;
            newCluster.centroid = el.first;
            newCluster.cluster_sequences = cluster_vector;
            clusters.push_back(newCluster);
            std::cout << el.second;
            i++;
            if (i == k) break;
        }

        for (int i = 0; i < clusters.size(); i++) {
            std::vector<std::string> consensusesVector;
            clusterConsensusesMap[i] = consensusesVector;
        }
        //or, given consensuses
        //clusters.clear();
        //for (int i = 0; i < k; i++) {
        //    std::vector<std::string> cluster_vector;
        //    kMeans::Cluster newCluster;
        //    newCluster.centroid = consensuses[i];
        //    newCluster.cluster_sequences = cluster_vector;
        //    clusters.push_back(newCluster);
        //}

        //the process is ongoing until all of the clusters converge
        do {
            iterations++;
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
                                                  1, -1, -1, alignment::AlignmentType::kSemiGlobal);
                    if (score > min_diff) {
                        min_diff = score;
                        belonging_cluster_index = j;
                    }
                }
                //std::cout << min_diff << std::endl;
                //adding the query to cluster
                clusters[belonging_cluster_index].cluster_sequences.push_back(sequences[i]);
            }
            //computing new clusters' centroids
            for (int i = 0; i < clusters.size(); i++) {

                if (clusters[i].cluster_sequences.size() != 0) {
                    auto alignment_engine = spoa::createAlignmentEngine(spoa::AlignmentType::kOV,
                                                                        1, -1, -1, -1);
                    auto graph = spoa::createGraph();

                    for (const auto &it: clusters[i].cluster_sequences) {
                        auto alignment = alignment_engine->align(it, graph);
                        graph->add_alignment(alignment, it);
                    }

                    std::string consensus = graph->generate_consensus();
                    std::cout << clusters[i].centroid << ' ' << consensus << std::endl;
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