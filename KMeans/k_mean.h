//
// Created by nera on 11. 05. 2020..
//

#include <string>
#include <vector>

#ifndef BDT2020_K_MEAN_H
#define BDT2020_K_MEAN_H

#endif //BDT2020_K_MEAN_H

namespace kMeans{
    enum CentroidInit {
        kFair,
        kWorst,
        kRandom
    };

    struct Cluster {
        std::string centroid;
        std::vector<std::string> cluster_sequences;

        bool operator==(const Cluster& a) const {
            return (centroid == a.centroid && cluster_sequences == a.cluster_sequences);
        }

        bool operator<(const Cluster& a) const {
            return (cluster_sequences < a.cluster_sequences);
        }

        bool operator>(const Cluster& a) const {
            return (cluster_sequences > a.cluster_sequences);
        }



    };

    std::vector<std::string> generateClusters(std::vector<std::string> sequences,
            int k, int numberOfIterations, int match, int mismatch, int gap,
            CentroidInit initParam, int newClusterScore, int clusterSize);

    std::vector<std::string> generateClusters(std::vector<std::string> sequences,
                                              int k, int numberOfIterations,
                                              std::vector<std::string> cents,
                                              int match, int mismatch, int gap);

}