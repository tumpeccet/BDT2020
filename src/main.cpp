//
// Created by nera on 17. 03. 2020..
//

#include "main.h"
#include <iostream>
#include <algorithm>
#include <map>
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

auto evaluating() {
    //evaluacija

    //for(const auto& elem : fasta_objects){
    //    std::cout << elem->sequence << std::endl;
    //}
//    map<string, map<int, int>> expectedDifsMap;
//    for (int j = 0; j < fasta_objects.size(); j++) {
//        map<int, int> newMap;
//        expectedDifsMap[fasta_objects[j]->sequence] = newMap;
//    }
//    for (int i = 0; i < fastq_objects.size(); i++) {
//        for (int j = 0; j < fasta_objects.size(); j++) {
//            //poravnaj i spremi u mapu mape
//            auto score = alignment::Align(fastq_objects[i]->sequence.c_str(), fastq_objects[i]->sequence.length(),
//                                          fasta_objects[j]->sequence.c_str(), fasta_objects[j]->sequence.length(), 1, 0, -1,
//                                          alignment::AlignmentType::kSemiGlobal);
//            auto min = fastq_objects[i]->sequence.length();
//            if (fasta_objects[j]->sequence.length() < fastq_objects[i]->sequence.length()) min = fasta_objects[j]->sequence.length();
//            if (expectedDifsMap[fasta_objects[j]->sequence].find(min - score) == expectedDifsMap[fasta_objects[j]->sequence].end())
//                expectedDifsMap[fasta_objects[j]->sequence][min - score] = 1;
//            else
//                expectedDifsMap[fasta_objects[j]->sequence][min - score] += 1;
//        }
//    }
//    map<string, map<int, int>>::iterator it2 = expectedDifsMap.begin();
//    evaluationFile.open("../evaluation.txt", ios::out | ios::app);
//    evaluationFile << "J29B_expected" << endl;
//    while (it2 != expectedDifsMap.end()) {
//        cout << it2->first << endl;
//        evaluationFile << "Alel (" << it2->first.size()
//                       << ", " << ")" << '\n';
//        evaluationFile << it2->first << '\n' << '\n';
//        map<int, int>::iterator it_inner = it2->second.begin();
//        while (it_inner != it2->second.end()) {
//            cout << it_inner->first << ' ' << it_inner->second << endl;
//            evaluationFile << it_inner->first << ' ' << it_inner->second << '\n';
//            it_inner++;
//        }
//        it2++;
//    }
//    evaluationFile << "--------------------------------------------" << '\n';
//
//    map<string, map<int, int>> expectedDifsMap_30;
//    for (int j = 0; j < fasta_objects_30.size(); j++) {
//        map<int, int> newMap;
//        expectedDifsMap_30[fasta_objects_30[j]->sequence] = newMap;
//    }
//    for (int i = 0; i < fastq_objects_30.size(); i++) {
//        for (int j = 0; j < fasta_objects_30.size(); j++) {
//            //poravnaj i spremi u mapu mape
//            auto score = alignment::Align(fastq_objects_30[i]->sequence.c_str(), fastq_objects_30[i]->sequence.length(),
//                                          fasta_objects_30[j]->sequence.c_str(), fasta_objects_30[j]->sequence.length(), 1, 0, -1,
//                                          alignment::AlignmentType::kSemiGlobal);
//            auto min = fastq_objects_30[i]->sequence.length();
//            if (fasta_objects_30[j]->sequence.length() < fastq_objects_30[i]->sequence.length()) min = fasta_objects_30[j]->sequence.length();
//            if (expectedDifsMap_30[fasta_objects_30[j]->sequence].find(min - score) == expectedDifsMap_30[fasta_objects_30[j]->sequence].end())
//                expectedDifsMap_30[fasta_objects_30[j]->sequence][min - score] = 1;
//            else
//                expectedDifsMap_30[fasta_objects_30[j]->sequence][min - score] += 1;
//        }
//    }
//    map<string, map<int, int>>::iterator it3 = expectedDifsMap_30.begin();
//    evaluationFile << "J30B_expected" << endl;
//    while (it3 != expectedDifsMap_30.end()) {
//        cout << it3->first << endl;
//        evaluationFile << "Alel (" << it3->first.size()
//                       << ", " << ")" << '\n';
//        evaluationFile << it3->first << '\n' << '\n';
//        map<int, int>::iterator it_inner_2 = it3->second.begin();
//        while (it_inner_2 != it3->second.end()) {
//            cout << it_inner_2->first << ' ' << it_inner_2->second << endl;
//            evaluationFile << it_inner_2->first << ' ' << it_inner_2->second << '\n';
//            it_inner_2++;
//        }
//        it3++;
//    }
//    evaluationFile << "--------------------------------------------" << '\n';
//    evaluationFile.close();
}

int main(int argc, char* argv[]){
    ofstream firstAlgorithmFile, secondAlgorithmFile, kMeansAlgorithmFile, evaluationFile;
    std::vector<std::unique_ptr<FASTQformat>> fastq_objects, fastq_objects_30;
    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FASTQformat>("../Bioinformatika-jeleni/fastq/J29_B_CE_IonXpress_005.fastq");
    fastq_parser->parse(fastq_objects, -1);
//    for (int i = 0; i < fastq_objects.size(); i++){
//        std::cout << fastq_objects[i]->sequence;
//    }
    auto fastq_parser_30 = bioparser::createParser<bioparser::FastqParser, FASTQformat>("../Bioinformatika-jeleni/fastq/J30_B_CE_IonXpress_006.fastq");
    fastq_parser_30->parse(fastq_objects_30, -1);

    std::vector<std::string> sequences = {
            "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
            "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
            "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
            "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
            "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
            "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
    };

    vector<string> kMeansSequences = {
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATACTACGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTACGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGTCCGCCAAGTACTGGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACAGGTACTGCAGACACAACTACGGGGTTCTTGACAGTTTCGCTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGGAGTATCATAAGAGCGAGTGTCATTTCTTCAACGGGACCGAGCGGGTGCGGTTCCTGGACAGATACTTCTATAATGGAGAAGAGTACGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGGCCGAGCTGGGGCGGCCGGCCGCCGAGCACTGGAACAGCCAGAAGGAGATTCTGGAGCAGAAGCGGGCCGCGGTGGACAGGTACTGCAGACACAACTACGGGGTCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA"

    };

    vector<string> kMeansSemiGlobalSequences = {
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTTTTTAACGGAGAAGAAAACGTGCCCCCAAAAAAGGGGGGCGTTTTTTTTAAAAAAAACCCCCGGGGGGTTTTTTACAGAACTGGGGGCCCCCGGGGGGTTTTTTTTCCCCCCCCTTTTTTTTGGGGGGGGAAAAAAAAGGGGGGGCCCCCCCCAAAAAAAATTTTTTTTCCCCCCCCGGGGGAAAAAAAATTTTTTCCCCGGGGAAAAAATTTTGGG",

            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACCCCCGGGGGGGTTTTTCCCCCTTTTTGGGGGGAAAAAATTTTTTTCCCCCCCGGGGGGAAAAAATTTTTGGGGGTTTTTAAAAACCCCCAAAAAGGGGGGCCCCCTTTTTAACCCCGGGGGGTTTTTAAAACCCCCGGGGGGTTTTTCCCCTTTTGGGGGGAAAAAGCCCAAAATTTCCGGGAAAATTCGAATACCATTTCCTTAATTAAGGGAAAAAT",

            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATACTACGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTACGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGTCCGCCAAGTACTGGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACAGGTACTGCAGACACAACTACGGGGTTCTTGACAGTTTCGCTGGTGCAGCGGTCGAGGTGACGCGAAGGTACAA",

            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAAAGGAGTACGT"
    };

    vector<string> secondSequences = {
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATACTACGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTACGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGTCCGCCAAGTACTGGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACAGGTACTGCAGACACAACTACGGGGTTCTTGACAGTTTCGCTGGTGCAGCGGTCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTTCACTGTGCAGCGGCGAGGTGACGCGAA"
    };

    vector<string> firstSequences = {
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA"
    };

    vector<string> firstSequencesLonger = {
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA",
            "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA"
    };


    std::vector<std::string> msa;

    std::vector<std::unique_ptr<FASTAformat>> fasta_objects;
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAformat>("../Bioinformatika-jeleni/J29B_expected.fasta");
    fasta_parser->parse(fasta_objects, -1);

    std::vector<std::unique_ptr<FASTAformat>> fasta_objects_30;
    auto fasta_parser_30 = bioparser::createParser<bioparser::FastaParser, FASTAformat>("../Bioinformatika-jeleni/J30B_expected.fasta");
    fasta_parser_30->parse(fasta_objects_30, -1);


    cout << "k - means semi-global shorter: " << endl;
    for (int i = 0; i < fasta_objects.size(); i++) {
        auto maxScore = 0;
        string theSequence;
        for (auto seq : kMeansSemiGlobalSequences) {
            auto score = alignment::Align(fasta_objects[i]->sequence.c_str(), fasta_objects[i]->sequence.length(),
                                          seq.c_str(), seq.length(), 1, 0, -1,
                                          alignment::AlignmentType::kSemiGlobal);
            if (score > maxScore) {
                maxScore = score;
                theSequence = seq;
            }

        }
        cout << fasta_objects[i]->sequence << ' ' << fasta_objects[i]->sequence.length() << endl;
        cout << theSequence << ' ' << theSequence.length() << endl;
        if (theSequence.length() < fasta_objects[i]->sequence.length()) cout << theSequence.length() - maxScore << endl;
        else cout << fasta_objects[i]->sequence.length() - maxScore << endl;
    }
    cout << "k - means: " << endl;
    for (int i = 0; i < fasta_objects.size(); i++) {
        auto maxScore = 0;
        for (auto seq : kMeansSequences) {
            auto score = alignment::Align(fasta_objects[i]->sequence.c_str(), fasta_objects[i]->sequence.length(),
                                          seq.c_str(), seq.length(), 1, 0, 0,
                                          alignment::AlignmentType::kSemiGlobal);
            if (score > maxScore) maxScore = score;
        }
        cout << fasta_objects[i]->sequence << endl;
        cout << fasta_objects[i]->sequence.length() - maxScore << endl;
    }
    cout << "second algorithm: " << endl;
    for (int i = 0; i < fasta_objects.size(); i++) {
        auto maxScore = 0;
        for (auto seq : secondSequences) {
            auto score = alignment::Align(fasta_objects[i]->sequence.c_str(), fasta_objects[i]->sequence.length(),
                                          seq.c_str(), seq.length(), 1, 0, 0,
                                          alignment::AlignmentType::kSemiGlobal);
            if (score > maxScore) maxScore = score;
        }
        cout << fasta_objects[i]->sequence << endl;
        cout << fasta_objects[i]->sequence.length() - maxScore << endl;
    }
    cout << "first algorithm: " << endl;
    for (int i = 0; i < fasta_objects.size(); i++) {
        auto maxScore = 0;
        for (auto seq : firstSequences) {
            auto score = alignment::Align(fasta_objects[i]->sequence.c_str(), fasta_objects[i]->sequence.length(),
                                          seq.c_str(), seq.length(), 1, 0, 0,
                                          alignment::AlignmentType::kSemiGlobal);
            if (score > maxScore) maxScore = score;
        }
        cout << fasta_objects[i]->sequence << endl;
        cout << fasta_objects[i]->sequence.length() - maxScore << endl;
    }
    cout << "first algorithm longer: " << endl;
    for (int i = 0; i < fasta_objects.size(); i++) {
        auto maxScore = 0;
        for (auto seq : firstSequencesLonger) {
            auto score = alignment::Align(fasta_objects[i]->sequence.c_str(), fasta_objects[i]->sequence.length(),
                                          seq.c_str(), seq.length(), 1, 0, 0,
                                          alignment::AlignmentType::kSemiGlobal);
            if (score > maxScore) maxScore = score;
        }
        cout << fasta_objects[i]->sequence << endl;
        cout << fasta_objects[i]->sequence.length() - maxScore << endl;
    }
    for (int i = 0; i < fasta_objects.size() - 1; i++) {
        for (int j = 1; j < fasta_objects.size(); j++) {
            if (i == j) continue;
            auto score = alignment::Align(fasta_objects[i]->sequence.c_str(), fasta_objects[i]->sequence.length(),
                fasta_objects[j]->sequence.c_str(), fasta_objects[j]->sequence.length(), 1, 0, -1,
                alignment::AlignmentType::kGlobal);
            int smaller = fasta_objects[i]->sequence.length();
            if (fasta_objects[j]->sequence.length() < smaller)
                smaller = fasta_objects[j]->sequence.length();
            std::cout << smaller << ',' << score << std::endl;
            std::cout << smaller - score << std::endl;
            if (fasta_objects[i]->sequence != fasta_objects[j]->sequence)
                std::cout << "Not equal (";
            else
                std::cout << "Equal (";
            std::cout << i+1 << '.' << ',' << j+1 << '.' << ")" << std::endl;
        }
    }

    std::cout << "-------------------------------------" << std::endl;
    //filtering by length
    std::map<int, int> lengths;
    std::vector<std::string> filtered_objects;
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
        if (it->sequence.length() == max_length->first) filtered_objects.push_back(it->sequence);
    }
    std::cout << "Length of " << max_length->first << " nucleotide bases is most common (" << max_length->second << " appearances)" << std::endl;
    auto alignment_engine_1 = spoa::createAlignmentEngine(spoa::AlignmentType::kNW,
                                                        0, -1, -100, -100);

    auto graph_1 = spoa::createGraph();

    for (const auto& it: filtered_objects) {
        auto alignment_1 = alignment_engine_1->align(it, graph_1);
        graph_1->add_alignment(alignment_1, it);
    }

    std::string consensus_1 = graph_1->generate_consensus();

    std::cout <<"Consensus (" << consensus_1.size() << ")" << std::endl;
    std::cout << consensus_1.c_str() << std::endl;

    vector<string> filtered_objects_longer;
    for (const auto& it : fastq_objects) {
        if (it->sequence.length() <= max_length->first + 5
            && it->sequence.length() >= max_length->first - 5) filtered_objects_longer.push_back(it->sequence);
    }

    vector<string> filtered_objects_shorter;
    for (const auto& it : fastq_objects) {
        if (it->sequence.length() <= max_length->first + 5
            && it->sequence.length() >= 200) filtered_objects_shorter.push_back(it->sequence);
    }
    // odrediti maksimalni broj iteracija
    auto kMeansClusters = kMeans::generateClusters(filtered_objects_shorter, 4, 100);
    kMeansAlgorithmFile.open("../kMeans.txt", ios::out | ios::app);
    std::cout << "K-means algorithm: " << std::endl;
    for (auto it : kMeansClusters){
        std::cout <<"Consensus (" << it.size() << ")" << std::endl;
        std::cout << it << std::endl;
        std::cout << std::endl;
        kMeansAlgorithmFile << "Consensus (" << it.size()
        << ", " << ")" << '\n';
        kMeansAlgorithmFile << it << '\n' << '\n';
    }
    kMeansAlgorithmFile << "results: " << '\n';
    for (int i = 0; i < fasta_objects.size(); i++) {
        auto maxScore = 0;
        string theSequence;
        for (auto seq : kMeansClusters) {
            auto score = alignment::Align(fasta_objects[i]->sequence.c_str(), fasta_objects[i]->sequence.length(),
                                          seq.c_str(), seq.length(), 1, 0, -1,
                                          alignment::AlignmentType::kSemiGlobal);
            if (score > maxScore) {
                maxScore = score;
                theSequence = seq;
            }

        }
        kMeansAlgorithmFile << fasta_objects[i]->sequence << ' ' << fasta_objects[i]->sequence.length() << endl;
        kMeansAlgorithmFile << theSequence << ' ' << theSequence.length() << endl;
        if (theSequence.length() < fasta_objects[i]->sequence.length()) kMeansAlgorithmFile << theSequence.length() - maxScore << endl;
        else kMeansAlgorithmFile << fasta_objects[i]->sequence.length() - maxScore << endl;
    }
    kMeansAlgorithmFile << "--------------------------------------------" << '\n';
    kMeansAlgorithmFile.close();

    //DRUGI ALGORITAM

    std::vector<std::vector<std::string>> clusters;
    std::map<std::string, std::vector<std::string>> cluster_consensus_map;
    int count = 1;
    int OK_score = -10;
    std::vector<std::string> vec;
    std::string first_allel = filtered_objects_longer.back();
    vec.push_back(first_allel);
    clusters.push_back(vec);
    cluster_consensus_map[first_allel] = vec;
    bool matched = false, this_cluster = true;
    std::string newKey;
    std::cout << "------------------------" << std::endl;
    std::vector<std::string> consensuses;
    int match = 0, mismatch = -1, gap = -1;
    for (int i = 0; i < filtered_objects_longer.size(); i++) {
        //std::cout << "tu smo " << std::endl;
        auto new_obj = filtered_objects_longer[i];
        std::map<std::string, std::vector<std::string>>::iterator it = cluster_consensus_map.begin();
        do {
            auto cluster_consensus = it->first;
            auto cluster_vector = it->second;
            //std::cout << it->first << std::endl;
            //for (auto it : cluster_vector) {
            auto score = alignment::Align(cluster_consensus.c_str(), cluster_consensus.length(), new_obj.c_str(), new_obj.length(), match, mismatch, gap,
                        alignment::AlignmentType::kGlobal);
            //std::cout << score << std::endl;
            if (score < OK_score) {
                this_cluster = false;
//                if (it == cluster_consensus_map.end()) break;
//                it++;
//                continue;
            }
            //}
            if (this_cluster == false) {
                this_cluster = true;
                if (it == cluster_consensus_map.end()) break;
                it++;
            } else {
                it->second.push_back(new_obj);
                //tu napraviti konsezusnu sekvencu
                matched = true;
                auto alignment_engine_3 = spoa::createAlignmentEngine(spoa::AlignmentType::kNW,
                                                                      match, mismatch, gap, gap);

                auto graph_3 = spoa::createGraph();

                for (const auto& alel : it->second) {
                    auto alignment_3 = alignment_engine_3->align(alel, graph_3);
                    graph_3->add_alignment(alignment_3, alel);
                }
                newKey = graph_3->generate_consensus();
                break;
            }
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
            std::cout << newKey << " zamjenjuje " << it->first << std::endl;
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

    int the_size = 5, i = 0;
    secondAlgorithmFile.open("../second_algorithm.txt", ios::out | ios::app);
    secondAlgorithmFile << "Minimal size of clusters: " << the_size << '\n';
    secondAlgorithmFile << "Parameters: " << '\n';
    secondAlgorithmFile << "match = " << match << ", mismatch = "
    << mismatch << ", gap = " << gap <<
    ", score for making a new cluster = " << OK_score << '\n';
    //UPISI TIP PORAVNANJA
    secondAlgorithmFile << "Global alignment" << '\n';
    secondAlgorithmFile << "Number of clusters: " << cluster_consensus_map.size() << '\n';
    std::map<std::string, std::vector<std::string>>::iterator it = cluster_consensus_map.begin();
    auto len = it->second.size();
    do {
        std::cout << it->second.size() << std::endl;
        if (it->second.size() < the_size) {
            if (it != cluster_consensus_map.end()){
                it++;
                continue;
            } else break;
        }
        std::cout << '!' << std::endl;
        std::cout <<"Consensus (" << it->first.size() <<
        ", " << it->second.size() << ")" << std::endl;
        std::cout << it->first << std::endl;
        consensuses.push_back(it->first);
        secondAlgorithmFile << "Consensus (" << it->first.size() <<
        ")" << '\n';
        secondAlgorithmFile << it->first << '\n' << '\n';
        it++;
        len = it->second.size();
    } while(it != cluster_consensus_map.end());
    secondAlgorithmFile << "----------------------------------" << '\n';
    secondAlgorithmFile.close();


// PRVI ALGORITAMMMMMM
    int OK_score2 = -10;
    int match2 = 0, mismatch2 = -1, gap2 = -100;
    for (int i = 0; i < filtered_objects_longer.size(); i++) {
        auto new_obj = filtered_objects_longer[i];
        for (int j = 0; j < clusters.size(); j++) {
            auto cluster_vector = clusters[j];
            for (auto it : cluster_vector) {
                auto score = alignment::Align(it.c_str(), it.length(), new_obj.c_str(), new_obj.length(), match2, mismatch2, gap2,
                        alignment::AlignmentType::kGlobal);
                //std::cout << score << std::endl;
                if (score < OK_score2) {
                    this_cluster = false;
                    break;
                }
            }
            if (this_cluster == false) {
                this_cluster = true;
                continue;
            } else {
                clusters[j].push_back(new_obj);
                //tu napraviti konsezusnu sekvencu
                matched = true;
                break;
            }
        }
        if (matched == false) {
            std::vector<std::string> vec;
            vec.push_back(new_obj);
            clusters.push_back(vec);
        }
        matched = false;
    }

    //for (auto it: clusters) {
    //    std::cout << it[0] << std::endl;
    //}
    std::cout << clusters.size() << std::endl;
    std::sort(clusters.begin(), clusters.end(),
            [](const std::vector<std::string>& v1, const std::vector<std::string>& v2){
        return v1.size() > v2.size();
    });
    int the_size2 = 10, i2 = 0;
    auto len2 = clusters[i2].size();
    firstAlgorithmFile.open("../first_algorithm.txt", ios::out | ios::app);
    firstAlgorithmFile << "Minimal size of clusters: " << the_size2 << '\n';
    firstAlgorithmFile << "Parameters: " << '\n';
    firstAlgorithmFile << "match = " << match2 << ", mismatch = "
                       << mismatch2 << ", gap = " << gap2 <<
                       ", score for making a new cluster = " << OK_score2 << '\n';
    firstAlgorithmFile << "Number of clusters: " << clusters.size() << '\n';
    while(len2 >= the_size2) {
        std::cout << '!' << std::endl;
        auto alignment_engine_2 = spoa::createAlignmentEngine(spoa::AlignmentType::kNW,
                                                              0, -1, -100, -100);

        auto graph_2 = spoa::createGraph();

        for (const auto& it: clusters[i]) {
            auto alignment_2 = alignment_engine_2->align(it, graph_2);
            graph_2->add_alignment(alignment_2, it);
        }

        std::string consensus_2 = graph_2->generate_consensus();

        std::cout <<"Consensus (" << consensus_2.size() << ")" << std::endl;
        std::cout << consensus_2.c_str() << std::endl;
        firstAlgorithmFile << "Consensus (" << consensus_2.size() << ")" << '\n';
        firstAlgorithmFile << consensus_2.c_str() << '\n' << '\n';

        i2++;
        len2 = clusters[i2].size();
    }
    firstAlgorithmFile << "-----------------------------------" << '\n';
    firstAlgorithmFile.close();



    return 0;
}