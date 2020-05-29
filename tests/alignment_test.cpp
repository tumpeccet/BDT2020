//
// Created by nera on 11. 05. 2020..
//
#include <gtest/gtest.h>
#include "../alignment/alignment.cpp"

TEST(PairwiseAlignment, globalType) {
EXPECT_EQ(alignment::Align("GCATGCU", 7, "GATTACA", 7, 1, -1, -1, alignment::AlignmentType::kGlobal), 0);
}

TEST(PairwiseAlignment, globalType2) {
::std::string cig2;
unsigned int tb2;
EXPECT_EQ(alignment::Align("GCAT", 4, "GCT", 3, 1, -1, -1, alignment::AlignmentType::kGlobal//, cig2, tb2
), 2);
//EXPECT_EQ(tb2, 0);
//EXPECT_EQ(cig2, "2M1I1M");
}

TEST(PairwiseAlignment, localType) {
EXPECT_EQ(alignment::Align("ACCTAAGG", 8, "GGCTCAATCA", 10, 2, -1, -2, alignment::AlignmentType::kLocal), 6);
}

TEST(PairwiseAlignment, semiGlobalType) {
std::string cig;
unsigned int tb;
EXPECT_EQ (alignment::Align("TCCG", 4, "ACTCCGAT", 8, 4, -1, -2, alignment::AlignmentType::kSemiGlobal
        //, cig, tb
        ), 16);
//EXPECT_EQ (cig, "4M");
//EXPECT_EQ (tb, 2);
}

TEST(PairwiseAlignment, globalTypeCigar) {
::std::string cigarStr;
unsigned int tBegin;
int value = alignment::Align("CAGCACTTGGATTCTCGG", 18,
                              "CAGCGTGG", 8,
                              1,
                              -1,
                              -1,
                             alignment::AlignmentType::kGlobal);//,
                              //cigarStr,
                              //tBegin);
//EXPECT_EQ (cigarStr, "3M2I1M3I1M4I1M1I2M");
//EXPECT_EQ (tBegin, 0);
}

TEST(PairwiseAlignment, semiGlobalType2) {

    EXPECT_EQ (alignment::Align("CAGCACTTGGATTCTCGG", 18, "CAGCGTGG", 8, 1, -1, -2, alignment::AlignmentType::kSemiGlobal
            //, cig, tb
    ), 3);
//EXPECT_EQ (cigarStr, "3M2I1M3I1M4I1M1I2M");
//EXPECT_EQ (tBegin, 0);
}


int main(int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
