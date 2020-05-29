//
// Created by nera on 17. 03. 2020..
//

#include "alignment.h"
#include <algorithm>

namespace alignment {

    int Align(const char *first, unsigned int first_length,
              const char *second, unsigned int second_length,
              int match,
              int mismatch,
              int gap,
              AlignmentType type) {

        int rowCount = first_length + 1;
        int colCount = second_length + 1;
        int alignmentRow = 0, alignmentCol = 0, optimalAlignment;
        int rowInd = 0, columnInd = 0, maxScore = 0;
        // creating a matrix consisted of Cell structure
        // declared in header file
        Cell **matrix = new Cell * [rowCount];
        for (int i = 0; i < rowCount; i++) {
            matrix[i] = new Cell[colCount];
        }
        //matrix initialization; 1st value is always (0,0)
        matrix[0][0] = {0, Direction::kNone};

        //setting first row and column values if the type is global
        if (type == AlignmentType::kGlobal) {
            for (rowInd = 1; rowInd < rowCount; rowInd++) {
                matrix[rowInd][0] = {gap * rowInd, Direction::kUp};
            }

            for (columnInd = 1; columnInd < colCount; columnInd++) {
                matrix[0][columnInd] = {gap * columnInd, Direction::kLeft};
            }
        } else { // for both local and semi-global type first row and column values are 0
            for (rowInd = 1; rowInd < rowCount; rowInd++) {
                matrix[rowInd][0] = {0, Direction::kNone};
            }

            for (columnInd = 1; columnInd < colCount; columnInd++) {
                matrix[0][columnInd] = {0, Direction::kNone};
            }
        }
        //initializing the rest of the matrix
        for (rowInd = 1; rowInd < rowCount; rowInd++) {
            for (columnInd = 1; columnInd < colCount; columnInd++) {
                //determining if the target is a match or a mismatch
                int w = first[rowInd - 1] == second[columnInd - 1] ? match : mismatch;
                int replacement = matrix[rowInd - 1][columnInd - 1].value + w;

                int insertion = matrix[rowInd][columnInd - 1].value + gap;
                int deletion = matrix[rowInd - 1][columnInd].value + gap;

                //computing the max value of parent + added cell values
                Cell maxValue = std::max({
                    Cell{replacement, Direction::kDiagonal},
                    Cell{insertion, Direction::kLeft},
                    Cell{deletion, Direction::kUp}
                });
                matrix[rowInd][columnInd] = maxValue;
                if (type == AlignmentType::kLocal) {
                    if (maxValue.value < 0) {
                        matrix[rowInd][columnInd].value = 0;
                    }
                    if (maxValue.value > maxScore) {
                        //maxScore = matrix[rowInd][columnInd].value;
                        maxScore = maxValue.value;
                        alignmentRow = rowInd;
                        alignmentCol = columnInd;
                    }
                }
            }
        }
        if (type == AlignmentType::kGlobal) {
            //global alignment score is the last cell value
            maxScore = matrix[rowCount - 1][colCount - 1].value;
            alignmentRow = rowCount - 1;
            alignmentCol = colCount - 1;
        } else if (type == AlignmentType::kSemiGlobal) {
            maxScore = matrix[1][colCount - 1].value;
            alignmentRow = 0, alignmentCol = colCount - 1;
            for (int i = 2; i < rowCount; i++) {
                if (matrix[i][colCount - 1].value > maxScore) {
                    maxScore = matrix[i][colCount - 1].value;
                    alignmentRow = i;
                }
            }
            for (int j = 1; j < colCount; j++) {
                if (matrix[rowCount - 1][j].value > maxScore) {
                    maxScore = matrix[rowCount - 1][j].value;
                    alignmentRow = rowCount - 1;
                    alignmentCol = j;
                }
            }
        }


        for (int i = 0; i < rowCount; i++) {
            delete[] matrix[i];
        }

        delete[] matrix;

        return maxScore;
    }
}