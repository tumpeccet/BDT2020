//
// Created by nera on 17. 03. 2020..
//

#ifndef BDT2020_ALIGNMENT_H
#define BDT2020_ALIGNMENT_H

#endif //BDT2020_ALIGNMENT_H

namespace alignment {
    enum AlignmentType {
        kGlobal,
        kLocal,
        kSemiGlobal
    };

    enum Direction {
        kUp,
        kLeft,
        kDiagonal,
        kNone
    };

    struct Cell {
        int value;
        Direction direction;

        bool operator<(const Cell &other) const {
            return value < other.value;
        };

        bool operator>(const Cell &other) const{
            return value > other.value;
        };

        bool operator==(const Cell &other) const{
            return value == other.value;
        };
    };

    int Align(const char *first, unsigned int first_length,
              const char *second, unsigned int second_length,
              int match,
              int mismatch,
              int gap,
              AlignmentType type);
}
