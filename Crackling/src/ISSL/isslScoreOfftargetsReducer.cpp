/*

 Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


To compile:

g++ -o reducer isslScoreOfftargetsReducer.cpp 
*/

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdint.h>

#include <bitset>
#include <fstream>
#include <iostream>
#include <climits>
#include <stdio.h>
#include <cstring>

#include <bits/stdc++.h>
using namespace std;


vector<char> signatureIndex(4);
size_t seqLength;

struct MapperResult {
    uint64_t querySignature;
    uint32_t targetId;
    double mitScore;
    double cfdScore;
};

struct Key {
    uint64_t query;
    uint32_t target;

    bool operator<(const Key& other) const {
        return (query < other.query) ||
               (query == other.query && target < other.target);
    }
};




/**
 * Binary encode genetic string `ptr`
 *
 * For example, 
 *   00 11 01 10 becomes (as 64-bit unsigned int)
 *    A  T  C  G  (without spaces)
 *
 * @param[in] signature the binary encoded genetic string
 */
string signatureToSequence(uint64_t signature)
{
    string sequence = string(seqLength, ' ');
    for (size_t j = 0; j < seqLength; j++) {
        sequence[j] = signatureIndex[(signature >> (j * 2)) & 0x3];
    }
    return sequence;
}


// The reducer for the AWS lambda may have to change, here we are loading
// everything into memory

int main(int argc, char** argv)
{
    if (argc < 3) {
        fprintf(stderr,
            "Usage: %s [output-file] [sequence-length] [mapper1] [mapper2]...\n",
            argv[0]);
        exit(1);
    }

    seqLength = atoi(argv[2]);

    /** Char to binary encoding */
    signatureIndex[0] = 'A';
    signatureIndex[1] = 'C';
    signatureIndex[2] = 'G';
    signatureIndex[3] = 'T';


    const char* outFile = argv[1];

    vector<MapperResult> all;
    all.reserve(10000000); // optional tuning - im assuming 10 million results, not scalable
    // this can be a start but we are essentially loading in ALL the duplicate pairs as well
    

    // Load all mapper outputs
    for (int i = 3; i < argc; i++) {
        ifstream in(argv[i], ios::binary);

        if (!in) {
            cerr << "Failed to open " << argv[i] << endl;
            continue;
        }

        MapperResult r;
        while (in.read(reinterpret_cast<char*>(&r), sizeof(r))) {
            all.push_back(r);
        }
    }

    if (all.empty()) {
        cerr << "No data found\n";
        return 0;
    }

    // Sort by (query, target) - place all identical query sigantures 
    // followed by target id next to eachother
    sort(all.begin(), all.end(),
        [](const MapperResult& a, const MapperResult& b) {
            if (a.querySignature != b.querySignature)
                return a.querySignature < b.querySignature;
            return a.targetId < b.targetId;
        });

    // Streaming reduction
    ofstream out(outFile);

    // uint64_t curQuery = all[0].querySignature;
    // uint32_t curTarget = all[0].targetId;

    // double mitSum = 0.0;
    // double cfdSum = 0.0;

    size_t i = 0;

    while (i < all.size()) {

        uint64_t query = all[i].querySignature;

        double mitSum = 0.0;
        double cfdSum = 0.0;

        std::unordered_set<uint32_t> seenTargets;
        while (i < all.size() &&
            all[i].querySignature == query) {

            auto &r = all[i];

            if (seenTargets.insert(r.targetId).second) {
                mitSum += r.mitScore;
                cfdSum += r.cfdScore;
            }

            i++;
        }

        double finalMit = 10000.0 / (100.0 + mitSum);
        double finalCfd = 10000.0 / (100.0 + cfdSum);

        out << signatureToSequence(query)
            << "\t" << finalMit
            << "\t" << finalCfd << "\n";
    }

    cout << "Reducer wrote final output to " << outFile << "\n";

    return 0;
}

/*

// Faster and better CRISPR guide RNA design with the Crackling method.
// Jacob Bradford, Timothy Chappell, Dimitri Perrin
// bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


// To compile:

// g++ -o reducer isslScoreOfftargetsReducer.cpp

*/

// #include <cstdio>
// #include <cstdlib>
// #include <cstdint>
// #include <vector>
// #include <string>

// #include <sys/types.h>
// #include <sys/stat.h>
// #include <unistd.h>
// #include <stdint.h>

// #include <bitset>
// #include <fstream>
// #include <iostream>
// #include <climits>
// #include <stdio.h>
// #include <cstring>

// #include <bits/stdc++.h>

// using namespace std;

// vector<char> signatureIndex(4);
// size_t seqLength;

// struct MapperResult {
//     uint64_t querySignature;
//     uint32_t targetId;
//     double mitScore;
//     double cfdScore;
// };

// struct Key {
//     uint64_t query;
//     uint32_t target;

//     bool operator<(const Key& other) const {
//         return (query < other.query) ||
//                (query == other.query && target < other.target);
//     }
// };

// /**
//  * Binary encode genetic string ptr
//  *
//  * For example,
//  * 00 11 01 10 becomes (as 64-bit unsigned int)
//  * A T C G (without spaces)
//  *
//  * @param[in] signature the binary encoded genetic string
//  */
// string signatureToSequence(uint64_t signature)
// {
//     string sequence = string(seqLength, ' ');
//     for (size_t j = 0; j < seqLength; j++) {
//         sequence[j] = signatureIndex[(signature >> (j * 2)) & 0x3];
//     }
//     return sequence;
// }

// // The reducer for the AWS lambda may have to change, here we are loading
// // everything into memory

// int main(int argc, char** argv)
// {
//     if (argc < 3) {
//         fprintf(stderr,
//             "Usage: %s [output-file] [sequence-length] [mapper1] [mapper2]...\n",
//             argv[0]);
//         exit(1);
//     }

//     seqLength = atoi(argv[2]);

//     /** Char to binary encoding */
//     signatureIndex[0] = 'A';
//     signatureIndex[1] = 'C';
//     signatureIndex[2] = 'G';
//     signatureIndex[3] = 'T';

//     const char* outFile = argv[1];

//     vector<MapperResult> all;
//     all.reserve(10000000); // optional tuning - im assuming 10 million results, not scalable
//     // this can be a start but we are essentially loading in ALL the duplicate pairs as well

//     // Load all mapper outputs
//     for (int i = 3; i < argc; i++) {
//         ifstream in(argv[i], ios::binary);

//         if (!in) {
//             cerr << "Failed to open " << argv[i] << endl;
//             continue;
//         }

//         MapperResult r;
//         while (in.read(reinterpret_cast<char*>(&r), sizeof(r))) {
//             all.push_back(r);
//         }
//     }

//     if (all.empty()) {
//         cerr << "No data found\n";
//         return 0;
//     }

//     // Sort by (query, target) - place all identical query sigantures
//     // followed by target id next to eachother
//     sort(all.begin(), all.end(),
//         [](const MapperResult& a, const MapperResult& b) {
//             if (a.querySignature != b.querySignature)
//                 return a.querySignature < b.querySignature;
//             return a.targetId < b.targetId;
//         });

//     // Streaming reduction
//     ofstream out(outFile);

//     uint64_t curQuery = all[0].querySignature;
//     uint32_t curTarget = all[0].targetId;

//     double mitSum = 0.0;
//     double cfdSum = 0.0;

//     //
//     for (size_t i = 0; i < all.size(); i++) {
//         auto &r = all[i];

//         // Skip duplicates
//         if (i > 0) {
//             if (r.querySignature == all[i - 1].querySignature &&
//                 r.targetId == all[i - 1].targetId)
//                 continue;
//         }

//         // new group → emit previous as list is sorted
//         if (r.querySignature != curQuery) {

//             double finalMit = 10000.0 / (100.0 + mitSum);
//             double finalCfd = 10000.0 / (100.0 + cfdSum);

//             auto querySequence = signatureToSequence(curQuery);

//             out << querySequence << "\t"
//                 << finalMit << "\t"
//                 << finalCfd << "\n";

//             // reset group
//             curQuery = r.querySignature;
//             curTarget = r.targetId;
//             mitSum = 0.0;
//             cfdSum = 0.0;
//         }

//         mitSum += r.mitScore;
//         cfdSum += r.cfdScore;
//     }

//     // flush last group
//     double finalMit = 10000.0 / (100.0 + mitSum);
//     double finalCfd = 10000.0 / (100.0 + cfdSum);

//     auto querySequence = signatureToSequence(curQuery);

//     out << querySequence << "\t"
//         << finalMit << "\t"
//         << finalCfd << "\n";

//     cout << "Reducer wrote final output to " << outFile << "\n";

//     return 0;
// }