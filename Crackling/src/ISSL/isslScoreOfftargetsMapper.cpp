
/*
DEVELOPMENT NOTES
-----------------
This distributed Mapper is derived from the original Crackling off-target
scorer at commit e247e9b:

https://github.com/bmds-lab/Crackling/blob/e247e9b48ca9ad8de04258d51295126250ef1c26/src/ISSL/isslScoreOfftargets.cpp

It retains the original MIT/CFD scoring logic, but consumes only the hydrated
candidates for the bucket selected by this Mapper's semantic slice. Each
candidate embeds its signature, original global ID and occurrence count, so the
Mapper neither maps the global off-target catalogue nor allocates global-ID
"seen" bitmaps. It emits per-off-target score contributions for the downstream
Reducer in the Coordinator/Mapper/Reducer pipeline.

Each OpenMP thread writes to a separate temporary binary file. The files are
concatenated into one shard result after scoring. Under the AWS Lambda wrapper,
these temporary files are created in the invocation's working directory under
/tmp.

Mapper output contract
----------------------
The shard output is a headerless sequence of binary MapperResult records. Each
record contains, in order:

    uint64_t querySignature
    uint32_t targetId
    4 bytes of structure padding
    double mitScore
    double cfdScore

On the current Linux x86_64 target, each record is 32 bytes and uses the host's
little-endian byte order. The AWS Lambda wrapper reads the same layout with the
Python struct format `<QI4xdd`. A targetId of UINT32_MAX is the sentinel for a
query that emitted no off-target contributions. This contract must remain in
sync with MapperResult and every downstream reader.


Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


To compile:

g++ -o mapper isslScoreOfftargetsMapper.cpp -O3 -std=c++11 -fopenmp -mpopcnt -Iinclude
*/

#include "cfdPenalties.h"

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/time.h>
#include <chrono>
#include <bitset>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <omp.h>
#include <phmap.h>
#include <map>

#include <fstream>

using namespace std;

#define BUF_SIZE 1000

size_t seqLength, seqCount, sliceWidth, sliceCount, offtargetsCount, scoresCount;

/** Binary record emitted by the Mapper and consumed by the downstream Reducer. */
struct MapperResult {
    uint64_t querySignature;
    uint32_t targetId;
    double mitScore;
    double cfdScore;
};

/** Self-contained candidate record produced by isslExtractCandidates. */
struct HydratedCandidate {
    uint64_t signature;
    uint32_t globalId;
    uint32_t occurrences;
};

static_assert(sizeof(HydratedCandidate) == 16,
              "HydratedCandidate must match Python struct format <QII");

vector<uint8_t> nucleotideIndex(256);
vector<char> signatureIndex(4);
enum ScoreMethod { unknown = 0, mit = 1, cfd = 2, mitAndCfd = 3, mitOrCfd = 4, avgMitCfd = 5 };

/// Returns the size (bytes) of the file at `path`
size_t getFileSize(const char *path) {
    struct stat64 statBuf;
    stat64(path, &statBuf);
    return statBuf.st_size;
}

/**
 * Binary encode genetic string `ptr`
 *
 * For example, 
 *   ATCG becomes
 *   00 11 01 10  (buffer with leading zeroes to encode as 64-bit unsigned int)
 *
 * @param[in] ptr the string containing ATCG to binary encode
 */
uint64_t sequenceToSignature(const char *ptr) {
    uint64_t signature = 0;
    for (size_t j = 0; j < seqLength; j++) {
        signature |= (uint64_t)(nucleotideIndex[*ptr]) << (j * 2);
        ptr++;
    }
    return signature;
}

/**
 * Binary encode genetic string `ptr`
 *
 * For example, 
 *   00 11 01 10 becomes (as 64-bit unsigned int)
 *    A  T  C  G  (without spaces)
 *
 * @param[in] signature the binary encoded genetic string
 */
string signatureToSequence(uint64_t signature) {
    string sequence = string(seqLength, ' ');
    for (size_t j = 0; j < seqLength; j++) {
        sequence[j] = signatureIndex[(signature >> (j * 2)) & 0x3];
    }
    return sequence;
}

auto toMB = [](size_t bytes) { return bytes / (1024.0 * 1024.0); };

int main(int argc, char **argv) {
    if (argc < 9) {
        fprintf(stderr,
                "Usage: %s [scoring-metadata] [query-file] [max-distance] "
                "[score-threshold] [score-method] [candidates-file] "
                "[bucket-plan] [prefix]\n",
                argv[0]);
        exit(1);
    }

    /** Char to binary encoding */
    nucleotideIndex['A'] = 0;
    nucleotideIndex['C'] = 1;
    nucleotideIndex['G'] = 2;
    nucleotideIndex['T'] = 3;
    signatureIndex[0] = 'A';
    signatureIndex[1] = 'C';
    signatureIndex[2] = 'G';
    signatureIndex[3] = 'T';

    std::string prefixID = argv[8];
    /** The maximum number of mismatches */
    int maxDist = atoi(argv[3]);
    if (maxDist != 4) {
        fprintf(stderr, "Bucket-selective scoring requires max distance 4\n");
        return 1;
    }

    /** The threshold used to exit scoring early */
    double threshold = atof(argv[4]);

    /** The following comment below is depracated - it is left for legacy purposes */
    /** Scoring methods. To exit early: 
     *      - only CFD must drop below `threshold`
     *      - only MIT must drop below `threshold`
     *      - both CFD and MIT must drop below `threshold`
     *      - CFD or MIT must drop below `threshold`
     *      - the average of CFD and MIT must below `threshold`
     */
    string argScoreMethod = argv[5];
    ScoreMethod scoreMethod = ScoreMethod::unknown;
    bool calcCfd = false;
    bool calcMit = false;
    if (!argScoreMethod.compare("and")) {
        scoreMethod = ScoreMethod::mitAndCfd;
        calcCfd = true;
        calcMit = true;
    } else if (!argScoreMethod.compare("or")) {
        scoreMethod = ScoreMethod::mitOrCfd;
        calcCfd = true;
        calcMit = true;
    } else if (!argScoreMethod.compare("avg")) {
        scoreMethod = ScoreMethod::avgMitCfd;
        calcCfd = true;
        calcMit = true;
    } else if (!argScoreMethod.compare("mit")) {
        scoreMethod = ScoreMethod::mit;
        calcMit = true;
    } else if (!argScoreMethod.compare("cfd")) {
        scoreMethod = ScoreMethod::cfd;
        calcCfd = true;
    }

    /** Read the selected buckets assigned to this single-slice Mapper invocation. */
    const char *shardFile = argv[7];
    ifstream shardRead(shardFile);
    size_t assignedSlice, selectedBucketCount;
    int shardId;

    shardRead >> shardId >> assignedSlice >> selectedBucketCount;
    if (!shardRead) {
        fprintf(stderr, "Error reading bucket plan header\n");
        return 1;
    }

    struct SelectedBucket {
        size_t compactOffsetBytes;
        size_t elementCount;
    };
    unordered_map<size_t, SelectedBucket> selectedBuckets;
    for (size_t i = 0; i < selectedBucketCount; i++) {
        size_t bucketId, compactOffsetBytes, elementCount;
        shardRead >> bucketId >> compactOffsetBytes >> elementCount;
        if (!shardRead || !selectedBuckets.emplace(
                bucketId,
                SelectedBucket{compactOffsetBytes, elementCount}).second) {
            fprintf(stderr, "Invalid or duplicate bucket-plan entry\n");
            return 1;
        }
    }
    shardRead.close();

    /** Read only the ISSL header and precalculated MIT score table. */
    FILE *fp = fopen(argv[1], "rb");
    if (fp == nullptr) {
        fprintf(stderr, "Error opening scoring metadata\n");
        return 1;
    }

    /** The index contains a fixed-sized header 
     *      - the number of off-targets in the index
     *      - the length of an off-target
     *      - 
     *      - chars per slice
     *      - the number of slices per sequence
     *      - the number of precalculated MIT scores
     */
    vector<size_t> slicelistHeader(6);

    if (fread(slicelistHeader.data(), sizeof(size_t), slicelistHeader.size(), fp) == 0) {
        fprintf(stderr, "Error reading index: header invalid\n");
        return 1;
    }

    offtargetsCount = slicelistHeader[0];
    seqLength = slicelistHeader[1];
    seqCount = slicelistHeader[2];
    sliceWidth = slicelistHeader[3];
    sliceCount = slicelistHeader[4];
    scoresCount = slicelistHeader[5];

    /** The maximum number of possibly slice identities
     *      4 chars per slice * each of A,T,C,G = limit of 16
     */
    size_t sliceLimit = 1 << sliceWidth;
    if (sliceCount != 5 || assignedSlice >= sliceCount) {
        fprintf(stderr, "Unsupported slice configuration in bucket plan\n");
        return 1;
    }

    /** Read in the precalculated MIT scores 
     *      - `mask` is a 2-bit encoding of mismatch positions
     *          For example,
     *              00 01 01 00 01  indicates mismatches in positions 1, 3 and 4
     *  
     *      - `score` is the local MIT score for this mismatch combination
     */
    phmap::flat_hash_map<uint64_t, double> precalculatedScores;

    for (int i = 0; i < scoresCount; i++) {
        uint64_t mask = 0;
        double score = 0.0;
        if (fread(&mask, sizeof(uint64_t), 1, fp) != 1 ||
            fread(&score, sizeof(double), 1, fp) != 1) {
            fprintf(stderr, "Error reading precalculated score metadata\n");
            fclose(fp);
            return 1;
        }

        precalculatedScores.insert(pair<uint64_t, double>(mask, score));
    }

    fclose(fp);

    /** Map the concatenated, self-contained hydrated candidate records. */
    const char *candidateFile = argv[6];
    size_t candidateFileSize = getFileSize(candidateFile);
    if (candidateFileSize % sizeof(HydratedCandidate) != 0) {
        fprintf(stderr, "Candidate file is not a multiple of 16 bytes\n");
        return 1;
    }
    FILE *candidateFp = fopen(candidateFile, "rb");
    if (candidateFp == nullptr) {
        fprintf(stderr, "Error opening candidate file\n");
        return 1;
    }
    void *candidateMapping = nullptr;
    if (candidateFileSize != 0) {
        candidateMapping = mmap(nullptr, candidateFileSize, PROT_READ, MAP_PRIVATE,
                                fileno(candidateFp), 0);
        if (candidateMapping == MAP_FAILED) {
            fprintf(stderr, "Error mapping candidate file\n");
            fclose(candidateFp);
            return 1;
        }
    }
    const HydratedCandidate *candidateRecords =
        static_cast<const HydratedCandidate *>(candidateMapping);
    vector<const HydratedCandidate *> bucketPointers(sliceLimit, nullptr);
    vector<size_t> selectedBucketSizes(sliceLimit, 0);
    vector<bool> bucketAvailable(sliceLimit, false);
    for (const auto &entry : selectedBuckets) {
        size_t bucketId = entry.first;
        size_t compactOffsetBytes = entry.second.compactOffsetBytes;
        size_t elementCount = entry.second.elementCount;
        if (bucketId >= sliceLimit ||
            compactOffsetBytes % sizeof(HydratedCandidate) != 0 ||
            compactOffsetBytes > candidateFileSize ||
            elementCount >
                (candidateFileSize - compactOffsetBytes) / sizeof(HydratedCandidate)) {
            fprintf(stderr, "Bucket plan points outside candidate payload\n");
            return 1;
        }
        selectedBucketSizes[bucketId] = elementCount;
        bucketAvailable[bucketId] = true;
        bucketPointers[bucketId] = entry.second.elementCount == 0
            ? nullptr
            : candidateRecords +
                (compactOffsetBytes / sizeof(HydratedCandidate));
    }

    /** Load query file (candidate guides)
     *      and prepare memory for calculated global scores
     */
    size_t seqLineLength = seqLength + 1;
    size_t fileSize = getFileSize(argv[2]);
    if (fileSize % seqLineLength != 0) {
        fprintf(stderr, "Error: query file is not a multiple of the expected line length (%zu)\n",
                seqLineLength);
        fprintf(stderr, "The sequence length may be incorrect; alternatively, the line endings\n");
        fprintf(stderr,
                "may be something other than LF, or there may be junk at the end of the file.\n");
        exit(1);
    }
    size_t queryCount = fileSize / seqLineLength;
    fp = fopen(argv[2], "rb");
    vector<char> queryDataSet(fileSize);
    vector<uint64_t> querySignatures(queryCount);

    if (fread(queryDataSet.data(), fileSize, 1, fp) < 1) {
        fprintf(stderr, "Failed to read in query file.\n");
        exit(1);
    }
    fclose(fp);

/** Binary encode query sequences */
#pragma omp parallel
    {
#pragma omp for
        for (size_t i = 0; i < queryCount; i++) {
            char *ptr = &queryDataSet[i * seqLineLength];
            uint64_t signature = sequenceToSignature(ptr);
            querySignatures[i] = signature;
        }
    }

/** Begin scoring */
#pragma omp parallel
    {
        // Give each thread a separate output stream to avoid synchronised writes.
        int tid = omp_get_thread_num();

        std::string tmpFile = "mapper_" + prefixID + "_thread_" + std::to_string(tid) + ".bin";

        std::ofstream threadOut(tmpFile, std::ios::binary | std::ios::out);

        // Buffer records so each thread writes in chunks with less I/O overhead.
        std::vector<MapperResult> writeBuffer;
        writeBuffer.reserve(BUF_SIZE); // May need to change based on testing

/** For each candidate guide */
#pragma omp for
        for (size_t searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
            bool emit = false;
            auto searchSignature = querySignatures[searchIdx];

            /** Global scores */
            double totScoreMit = 0.0;
            double totScoreCfd = 0.0;

            int numOffTargetSitesScored = 0;
            double maximum_sum = (10000.0 - threshold * 100) / threshold;
            // bool checkNextSlice = true;

            /** Process only the semantic slice assigned to this Mapper. */
            {
                size_t i = assignedSlice;
                uint64_t sliceMask = sliceLimit - 1;
                int sliceShift = sliceWidth * i;
                sliceMask = sliceMask << sliceShift;

                uint64_t searchSlice = (searchSignature & sliceMask) >> sliceShift;
                if (!bucketAvailable[searchSlice]) {
                    fprintf(stderr, "Query requires a bucket absent from the bucket plan\n");
                    exit(1);
                }
                size_t signaturesInSlice = selectedBucketSizes[searchSlice];
                const HydratedCandidate *sliceOffset = bucketPointers[searchSlice];

                /** For each off-target signature in slice */
                for (size_t j = 0; j < signaturesInSlice; j++) {
                    const HydratedCandidate &candidate = sliceOffset[j];
                    uint64_t offtargetSignature = candidate.signature;
                    uint32_t signatureId = candidate.globalId;
                    uint32_t occurrences = candidate.occurrences;

                    /** Find the positions of mismatches 
                     *
                     *  Search signature (SS):    A  A  T  T    G  C  A  T
                     *                           00 00 11 11   10 01 00 11
                     *              
                     *        Off-target (OT):    A  T  A  T    C  G  A  T
                     *                           00 11 00 11   01 10 00 11
                     *                           
                     *                SS ^ OT:   00 00 11 11   10 01 00 11
                     *                         ^ 00 11 00 11   01 10 00 11
                     *                  (XORd) = 00 11 11 00   11 11 00 00
                     *
                     *        XORd & evenBits:   00 11 11 00   11 11 00 00
                     *                         & 10 10 10 10   10 10 10 10
                     *                   (eX)  = 00 10 10 00   10 10 00 00
                     *
                     *         XORd & oddBits:   00 11 11 00   11 11 00 00
                     *                         & 01 01 01 01   01 01 01 01
                     *                   (oX)  = 00 01 01 00   01 01 00 00
                     *
                     *         (eX >> 1) | oX:   00 01 01 00   01 01 00 00 (>>1)
                     *                         | 00 01 01 00   01 01 00 00
                     *            mismatches   = 00 01 01 00   01 01 00 00
                     *
                     *   popcount(mismatches):   4
                     */
                    uint64_t xoredSignatures = searchSignature ^ offtargetSignature;
                    uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAull;
                    uint64_t oddBits = xoredSignatures & 0x5555555555555555ull;
                    uint64_t mismatches = (evenBits >> 1) | oddBits;
                    int dist = __builtin_popcountll(mismatches);
                    // --- Distance is computed per unique (candidate guide (full sequence), off-target id) pair
                    if (dist >= 0 && dist <= maxDist) {

                        double mitContribution = 0.0;
                        double cfdContribution = 0.0;

                        {
                            // Begin calculating MIT score
                            if (calcMit) {
                                if (dist > 0) {
                                    mitContribution =
                                        precalculatedScores[mismatches] * (double)occurrences;
                                    totScoreMit += mitContribution;
                                }
                            }

                            // Begin calculating CFD score
                            if (calcCfd) {
                                /** "In other words, for the CFD score, a value of 0 
								 *      indicates no predicted off-target activity whereas 
								 *      a value of 1 indicates a perfect match"
								 *      John Doench, 2016. 
								 *      https://www.nature.com/articles/nbt.3437
								*/
                                double cfdScore = 0;
                                if (dist == 0) {
                                    cfdScore = 1;
                                } else if (dist > 0 && dist <= maxDist) {
                                    cfdScore = cfdPamPenalties
                                        [0b1010]; // PAM: NGG, TODO: do not hard-code the PAM

                                    for (size_t pos = 0; pos < 20; pos++) {
                                        size_t mask = pos << 4;

                                        /** Create the mask to look up the position-identity score
										 *      In Python... c2b is char to bit
										 *       mask = pos << 4
										 *       mask |= c2b[sgRNA[pos]] << 2
										 *       mask |= c2b[revcom(offTaret[pos])]
										 *      
										 *      Find identity at `pos` for search signature
										 *      example: find identity in pos=2
										 *       Recall ISSL is inverted, hence:
										 *                   3'-  T  G  C  C  G  A -5'
										 *       start           11 10 01 01 10 00   
										 *       3UL << pos*2    00 00 00 11 00 00 
										 *       and             00 00 00 01 00 00
										 *       shift           00 00 00 00 01 00
										 */
                                        uint64_t searchSigIdentityPos = searchSignature;
                                        searchSigIdentityPos &= (3UL << (pos * 2));
                                        searchSigIdentityPos = searchSigIdentityPos >> (pos * 2);
                                        searchSigIdentityPos = searchSigIdentityPos << 2;

                                        /** Find identity at `pos` for offtarget
										 *      Example: find identity in pos=2
										 *      Recall ISSL is inverted, hence:
										 *                  3'-  T  G  C  C  G  A -5'
										 *      start           11 10 01 01 10 00   
										 *      3UL<<pos*2      00 00 00 11 00 00 
										 *      and             00 00 00 01 00 00
										 *      shift           00 00 00 00 00 01
										 *      rev comp 3UL    00 00 00 00 00 10 (done below)
										 */
                                        uint64_t offtargetIdentityPos = offtargetSignature;
                                        offtargetIdentityPos &= (3UL << (pos * 2));
                                        offtargetIdentityPos = offtargetIdentityPos >> (pos * 2);

                                        /** Complete the mask
										 *      reverse complement (^3UL) `offtargetIdentityPos` here
										 */
                                        mask = (mask | searchSigIdentityPos |
                                                (offtargetIdentityPos ^ 3UL));

                                        if (searchSigIdentityPos >> 2 != offtargetIdentityPos) {
                                            cfdScore *= cfdPosPenalties[mask];
                                        }
                                    }
                                }
                                cfdContribution = cfdScore * (double)occurrences;
                                totScoreCfd += cfdContribution;
                            }

                            MapperResult result;
                            result.querySignature = searchSignature;
                            result.targetId = (uint32_t)signatureId; // caution check
                            result.mitScore = mitContribution;
                            result.cfdScore = cfdContribution;

                            writeBuffer.push_back(result);
                            emit = true;
                            // Flush this thread's buffered records to its temporary file.
                            if (writeBuffer.size() >= BUF_SIZE) {
                                threadOut.write(reinterpret_cast<char *>(writeBuffer.data()),
                                                writeBuffer.size() * sizeof(MapperResult));

                                writeBuffer.clear();
                            }
                            numOffTargetSitesScored += occurrences;
                        }
                    }
                }
            }
            if (!emit) {
                MapperResult q;
                q.querySignature = searchSignature;
                q.targetId = UINT32_MAX;
                q.mitScore = 0.0;
                q.cfdScore = 0.0;
                writeBuffer.push_back(q);
            }
        }

        // Flush any records remaining after the final iteration.
        if (!writeBuffer.empty()) {
            threadOut.write(reinterpret_cast<char *>(writeBuffer.data()),
                            writeBuffer.size() * sizeof(MapperResult));
            writeBuffer.clear();
        }

        threadOut.flush();
        threadOut.close();
    }

    // Create the combined output file for this shard.
    std::string shardFILEOUT = prefixID + "_shard_" + std::to_string(shardId) + ".bin";
    std::ofstream shardOut(shardFILEOUT, std::ios::binary | std::ios::out);

    int maxThreads = omp_get_max_threads();

    /** Merge and remove the temporary files created by worker threads. */
    for (int tid = 0; tid < maxThreads; tid++) {
        std::string tmpFile =
            std::string("mapper_") + prefixID + "_thread_" + std::to_string(tid) + ".bin";
        std::ifstream in(tmpFile, std::ios::binary);

        shardOut << in.rdbuf();
        in.close();
        remove(tmpFile.c_str());
    }
    shardOut.close();

    cout << "\n=== Memory Usage ===\n";
    cout << "hydratedCandidates: " << toMB(candidateFileSize) << " MB\n";
    cout << "queryDataSet:  " << toMB(queryDataSet.size()) * sizeof(char) << " MB\n";
    cout << "====================\n";
    cout << "Mapper finished, results written to: " << shardFILEOUT << endl;

    if (candidateFileSize != 0) {
        munmap(candidateMapping, candidateFileSize);
    }
    fclose(candidateFp);
    return 0;
}
