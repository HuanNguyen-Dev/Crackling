
/*

Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


To compile:

g++ -o isslScoreOfftargets isslScoreOfftargets.cpp -O3 -std=c++11 -fopenmp -mpopcnt -Iparallel_hashmap
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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/time.h>
#include <chrono>
#include <bitset>
#include <iostream>
#include <climits>
#include <stdio.h>
#include <cstring>
#include <omp.h>
#include <phmap.h>
#include <map>

#include <fstream>

using namespace std;

#define BUF_SIZE 1000

size_t seqLength, seqCount, sliceWidth, sliceCount, offtargetsCount, scoresCount;
struct MapperResult{
    uint64_t querySignature;
    uint32_t targetId; // -- Ignore number of occurences for the results
    double mitScore;
    double cfdScore; // -- Is this amount of precision needed? Float halves the space
};


vector<uint8_t> nucleotideIndex(256);
vector<char> signatureIndex(4);
enum ScoreMethod { unknown = 0, mit = 1, cfd = 2, mitAndCfd = 3, mitOrCfd = 4, avgMitCfd = 5 };

/// Returns the size (bytes) of the file at `path`
size_t getFileSize(const char *path)
{
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
uint64_t sequenceToSignature(const char *ptr)
{
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
string signatureToSequence(uint64_t signature)
{
    string sequence = string(seqLength, ' ');
    for (size_t j = 0; j < seqLength; j++) {
        sequence[j] = signatureIndex[(signature >> (j * 2)) & 0x3];
    }
    return sequence;
}

auto toMB = [](size_t bytes) {
    return bytes / (1024.0 * 1024.0);
};

int main(int argc, char **argv)
{
    if (argc < 8) {
        fprintf(stderr, "Usage: %s [issltable] [query file] [max distance] [score-threshold] [score-method] [shard-file] [prefix]\n", argv[0]);
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

    std::string prefixID = argv[7];
    /** The maximum number of mismatches */
    int maxDist = atoi(argv[3]);
    
    /** The threshold used to exit scoring early */
    double threshold = atof(argv[4]);
    
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
	
    /** New - Retrieve the start and end bytes of the slices this mapper will process */
    const char* shardFile= argv[6];
    ifstream shardRead(shardFile);
    size_t startSlice, endSlice, startByte, endByte;
    int shardId;
    // might not need the rest here
    // size_t offtargetsCount, seqLength, seqCount, sliceWidth, sliceCount, scoresCount;

    shardRead >> shardId >> startSlice >> endSlice >> startByte >> endByte;
    // shardRead >> offtargetsCount >> seqLength >> seqCount >> sliceWidth >> sliceCount >> scoresCount
    shardRead.close();

    /** Begin reading the binary encoded ISSL, structured as:
     *      - a header (6 items)
     *      - precalcuated local MIT scores
     *      - all binary-encoded off-target sites
     *      - slice list sizes
     *      - slice contents
     */
    FILE *fp = fopen(argv[1], "rb");
    
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
    seqLength       = slicelistHeader[1]; 
    seqCount        = slicelistHeader[2]; 
    sliceWidth      = slicelistHeader[3]; 
    sliceCount      = slicelistHeader[4]; 
    scoresCount     = slicelistHeader[5]; 
    
    cout << "offtargetsCount = " << offtargetsCount << endl;
    /** The maximum number of possibly slice identities
     *      4 chars per slice * each of A,T,C,G = limit of 16
     */
    size_t sliceLimit = 1 << sliceWidth;
    
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
        fread(&mask, sizeof(uint64_t), 1, fp);
        fread(&score, sizeof(double), 1, fp);
        
        precalculatedScores.insert(pair<uint64_t, double>(mask, score));
    }
    
    /** Load in all of the off-target sites */
    vector<uint64_t> offtargets(offtargetsCount);
    if (fread(offtargets.data(), sizeof(uint64_t), offtargetsCount, fp) == 0) {
        fprintf(stderr, "Error reading index: loading off-target sequences failed\n");
        return 1;
    }
    
    /** Prevent assessing an off-target site for multiple slices
     *
     *      Create enough 1-bit "seen" flags for the off-targets
     *      We only want to score a candidate guide against an off-target once.
     *      The least-significant bit represents the first off-target
     *      0 0 0 1   0 1 0 0   would indicate that the 3rd and 5th off-target have been seen.
     *      The CHAR_BIT macro tells us how many bits are in a byte (C++ >= 8 bits per byte)
     */
    uint64_t numOfftargetToggles = (offtargetsCount / ((size_t)sizeof(uint64_t) * (size_t)CHAR_BIT)) + 1;

    /** The number of signatures embedded per slice
     *
     *      These counts are stored contiguously
     *
     */
    vector<size_t> allSlicelistSizes(sliceCount * sliceLimit);
    
    if (fread(allSlicelistSizes.data(), sizeof(size_t), allSlicelistSizes.size(), fp) == 0) {
        fprintf(stderr, "Error reading index: reading slice list sizes failed\n");
        return 1;
    }
    
    /** New - The contents of the slices from the sharded file
     *
     *      Stored contiguously
     *
     *      Each signature (64-bit) is structured as:
     *          <occurrences 32-bit><off-target-id 32-bit>
     */
    // vector<uint64_t> allSignatures(seqCount * sliceCount);

    // Skip to the start of the slice to process
    fseek(fp, startByte, SEEK_SET);

    size_t sliceBytes = endByte - startByte;
    size_t numElements = sliceBytes / sizeof(uint64_t);

    vector<uint64_t> sliceSignatureBuffer(numElements);

    if (fread(sliceSignatureBuffer.data(), sizeof(uint64_t), numElements, fp) != numElements) {
        fprintf(stderr, "Error reading index: reading slice contents failed, Could not read all %zu elements\n", numElements);
        return 1;
    }
    
    /** End reading the index */
    fclose(fp);
    size_t localSliceCount = endSlice - startSlice;
    /** Start constructing index in memory
     *
     *      To begin, reverse the contiguous storage of the slices,
     *         into the following:
     *
     *         + Slice 0 :
     *         |---- AAAA : <slice contents>
     *         |---- AAAC : <slice contents>
     *         |----  ...
     *         | 
     *         + Slice 1 :
     *         |---- AAAA : <slice contents>
     *         |---- AAAC : <slice contents>
     *         |---- ...
     *         | ...
     */
    // Note: To Do: this should be using a new slice count
    vector<vector<uint64_t *>> sliceLists(localSliceCount, vector<uint64_t *>(sliceLimit));

    uint64_t *offset = sliceSignatureBuffer.data();
 
    for (size_t i = startSlice; i < endSlice; i++) {
        size_t local_i = i - startSlice;
        for (size_t j = 0; j < sliceLimit; j++) {
            size_t idx = i * sliceLimit + j;
            sliceLists[local_i][j] = offset;
            offset += allSlicelistSizes[idx];
        }
    }   
    if (offset != sliceSignatureBuffer.data() + numElements)
    {
        fprintf(stderr, "Shard reconstruction mismatch\n");
    }
    
    /** Load query file (candidate guides)
     *      and prepare memory for calculated global scores
     */
    size_t seqLineLength = seqLength + 1;
    size_t fileSize = getFileSize(argv[2]);
    if (fileSize % seqLineLength != 0) {
        fprintf(stderr, "Error: query file is not a multiple of the expected line length (%zu)\n", seqLineLength);
        fprintf(stderr, "The sequence length may be incorrect; alternatively, the line endings\n");
        fprintf(stderr, "may be something other than LF, or there may be junk at the end of the file.\n");
        exit(1);
    }
    size_t queryCount = fileSize / seqLineLength;
    fp = fopen(argv[2], "rb");
    vector<char> queryDataSet(fileSize);
    vector<uint64_t> querySignatures(queryCount);
    // vector<double> querySignatureMitScores(queryCount);
    // vector<double> querySignatureCfdScores(queryCount);

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
        // Create temporary binary stream for each thread
        int tid = omp_get_thread_num();

        std::string tmpFile = "mapper_" +prefixID + 
                                "_thread_" + std::to_string(tid) + ".bin";


        std::ofstream threadOut(
            tmpFile,
            std::ios::binary | std::ios::out
        );
        
        // Buffer for thread to write to file in chunks for less IO overhead
        std::vector<MapperResult> writeBuffer;
        // 32 Bytes / record -> 3.2 MB/Thread for 100000  -> 8 bytes + 4 + 8 + 8 = 28 bytes + 4 padding to get it divisible by 8
        writeBuffer.reserve(BUF_SIZE); // May need to change based on testing

        // searchResults is unused ??
        // unordered_map<uint64_t, unordered_set<uint64_t>> searchResults;
        vector<uint64_t> offtargetToggles(numOfftargetToggles);
    
        uint64_t * offtargetTogglesTail = offtargetToggles.data() + numOfftargetToggles - 1;

        /** For each candidate guide */
        #pragma omp for
        for (size_t searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
            bool emit = false;
            // --- Represents the candidate guide sequence
            auto searchSignature = querySignatures[searchIdx];

            /** Global scores */
            double totScoreMit = 0.0;
            double totScoreCfd = 0.0;
            
            int numOffTargetSitesScored = 0;
            double maximum_sum = (10000.0 - threshold*100) / threshold;
            // bool checkNextSlice = true;
            
            /** For each ISSL slice */
            /** We want to use the global i as bit extraction depends on absolute position
             *  and the 'slicelistsizes' are laid out for all slices globally
             */
            for (size_t i = startSlice; i < endSlice; i++) {
                // --- build a mask for one slice - given slice limit = 16, slice mask 4 
                // --- bits (e.g 2 bases -> (A,C,G,T) encoded into binary would be something like to AG -> 0011 which is 4 bits)
                uint64_t sliceMask = sliceLimit - 1;
                // --- We want to go down the query signature sequence for every mask
                int sliceShift = sliceWidth * i;
                // --- This essentially targets x bits of the sequence 
                // --- e.g is slice shift = 8, then the mask targets bits 8-11
                sliceMask = sliceMask << sliceShift;
                auto &sliceList = sliceLists[i - startSlice];
                
                // --- Extract the slice from the query signature sequence and shift it back down
                // --- to a value between 0-15
                uint64_t searchSlice = (searchSignature & sliceMask) >> sliceShift;
                
                // --- This is the index for the contiguously stored slicelist sizes
                size_t idx = i * sliceLimit + searchSlice;
                // --- Get how many sequences are in the bucket
                size_t signaturesInSlice = allSlicelistSizes[idx];
                // --- Look up the bucket. This means that all sequences have this slice value at position i
                uint64_t *sliceOffset = sliceList[searchSlice];
                
                // --- We now know where the list starts (sliceoffset) and how many entries it has 
                // --- signatures in slice
                // --- For each i, we extract slice value from query, use it as an index, jump directly to matching 
                // --- sequences

                /** For each off-target signature in slice */
                for (size_t j = 0; j < signaturesInSlice; j++) {
                    // --- [32-bit occurances | 32-bit signatureId]
                    // --- equivalent to off-target IDs
                    auto signatureWithOccurrencesAndId = sliceOffset[j];
                    auto signatureId = signatureWithOccurrencesAndId & 0xFFFFFFFFull;
                    uint32_t occurrences = (signatureWithOccurrencesAndId >> (32));

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
                    uint64_t xoredSignatures = searchSignature ^ offtargets[signatureId];
                    uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAull;
                    uint64_t oddBits = xoredSignatures & 0x5555555555555555ull;
                    uint64_t mismatches = (evenBits >> 1) | oddBits;
                    int dist = __builtin_popcountll(mismatches);
					// --- Distance is computed per unique (candidate guide (full sequence), off-target id) pair
					if (dist >= 0 && dist <= maxDist) {

						/** Prevent assessing the same off-target for multiple slices */
						uint64_t seenOfftargetAlready = 0;
						uint64_t * ptrOfftargetFlag = (offtargetTogglesTail - (signatureId / 64));
						seenOfftargetAlready = (*ptrOfftargetFlag >> (signatureId % 64)) & 1ULL;
                        double mitContribution = 0.0;
                        double cfdContribution = 0.0;
										
						if (!seenOfftargetAlready) {
							// Begin calculating MIT score
							if (calcMit) {
								if (dist > 0) {
                                    mitContribution = precalculatedScores[mismatches] * (double)occurrences;
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
								}
								else if (dist > 0 && dist <= maxDist) {
									cfdScore = cfdPamPenalties[0b1010]; // PAM: NGG, TODO: do not hard-code the PAM
									
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
										uint64_t offtargetIdentityPos = offtargets[signatureId];
										offtargetIdentityPos &= (3UL << (pos * 2));
										offtargetIdentityPos = offtargetIdentityPos >> (pos * 2); 

										/** Complete the mask
										 *      reverse complement (^3UL) `offtargetIdentityPos` here
										 */
										mask = (mask | searchSigIdentityPos | (offtargetIdentityPos ^ 3UL));

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
                            result.targetId = (uint32_t) signatureId; // caution check
                            result.mitScore = mitContribution;
                            result.cfdScore = cfdContribution;

                            writeBuffer.push_back(result);
                            emit = true;
                            // Flush the chunk held in the buffer
                            if (writeBuffer.size() >= BUF_SIZE) {
                                threadOut.write(
                                    reinterpret_cast<char*>(writeBuffer.data()),
                                    writeBuffer.size() * sizeof(MapperResult)
                                );

                                writeBuffer.clear();
                            }
							*ptrOfftargetFlag |= (1ULL << (signatureId % 64));
							numOffTargetSitesScored += occurrences;

							/** Stop calculating global score early if possible */
				// 			if (scoreMethod == ScoreMethod::mitAndCfd) {
				// 				if (totScoreMit > maximum_sum && totScoreCfd > maximum_sum) {
				// 					checkNextSlice = false;
				// 					break;
				// 				}
				// 			}
				// 			if (scoreMethod == ScoreMethod::mitOrCfd) {
				// 				if (totScoreMit > maximum_sum || totScoreCfd > maximum_sum) {
				// 					checkNextSlice = false;
				// 					break;
				// 				}
				// 			}
				// 			if (scoreMethod == ScoreMethod::avgMitCfd) {
				// 				if (((totScoreMit + totScoreCfd) / 2.0) > maximum_sum) {
				// 					checkNextSlice = false;
				// 					break;
				// 				}
				// 			}
				// 			if (scoreMethod == ScoreMethod::mit) {
				// 				if (totScoreMit > maximum_sum) {
				// 					checkNextSlice = false;
				// 					break;
				// 				}
				// 			}
				// 			if (scoreMethod == ScoreMethod::cfd) {
				// 				if (totScoreCfd > maximum_sum) {
				// 					checkNextSlice = false;
				// 					break;
				// 				}
				// 			}
						}
                    }
                }

                // if (!checkNextSlice)
                //     break;
            }
            if (!emit){
                MapperResult q;
                q.querySignature = searchSignature;
                q.targetId = UINT32_MAX;
                q.mitScore = 0.0;
                q.cfdScore = 0.0;
                writeBuffer.push_back(q);
            }
            // probably do not need this
            // querySignatureMitScores[searchIdx] = 10000.0 / (100.0 + totScoreMit);
            // querySignatureCfdScores[searchIdx] = 10000.0 / (100.0 + totScoreCfd);
            // MapperResult q;
            // q.querySignature = searchSignature;
            // q.targetId = UINT32_MAX;
            // q.mitScore = 10000.0 / (100.0 + totScoreMit);
            // q.cfdScore = 10000.0 / (100.0 + totScoreCfd);
            // writeBuffer.push_back(q);
            memset(offtargetToggles.data(), 0, sizeof(uint64_t)*offtargetToggles.size());
        }
        // This preserves: “queries outside threshold still appear”
        // MapperResult q;
        // q.querySignature = searchSignature;
        // q.targetId = UINT32_MAX;
        // q.mitScore = 10000.0 / (100.0 + totScoreMit);
        // q.cfdScore = 10000.0 / (100.0 + totScoreCfd);

        // writeBuffer.push_back(q);

        // Empty the buffer out at the last iteration
        if (!writeBuffer.empty()) {
            threadOut.write(
                reinterpret_cast<char*>(writeBuffer.data()),
                writeBuffer.size() * sizeof(MapperResult)
            );
            writeBuffer.clear();
        }

        threadOut.flush();
        threadOut.close();

    }

    // New - Print output to mapper file
    std::string shardFILEOUT = prefixID + "_shard_" + std::to_string(shardId) + ".bin";
    std::ofstream shardOut(
        shardFILEOUT,
        std::ios::binary | std::ios::out
    );

    int maxThreads = omp_get_max_threads();


    /** Merge temporary files created by threads */
    // To Do: Parallelise this for the future
    for (int tid = 0; tid < maxThreads; tid ++) {
        std::string tmpFile =
        std::string("mapper_") + prefixID + "_thread_" + std::to_string(tid) + ".bin";
        std::ifstream in(tmpFile, std::ios::binary);

        shardOut << in.rdbuf();
        in.close();
        remove(tmpFile.c_str());
    }
    shardOut.close();


    cout << "\n=== Memory Usage ===\n";
    cout << "allSignatures: " << toMB(sliceSignatureBuffer.size() * sizeof(uint64_t)) << " MB\n";
    cout << "queryDataSet:  " << toMB(queryDataSet.size()) * sizeof(char) << " MB\n";
    cout << "offtargets:    " << toMB(offtargets.size() * sizeof(uint64_t)) << " MB\n";
    cout << "====================\n";
    cout << "Mapper finished, results written to: " << shardFILEOUT << endl;

    return 0;
}