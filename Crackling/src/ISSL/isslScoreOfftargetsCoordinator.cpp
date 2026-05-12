/*

Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


To compile:

g++ -o coordinator isslScoreOfftargetsCoordinator.cpp 
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


using namespace std;

size_t seqLength, seqCount, sliceWidth, sliceCount, offtargetsCount, scoresCount;


int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "Usage: %s [issltable] [num-shards] [output-prefix]\n", argv[0]);
        exit(1);
    }
    
    /** The number of shards to chunk the ISSL - Assume evenly distributed */
    size_t numShards = atoi(argv[2]);

    /** The output prefix for the file being written to */
    string prefix = argv[3];
    

	
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
    /** The maximum number of possibly slice identities
     *      4 chars per slice * each of A,T,C,G = limit of 16
     */
    size_t sliceLimit = 1 << sliceWidth;
    size_t slicesPerShard = (sliceCount + numShards - 1)/numShards;
    size_t precalculatedScoresBytes = (scoresCount * 2) * sizeof(uint64_t);

    // Skip the precalculated scores
    fseek(fp, precalculatedScoresBytes, SEEK_CUR);

    size_t offTargetsBytes = (offtargetsCount * sizeof(uint64_t));
    
    // Skip the offtarget array
    fseek(fp, offTargetsBytes, SEEK_CUR);

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

    /** End reading the index */
    fclose(fp);

    // Precompute cumulative sums of bucket sizes up to but not including i
    // Total bytes before a slice = sum of all previous buckets * sizeof(size_t)
    vector<size_t> cumulativeBytes(sliceCount * sliceLimit + 1, 0);
    for (size_t i = 0; i < allSlicelistSizes.size(); i++){
        // note cumulativeBytes[0] = 0
        cumulativeBytes[i+1] = cumulativeBytes[i] + (allSlicelistSizes[i] * sizeof(uint64_t));
    }

    // Base offset in bytes from start of ISSL to just before all signatures array
    size_t baseOffsetBytes = (6 + (scoresCount * 2) + offtargetsCount + (sliceCount * sliceLimit)) * sizeof(uint64_t);

    // Determine start and end bytes for each shard
    for (int shard = 0; shard < numShards; shard ++){
        size_t startSlice = shard * slicesPerShard;
        size_t endSlice = std::min((shard + 1) * slicesPerShard, sliceCount);

        if (startSlice >= sliceCount) continue; // Skip empty shards

        size_t startByte = baseOffsetBytes + cumulativeBytes[startSlice * sliceLimit];
        
        size_t endByte = baseOffsetBytes + cumulativeBytes[endSlice * sliceLimit];

        string shardFile = prefix + "_shard_" + to_string(shard) + ".txt";
        ofstream out(shardFile);
        // might not need all that, better to keep it lean maybe
        out << shard << " "<< startSlice << " " << endSlice << " " << startByte << " " << endByte << "\n";
        out.close();

        cout << "Created shard file: " << shardFile << "\n";
        cout << "Contents: " << "Shard: " << shard << " startSlice: " << startSlice << " endSlice: " << endSlice << 
        " startByte -> endByte: " << startByte << " -> " << endByte << "\n";

    }

    return 0;
} 

// notes include a new slice count
// notes shard file should include a prefix to identify which shard
// also a job id will be required