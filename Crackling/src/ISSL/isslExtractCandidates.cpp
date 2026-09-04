/*
Hydrates one raw ISSL bucket for one contiguous global off-target ID interval.

Input contracts (little-endian on the Lambda x86_64 build):
  catalogue-part.bin: headerless uint64_t signatures for [startId, endId)
  bucket.bin:         headerless uint64_t entries [occurrences:32][globalId:32]

Output contract:
  candidates.bin: headerless HydratedCandidate records compatible with <QII

Build:
  g++ -o extractor isslExtractCandidates.cpp -O3 -std=c++11
*/

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

struct HydratedCandidate {
    uint64_t signature;
    uint32_t globalId;
    uint32_t occurrences;
};

static_assert(sizeof(HydratedCandidate) == 16,
              "HydratedCandidate must match Python struct format <QII");

static bool parseUint32(const char *text, uint32_t &value) {
    if (text == nullptr || *text == '\0' || *text == '-') return false;
    errno = 0;
    char *end = nullptr;
    unsigned long long parsed = std::strtoull(text, &end, 10);
    if (errno != 0 || *end != '\0' ||
        parsed > std::numeric_limits<uint32_t>::max()) return false;
    value = static_cast<uint32_t>(parsed);
    return true;
}

static bool fileSize(FILE *file, size_t &size) {
    struct stat statBuffer;
    if (fstat(fileno(file), &statBuffer) != 0 || statBuffer.st_size < 0) return false;
    size = static_cast<size_t>(statBuffer.st_size);
    return true;
}

int main(int argc, char **argv) {
    if (argc != 6) {
        std::fprintf(stderr,
            "Usage: %s <catalogue-part.bin> <bucket.bin> <start-id> <end-id> <output.bin>\n",
            argv[0]);
        return 1;
    }

    uint32_t startId = 0, endId = 0;
    if (!parseUint32(argv[3], startId) || !parseUint32(argv[4], endId) ||
        endId < startId) {
        std::fprintf(stderr, "Invalid half-open global ID interval\n");
        return 1;
    }

    FILE *catalogueFile = std::fopen(argv[1], "rb");
    FILE *bucketFile = std::fopen(argv[2], "rb");
    if (catalogueFile == nullptr || bucketFile == nullptr) {
        std::fprintf(stderr, "Unable to open extractor input file\n");
        if (catalogueFile) std::fclose(catalogueFile);
        if (bucketFile) std::fclose(bucketFile);
        return 1;
    }

    size_t catalogueBytes = 0, bucketBytes = 0;
    if (!fileSize(catalogueFile, catalogueBytes) ||
        !fileSize(bucketFile, bucketBytes) ||
        catalogueBytes % sizeof(uint64_t) != 0 ||
        bucketBytes % sizeof(uint64_t) != 0 ||
        catalogueBytes / sizeof(uint64_t) !=
            static_cast<uint64_t>(endId) - static_cast<uint64_t>(startId)) {
        std::fprintf(stderr, "Extractor input size does not match its record contract\n");
        std::fclose(catalogueFile);
        std::fclose(bucketFile);
        return 1;
    }

    const uint64_t *catalogue = nullptr;
    const uint64_t *bucket = nullptr;
    void *catalogueMapping = nullptr, *bucketMapping = nullptr;
    if (catalogueBytes != 0) {
        catalogueMapping = mmap(nullptr, catalogueBytes, PROT_READ, MAP_PRIVATE,
                                fileno(catalogueFile), 0);
        if (catalogueMapping == MAP_FAILED) catalogueMapping = nullptr;
        catalogue = static_cast<const uint64_t *>(catalogueMapping);
    }
    if (bucketBytes != 0) {
        bucketMapping = mmap(nullptr, bucketBytes, PROT_READ, MAP_PRIVATE,
                            fileno(bucketFile), 0);
        if (bucketMapping == MAP_FAILED) bucketMapping = nullptr;
        bucket = static_cast<const uint64_t *>(bucketMapping);
    }
    if ((catalogueBytes != 0 && catalogueMapping == nullptr) ||
        (bucketBytes != 0 && bucketMapping == nullptr)) {
        std::fprintf(stderr, "Unable to memory-map extractor input\n");
        std::fclose(catalogueFile);
        std::fclose(bucketFile);
        return 1;
    }

    std::ofstream output(argv[5], std::ios::binary | std::ios::trunc);
    if (!output) {
        std::fprintf(stderr, "Unable to create candidate output\n");
        return 1;
    }

    const size_t bucketCount = bucketBytes / sizeof(uint64_t);
    uint64_t previousId = 0;
    bool havePrevious = false;
    size_t written = 0;
    for (size_t i = 0; i < bucketCount; ++i) {
        uint64_t packed = bucket[i];
        uint32_t globalId = static_cast<uint32_t>(packed & 0xffffffffULL);
        if (havePrevious && globalId <= previousId) {
            std::fprintf(stderr, "Bucket global IDs are not strictly increasing\n");
            return 1;
        }
        previousId = globalId;
        havePrevious = true;

        if (globalId < startId) continue;
        if (globalId >= endId) break;

        HydratedCandidate candidate;
        candidate.signature = catalogue[static_cast<size_t>(globalId - startId)];
        candidate.globalId = globalId;
        candidate.occurrences = static_cast<uint32_t>(packed >> 32);
        output.write(reinterpret_cast<const char *>(&candidate), sizeof(candidate));
        if (!output) {
            std::fprintf(stderr, "Failed while writing candidate output\n");
            return 1;
        }
        ++written;
    }

    output.close();
    if (catalogueBytes != 0) munmap(catalogueMapping, catalogueBytes);
    if (bucketBytes != 0) munmap(bucketMapping, bucketBytes);
    std::fclose(catalogueFile);
    std::fclose(bucketFile);

    std::cout << "candidateRecords=" << written << "\n";
    return 0;
}
