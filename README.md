# Crackling ISSL tools

This directory contains the C++ tools used to build an inverted signature slice
list (ISSL) index and score CRISPR guide candidates against it. It includes the
original single-process scorer and an experimental Coordinator/Mapper/Reducer
pipeline for distributing off-target scoring across shards.

The distributed Mapper is derived from the original Crackling scorer at commit
[`e247e9b`](https://github.com/bmds-lab/Crackling/blob/e247e9b48ca9ad8de04258d51295126250ef1c26/src/ISSL/isslScoreOfftargets.cpp).

## Components

| Source | Executable | Purpose |
| --- | --- | --- |
| `isslCreateIndex.cpp` | `createIsslIndex` | Builds a binary ISSL index from an off-target sequence file. |
| `isslScoreOfftargets.cpp` | `isslScoreOfftargets` | Runs the original, non-distributed off-target scorer. |
| `isslScoreOfftargetsCoordinator.cpp` | `coordinator` | Divides the ISSL slice-data region into shard descriptions. |
| `isslScoreOfftargetsMapper.cpp` | `mapper` | Scores one assigned shard and emits per-off-target MIT/CFD contributions. |
| `isslScoreOfftargetsReducer.cpp` | `reducer` | Deduplicates Mapper records and produces final guide scores. |
| `difference.cpp` | `difference` | Compares two text score outputs for validation. |
| `isslScoreOfftargetsV2.cpp` | `isslScoreOfftargetsV2` | An alternate local scorer retained for development and comparison. |

Precompiled Linux binaries may also be present. Rebuild them after changing
their corresponding source files.

## Requirements

- Linux x86_64
- A C++11-compatible `g++`
- OpenMP support
- POPCNT support on the build and runtime CPU
- Headers under `include/`, including `phmap.h` and `cfdPenalties.h`

## Building

Run these commands from this directory:

```bash
g++ -o createIsslIndex isslCreateIndex.cpp -O3 -std=c++11 -fopenmp -mpopcnt
g++ -o isslScoreOfftargets isslScoreOfftargets.cpp -O3 -std=c++11 -fopenmp -mpopcnt -Iinclude
g++ -o coordinator isslScoreOfftargetsCoordinator.cpp
g++ -o mapper isslScoreOfftargetsMapper.cpp -O3 -std=c++11 -fopenmp -mpopcnt -Iinclude
g++ -o reducer isslScoreOfftargetsReducer.cpp -std=c++17
g++ -o difference difference.cpp
g++ -o isslScoreOfftargetsV2 isslScoreOfftargetsV2.cpp -O3 -std=c++11 -fopenmp -mpopcnt -Iinclude
```

## Input conventions

Candidate-guide and off-target input files contain one fixed-length uppercase
A/C/G/T sequence per line. The programs expect LF line endings and validate
file size against `sequence length + 1`. The off-target file supplied to the
index builder must be sorted because adjacent identical sequences are counted
as occurrences. Sequence length is limited to 32 bases by the index builder.

The score method accepted by the scorers is one of:

- `mit`
- `cfd`
- `and`
- `or`
- `avg`

## Creating an ISSL index

```bash
./createIsslIndex \
  offtargetSites.txt \
  <sequence-length> \
  <slice-width-bits> \
  output.issl
```

Example:

```bash
./createIsslIndex offtargetSites.txt 20 4 genome.issl
```

The binary index contains a six-value header, precalculated MIT scores,
binary-encoded off-target sequences, slice-list sizes, and slice contents. It
uses native C++ integer sizes and byte order, so it is not a portable interchange
format across arbitrary architectures or ABIs.

## Original scorer

```bash
./isslScoreOfftargets \
  <issl-index> \
  <query-file> \
  <max-distance> \
  <score-threshold> \
  <score-method>
```

Example:

```bash
./isslScoreOfftargets genome.issl queries.txt 4 75 and
```

The original scorer writes final guide sequences and MIT/CFD scores to standard
output.

## Distributed scoring pipeline

The distributed path separates the work into three stages:

```text
ISSL index -> Coordinator -> shard descriptions -> Mappers
                                                    |
                                                    v
                                           binary shard results
                                                    |
                                                    v
                                                 Reducer
                                                    |
                                                    v
                                           final text scores
```

### 1. Coordinator

```bash
./coordinator <issl-index> <number-of-shards> <output-prefix>
```

Example:

```bash
./coordinator genome.issl 5 genome
```

For each shard, the Coordinator creates:

```text
<output-prefix>_shard_<shard-id>.txt
```

Each file contains one line:

```text
<shard-id> <start-slice> <end-slice> <start-byte> <end-byte>
```

Slice ranges use an inclusive start and exclusive end. Byte ranges identify the
assigned slice-data region in the ISSL index.

### 2. Mapper

Run one Mapper for each shard description:

```bash
./mapper \
  <issl-index> \
  <query-file> \
  <max-distance> \
  <score-threshold> \
  <score-method> \
  <shard-file> \
  <output-prefix>
```

Example:

```bash
./mapper genome.issl queries.txt 4 75 and genome_shard_0.txt result
```

The Mapper creates:

```text
<output-prefix>_shard_<shard-id>.bin
```

OpenMP threads initially write independent files named:

```text
mapper_<output-prefix>_thread_<thread-id>.bin
```

The Mapper concatenates these files into the shard output and removes them.
Because the names are relative, they are created in the process's current
working directory. The AWS Lambda wrapper launches the Mapper in a temporary
directory under `/tmp`, uploads the durable results to S3, and then removes the
temporary directory.

During the merge, the per-thread files coexist with the growing combined file,
so peak temporary-disk usage can approach twice the total Mapper output size.
Set `OMP_NUM_THREADS` to control Mapper thread count when required.

#### Mapper binary output contract

The shard output is a headerless sequence of `MapperResult` records. On the
current Linux x86_64 target, each record is 32 bytes:

| Offset | Size | Type | Field |
| ---: | ---: | --- | --- |
| 0 | 8 | `uint64_t` | `querySignature` |
| 8 | 4 | `uint32_t` | `targetId` |
| 12 | 4 | padding | Structure alignment padding |
| 16 | 8 | `double` | `mitScore` |
| 24 | 8 | `double` | `cfdScore` |

Records use the host's little-endian byte order. The AWS Lambda wrapper reads
this layout using Python's `struct` format `<QI4xdd`. A `targetId` of
`UINT32_MAX` indicates that the query emitted no off-target contributions.

Keep `MapperResult`, the Python unpacking format, and every Reducer reader in
sync. This raw structure format is architecture-dependent; changing compiler,
ABI, field order, or target architecture requires compatibility validation.

### 3. Reducer

```bash
./reducer \
  <output-file> \
  <sequence-length> \
  <mapper-output-1> \
  <mapper-output-2> ...
```

Example:

```bash
./reducer final-scores.txt 20 result_shard_0.bin result_shard_1.bin \
  result_shard_2.bin result_shard_3.bin result_shard_4.bin
```

The local C++ Reducer loads the supplied Mapper records, sorts them by query
signature and target ID, ignores `UINT32_MAX` sentinels, and deduplicates target
IDs before summing their MIT and CFD contributions. It writes tab-separated
text records:

```text
<guide-sequence>\t<final-mit-score>\t<final-cfd-score>
```

The local Reducer loads all Mapper results into memory. The AWS implementation
uses a separate Python Reducer designed for Lambda's resource constraints.

## Comparing results

Use `difference` to compare two tab-separated scorer outputs:

```bash
./difference <file-1> <file-2> <sequence-length>
```

Differences greater than `1e-4` are printed for inspection.

## AWS integration notes

The AWS deployment uses the compiled Linux x86_64 `mapper` binary through a
Lambda layer. The Lambda wrapper prepares the assigned ISSL data and query/shard
files, runs the binary in `/tmp`, splits the combined Mapper output into MIT and
CFD streams, and uploads those streams and a completion marker to S3.

The local C++ Coordinator and Reducer are useful for development, validation,
and benchmarking. They are not the binaries executed by the current AWS
Coordinator and Reducer Lambdas.

## Attribution

The scoring code is based on:

> Faster and better CRISPR guide RNA design with the Crackling method.  
> Jacob Bradford, Timothy Chappell, Dimitri Perrin.  
> bioRxiv 2020.02.14.950261. https://doi.org/10.1101/2020.02.14.950261
