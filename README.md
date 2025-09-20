# Streaming Masked Superstrings

## Introduction

Streaming Masked Superstrings is a tool for computing [masked
superstrings](https://www.biorxiv.org/content/10.1101/2023.02.01.526717v1) of
genomic data.

The main goal of this project is to provide a method of computing masked
superstrings with small memory overhead. This is achieved by using
probabilistic data structures to represent the set of k-mers present in the
input sequence. The trade-off of this approach is that the output is not exact.
It can incorrectly mark some k-mers as not present and represent some k-mers
more than once.

## Software requirements

To compile this project, you need to have the following software installed:
- CMake (version 3.10 or higher)
- C++ compiler that supports the C++23 standard

## Build instructions

After cloning the repository, navigate to the project directory and run the following commands:

1. Create a build directory:
```bash
mkdir build && cd build
```
2. Run CMake to configure the project:
```bash
cmake ..
```
3. Build the project:
```bash
make
```

This will create the `streaming-masked-superstring` executable in the
`build/src` directory.

## Usage

### Streaming algorithm

Compute the masked superstring with the approximate algorithm by using the
`compute` subcommand:
```bash
streaming-masked-superstring compute -k 31 -bpk 10 <input-fasta> <output-fasta> # Compute masked superstring with k-mer size 31 and 10 bits-per-kmer
streaming-masked-superstring compute -f <input-fasta> <output-fasta> # Run only the first phase of the streaming algorithm
streaming-masked-superstring compute -t tmp.fa --no-splice <input-fasta> <output-fasta> # Do not use splicing in the final output and write intermediate result to tmp.fa
```

### Exact algorithm

Compute the exact masked superstring by using the `exact` subcommand:
```bash
streaming-masked-superstring exact -k 31 <input-fasta> <output-fasta> # Compute exact masked superstring with k-mer size 31
streaming-masked-superstring exact -u --no-splice <input-fasta> <output-fasta> # Do not use splicing and ignore reverse complements
```

### Testing accuracy

To compare and report the accuracy of the approximate algorithm use the
`compare` subcommand:
```bash
streaming-masked-superstring compare -k 31 <approximate-output> <exact-output> # Compute the accuracy of the approximate output with k-mer size 31
```

To view all options for a particular subcommand, run `streaming-masked-superstring <subcommand> --help`. The maximum supported value of `k` for all subcommands is 32.

## How it works

For details about the algorithm, see the [algorithms.md](./docs/algorithms.md) document.

For details about design and implementation, see the [design.md](./docs/design.md) document.
