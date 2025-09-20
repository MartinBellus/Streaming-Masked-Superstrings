# Programmers Documentation

## Introduction

This document serves as an entry point for understanding the architecture and
design of the Streaming‑Masked‑Superstrings project.

For algorithm overview and theoretical background, see the
[algorithms.md](./algorithms.md) document or the other references listed at
the end.

---

## Module overview

### Algorithm Module

The algorithm module contains the implementations of both the exact and
streaming algorithms for computing masked superstrings.

#### Streaming algorithm

The streaming algorithm is split into two phases, both using Bloom Filters to
achieve memory efficiency while maintaining reasonable accuracy.

The sizes of the Bloom filters are automatically determined based on the
approximate number of unique k-mers present in the input sequence and a
configurable bits-per-kmer parameter.

For more details, see the [algorithms.md](./algorithms.md) document.

#### Exact algorithm

The exact algorithm uses a `std::unordered_set` to store all present
k-mers in the sequence. It processes the input in a single pass, adding each
k-mer to the set as it is encountered.

The output of the exact algorithm can then be used to determine the accuracy of
the streaming algorithm.

### Hash Module

The Hash module currently provides implementations of two non-cryptographic hash
families: Polynomial hash and [Murmur hash][murmurhash].

For both of the hash families, we use the [double
hashing][kirsch-mitzmacker-2008] technique for better performance.

### Helper Module

The Helper module contains various utility functions and classes used across the
project.

#### `ComputeArgs`, `ExactArgs` and `CompareArgs`

These classes encapsulate the command-line arguments for the `compute`, `exact`
and `compare` subcommands, respectively.

#### `DynamicBitset`

A simple bitset implementation that holds its data on the heap. It is used by
the `BloomFilter` data structure.

#### `CountingBitset`

CountingBitset represents an array of counters. It has a template parameter
`BPC` for the number of bits per counter. The value of `BPC` should divide 32 to
minimize internal fragmentation of the array.

The `CountingBitset` is used in the `CountingBloomFilter` and `HyperLogLog` data
structures.

#### `Kmer`

The `Kmer` class represents a k-mer and provides methods for accessing the
integer representation of that k-mer. K-mers are represented as 64-bit unsigned
integers, with 2 bits per nucleotide.

The `Kmer` class also maintains the representation of the reverse complement of
the k-mer. The *canonical* representation of a k-mer is then defined as the
smaller of the two representations.

### IO Module

The IO module is responsible for reading and writing nucleotide sequence data
from FASTA files. It provides components for parsing input streams and
formatting output.

#### `FastaReader`

Example usage:
```cpp
FastaReader reader("input.fasta");
while (reader.next_sequence()) {
    std::string header = reader.get_header();
    char nucleotide;
    while (reader.next_nucleotide(nucleotide)) {
        // Process each nucleotide
    }
}
```

The `FastaReader` can handle multiple sequences in a single FASTA file.

#### `KmerWriter`
- Handles output of sequences with optional k-mer splicing
- The presence/absence information of each k-mer is indicated using uppercase (present) and lowercase (absent) letters

Example usage:
```cpp
KmerWriter writer("output.fasta", 5, true);  // k=5, splicing enabled
writer.write_header(">Sequence1");
for (char c : "ACGTACGT") {
    writer.add_nucleotide(c);
    writer.print_nucleotide(PRESENT);
}
writer.flush();
```

### Math Module

The Math module contains implementations of fast modular arithmetic operations.

The `Modulus` class provides efficient computation of modulo operations for a
fixed modulus value using the [Barrett reduction
algorithm](https://en.wikipedia.org/wiki/Barrett_reduction).

The `Modulus` class has two main methods:

- `uint64_t reduce(uint128_t)`

  Returns the result of `value % modulus` for a 128-bit integer.

- `uint64_t reduce2(uint128_t)`

  Returns a number `x` such that `x = value % modulus` and `0 <= x < 2 *
  modulus`. This method is slightly faster than `reduce` because it contains no
  branching.

### Sketch Module

The Sketch module contains implementations of `BloomFilter`,
`CountingBloomFilter` and `HyperLogLog` data structures.

All of these data structures have a template parameter for the hash function
family, which must satisfy the `HashFamily` concept. For the rolling variants
of `BloomFilter` and `CountingBloomFilter`, the hash family must satisfy the
`RollingHashFamily` concept.

---

## References

1. Cormode, Graham, H., & others. *Small Summaries for Big Data*, DIMACS Series in Discrete Mathematics and Theoretical Computer Science, Rutgers University / DIMACS. Available at: http://dimacs.rutgers.edu/~graham/ssbd/ssbd2.pdf

[small-summaries]: http://dimacs.rutgers.edu/~graham/ssbd/ssbd2.pdf

2. Austin Appleby. *MurmurHash3*. Available at: https://github.com/aappleby/smhasher

[murmurhash]: https://github.com/aappleby/smhasher

3. Kirsch, A., Mitzenmacher, M. (2006). Less Hashing, Same Performance: Building a Better Bloom Filter. In: Azar, Y., Erlebach, T. (eds) Algorithms – ESA 2006. ESA 2006. Lecture Notes in Computer Science, vol 4168. Springer, Berlin, Heidelberg. https://doi.org/10.1007/11841036_42

[kirsch-mitzmacker-2008]: https://doi.org/10.1007/11841036_42
