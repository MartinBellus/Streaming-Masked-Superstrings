# Algorithm overview

## Streaming algorithm

The main data structure used in the streaming algorithm is a Bloom Filter[^1]. A
Bloom Filter is a space-efficient representation of a set, supporting insertions
and membership queries. The queries have a guarantee of no false negatives, but
can have false positives.

[^1]: See chapter 2, section 2.7 of [Small Summaries for Big Data](http://dimacs.rutgers.edu/~graham/ssbd/ssbd2.pdf) for more details.

### First phase (Computation)

An example pseudocode for the first phase of the algorithm:
```py
filter = BloomFilter(optimal_size)

for i in range(len(sequence) - K):
    kmer = sequence[i:i+K]
    if not filter.contains(kmer):
        filter.insert(kmer)
        # Mark kmer as present (uppercase)
    else:
        # Mark kmer as not present (lowercase)
```

Because of the probabilistic nature of Bloom Filters, some k-mers might not be
represented in the final output.

### Second phase (Correction)

The input for this phase is a sequence produced by the first phase, where some
of the k-mers are incorrectly marked as not present.

An example pseudocode for the second phase of the algorithm:
```py
filter = CountingBloomFilter(optimal_size)

# First pass: Insert all possibly missing k-mers
for i in range(len(sequence) - K):
    kmer = sequence[i:i+K]
    if sequence[i].islower():  # k-mer marked as not present
        filter.insert(kmer)

# Second pass: Remove all present k-mers
for i in range(len(sequence) - K):
    kmer = sequence[i:i+K]
    if sequence[i].isupper():  # k-mer marked as present
        filter.remove(kmer)

# Third pass: Correct the output
for i in range(len(sequence) - K):
    kmer = sequence[i:i+K]
    if sequence[i].isupper() or filter.contains(kmer):
        filter.remove(kmer)
        # Mark kmer as present (uppercase)
    else:
        # Mark kmer as not present (lowercase)
```

The data structure used in this phase is a Counting Bloom Filter, allowing for
both insertions and deletions of k-mers. The filter represents a set, therefore
we only insert k-mers that are not already represented in the filter.

After the first pass, the Counting Bloom Filter contains missing and repeating
k-mers (but possibly not all of them).

In the second pass, we delete all of the k-mers marked as present in the first phase.

After that, the filter should contain mostly missing k-mers. Because of the
false positives, it can contain multiple copies of the same k-mer. Therefore the
output of the second phase can contain k-mers that are represented more than
once.

### Estimating k-mer count

To set the size of the Bloom Filters, we need to estimate the number of elements
$N$ that will be inserted into them. Given the fixed parameter of bits per
k-mer $b$, the size of the Bloom filter will be $M = b \cdot N$ bits.

In the first phase, the number of inserted elements is exactly the number of
unique k-mers in the input sequence. This can be estimated using the
HyperLogLog[^2] algorithm, which provides a good balance between accuracy and
memory efficiency.

In the second phase, the number of k-mers inserted into the Counting Bloom
Filter is proportional to the number of distinct k-mers marked as not present
during the first phase. We upper bound this by the number of repeated k-mers in
the input sequence.

Given the total length of the input sequences $L$, the number of sequences $m$
and the number of unique k-mers $N$, the number of repeated k-mers is $L - m
\cdot (K - 1) - N$.

[^2]: See chapter 2, section 2.6 of [Small Summaries for Big Data](http://dimacs.rutgers.edu/~graham/ssbd/ssbd2.pdf) for more details.

---

## References

1. Cormode, Graham, H., & others. *Small Summaries for Big Data*, DIMACS Series in Discrete Mathematics and Theoretical Computer Science, Rutgers University / DIMACS. Available at: http://dimacs.rutgers.edu/~graham/ssbd.html
