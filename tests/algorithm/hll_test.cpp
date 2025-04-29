#include "algorithm/approximate_count.hpp"
#include "hash/murmur_hash.hpp"
#include "hash/poly_hash.hpp"
#include <iostream>
#include <unordered_set>

class KmerKey {
  public:
    KmerKey(const Kmer &kmer, std::size_t hash) : kmer(kmer), _hash(hash) {}
    bool operator==(const KmerKey &other) const { return kmer == other.kmer; }
    std::size_t hash() const { return _hash; }

  private:
    Kmer kmer;
    std::size_t _hash;
};

struct KmerHash {
    std::size_t operator()(const KmerKey &kmer_key) const {
        return kmer_key.hash();
    }
};

std::size_t exact_count(io::FastaReader &in, std::size_t K) {
    std::unordered_set<KmerKey, KmerHash> kmer_set;
    Kmer kmer(K);
    poly_hash hash(K);
    char c;
    while (in.next_sequence()) {
        kmer.reset();
        for (std::size_t i = 0; i < K && in.next_nucleotide(c); i++) {
            kmer.roll(c);
        }
        hash.init(kmer);
        while (in.next_nucleotide(c)) {
            auto nout = kmer.last();
            kmer.roll(c);
            hash.roll(char_to_nucleotide(c), nout);
            KmerKey kmer_key(kmer, hash.get_hash());
            kmer_set.insert(kmer_key);
        }
    }
    return kmer_set.size();
}

template <typename H>
void hll_test(const std::string &path, std::size_t K) {
    io::FastaReader in(path);
    auto count = approximate_count<H>(in, K);
    in.reset();
    auto count_real = exact_count(in, K);

    std::cout << "Approximate count: \t" << count << "\n";
    std::cout << "Exact count: \t\t" << count_real << "\n";
    std::cout << "Error: "
              << static_cast<double>(count_real - count) /
                         static_cast<double>(count_real) * 100
              << "%\n";
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file>\n";
        return 1;
    }

    std::string path = argv[1];
    std::size_t K = 31;
    hll_test<murmur_hash_family>(path, K);
}
