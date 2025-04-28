#include "algorithm/exact.hpp"
#include "hash/poly_hash.hpp"
#include "helper/kmer.hpp"
#include <unordered_set>

using namespace exact;

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

int exact::compute_superstring(std::size_t K, io::FastaReader &in,
                               io::KmerWriter &out) {
    std::unordered_set<KmerKey, KmerHash> kmer_set;
    Kmer kmer(K);
    poly_hash hash(K);
    char c;
    while (in.next_sequence()) {
        kmer.reset();
        for (std::size_t i = 0; i < K && in.next_nucleotide(c); i++) {
            kmer.roll(c);
            out.add_nucleotide(c);
        }
        hash.init(kmer);
        while (in.next_nucleotide(c)) {
            auto nout = kmer.last();
            kmer.roll(c);
            hash.roll(char_to_nucleotide(c), nout);
            KmerKey kmer_key(kmer, hash.get_hash());
            if (kmer_set.contains(kmer_key)) {
                out.print_nucleotide(io::NOT_PRESENT);
            } else {
                out.print_nucleotide(io::PRESENT);
                kmer_set.insert(kmer_key);
            }
            out.add_nucleotide(c);
        }
        out.flush();
    }
    return 0;
}
