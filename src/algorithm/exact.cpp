#include "algorithm/exact.hpp"
#include "hash/poly_hash.hpp"
#include "helper/hashed_kmer.hpp"
#include "io/fasta.hpp"
#include <unordered_set>

using namespace exact;

std::ostream &operator<<(std::ostream &os, const Accuracy &acc) {
    double missing_percent = (double)acc.missing_kmers / acc.length * 100;
    double additional_percent = (double)acc.additional_kmers / acc.length * 100;
    os << "Missing kmers: " << acc.missing_kmers << " / " << acc.length << " ("
       << missing_percent << "%)\n"
       << "Additional kmers: " << acc.additional_kmers << " / " << acc.length
       << " (" << additional_percent << "%)\n";
    return os;
}

Accuracy exact::compute_accuracy(const CompareArgs &args) {
    io::FastaReader output(args.output());
    io::FastaReader golden_output(args.golden());
    Accuracy result;
    char out_c, golden_c;
    if (!output.next_sequence() || !golden_output.next_sequence()) {
        return result;
    }
    while (output.next_nucleotide(out_c) &&
           golden_output.next_nucleotide(golden_c)) {
        bool out = std::isupper(out_c);
        bool golden = std::isupper(golden_c);
        if (out && !golden) {
            result.additional_kmers++;
        } else if (!out && golden) {
            result.missing_kmers++;
        }
        result.length++;
    }
    return result;
}

int exact::compute_superstring(const ExactArgs &args) {
    using kmer_t = HashedKmer<poly_hash>;
    auto K = args.k();
    io::FastaReader in(args.dataset());
    io::KmerWriter out(args.output(), K, args.splice());
    auto kmer_repr =
            args.unidirectional() ? KmerRepr::FORWARD : KmerRepr::CANON;
    out.write_header(args.fasta_header());
    std::unordered_set<kmer_t, KmerHash> kmer_set;
    kmer_t kmer(K, 0, kmer_repr);
    char c;
    while (in.next_sequence()) {
        kmer.reset();
        for (std::size_t i = 0; i < K - 1 && in.next_nucleotide(c); i++) {
            kmer.roll(c);
            out.add_nucleotide(c);
        }
        while (in.next_nucleotide(c)) {
            kmer.roll(c);
            out.add_nucleotide(c);
            if (kmer_set.contains(kmer)) {
                out.print_nucleotide(io::NOT_PRESENT);
            } else {
                out.print_nucleotide(io::PRESENT);
                kmer_set.insert(kmer);
            }
        }
        out.flush();
    }
    return 0;
}
