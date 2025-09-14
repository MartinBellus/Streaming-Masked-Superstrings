#include "algorithm/exact.hpp"
#include "io/fasta.hpp"
#include <unordered_set>

using namespace exact;

std::ostream &operator<<(std::ostream &os, const Accuracy &acc) {
    double missing_percent =
            (double)acc.missing_kmers / acc.present_kmers * 100;
    double additional_percent =
            (double)acc.additional_kmers / acc.present_kmers * 100;
    os << "Missing kmers: " << acc.missing_kmers << " / " << acc.present_kmers
       << " (" << missing_percent << "%)\n"
       << "Additional kmers: " << acc.additional_kmers << " / "
       << acc.present_kmers << " (" << additional_percent << "%)\n";
    return os;
}

Accuracy exact::compute_accuracy(const CompareArgs &args) {
    auto K = args.k();
    io::FastaReader output(args.output());
    io::FastaReader golden_output(args.golden());
    auto kmer_repr =
            args.unidirectional() ? KmerRepr::FORWARD : KmerRepr::CANON;
    std::unordered_set<Kmer::data_t> kmer_set;
    Accuracy result;

    while (golden_output.next_sequence()) {
        Kmer kmer(K);
        char c;
        std::size_t mask = 0;
        while (golden_output.next_nucleotide(c)) {
            kmer.roll(c);
            mask = (mask << 1) | (bool)std::isupper(c);
            mask &= (1ULL << K) - 1;
            bool present = (mask & (1ULL << (K - 1))) != 0;
            if (kmer.available() >= K && present) {
                kmer_set.insert(kmer.data(kmer_repr));
            }
        }
    }

    result.present_kmers = kmer_set.size();

    while (output.next_sequence()) {
        Kmer kmer(K);
        char c;
        std::size_t mask = 0;
        while (output.next_nucleotide(c)) {
            kmer.roll(c);
            mask = (mask << 1) | (bool)std::isupper(c);
            mask &= (1ULL << K) - 1;
            bool present = (mask & (1ULL << (K - 1))) != 0;
            if (kmer.available() >= K && present) {
                auto kmer_data = kmer.data(kmer_repr);
                if (kmer_set.contains(kmer_data)) {
                    kmer_set.erase(kmer_data);
                } else {
                    result.additional_kmers++;
                }
            }
        }
    }

    result.missing_kmers = kmer_set.size();
    return result;
}

int exact::compute_superstring(const ExactArgs &args) {
    auto K = args.k();
    io::FastaReader in(args.dataset());
    io::KmerWriter out(args.output(), K, args.splice());
    auto kmer_repr =
            args.unidirectional() ? KmerRepr::FORWARD : KmerRepr::CANON;
    out.write_header(args.fasta_header());
    std::unordered_set<Kmer::data_t> kmer_set;

    while (in.next_sequence()) {
        Kmer kmer(K);
        char c;
        while (in.next_nucleotide(c)) {
            kmer.roll(c);
            out.add_nucleotide(c);
            if (kmer.available() < K) {
                continue;
            }
            auto kmer_data = kmer.data(kmer_repr);
            if (kmer_set.contains(kmer_data)) {
                out.print_nucleotide(io::NOT_PRESENT);
            } else {
                out.print_nucleotide(io::PRESENT);
                kmer_set.insert(kmer_data);
            }
        }
        out.flush();
    }
    return 0;
}
