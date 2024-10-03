import unittest

from main import find_palindromes, find_top_k_mers, find_invalid_bases


VALID_BASES = ["A", "B", "C"]
TOY_DATASET = {
    "num_sequences": 5,
    "sequence_length": 10,
    "sequences": [
        "AAAAAAAAAA",  # long palindromes but no diversity
        "ABBBBBBBBB",
        "CCCCCCCCCX",
        "ABCABCABCA",  # has no palindromes
        "ABACACBBCA",  # has odd palindromes at pos 0 and 2 and 3, even ones at pos 4 and 5
    ],
}


class TestPalindromes(unittest.TestCase):
    def setUp(self) -> None:
        self.min_length = 2
        self.min_diversity = 2

    def test_even_palindromes(self):
        idx_to_check = 4
        seq = TOY_DATASET["sequences"][idx_to_check]
        palindromes = find_palindromes(
            seq, min_substring_len=self.min_length, min_bases=self.min_diversity
        )
        even_palindromes = [
            p for p in palindromes["palindromes"] if p["length"] % 2 == 0
        ]
        self.assertEqual(len(even_palindromes), 2)

    def test_odd_palindromes(self):
        idx_to_check = 4
        seq = TOY_DATASET["sequences"][idx_to_check]
        palindromes = find_palindromes(
            seq, min_substring_len=self.min_length, min_bases=self.min_diversity
        )
        odd_palindromes = [
            p for p in palindromes["palindromes"] if p["length"] % 2 == 1
        ]
        self.assertEqual(len(odd_palindromes), 3)

    def test_no_palindromes(self):
        for idx_to_check in range(3):
            seq = TOY_DATASET["sequences"][idx_to_check]
            palindromes = find_palindromes(
                seq, min_substring_len=self.min_length, min_bases=self.min_diversity
            )
            self.assertEqual(palindromes["num_palindromes"], 0)

    def test_palindrome_min_len(self):
        # There is only 1 palindrome of length > 5 in this sequence:
        # ACBBCA (len 6 at pos 5)
        self.min_length = 5
        idx_to_check = 4
        seq = TOY_DATASET["sequences"][idx_to_check]
        palindromes = find_palindromes(
            seq, min_substring_len=self.min_length, min_bases=self.min_diversity
        )
        self.assertEqual(palindromes["num_palindromes"], 1)

    def test_palindrome_min_diversity(self):
        self.min_length = 2
        self.min_diversity = 3
        idx_to_check = 4
        seq = TOY_DATASET["sequences"][idx_to_check]
        palindromes = find_palindromes(
            seq, min_substring_len=self.min_length, min_bases=self.min_diversity
        )
        self.assertEqual(palindromes["num_palindromes"], 1)


class TestInvalid(unittest.TestCase):

    def test_all_valid(self):
        indices_to_check = [0, 1, 3, 4]
        for idx in indices_to_check:
            seq = TOY_DATASET["sequences"][idx]
            invalid = find_invalid_bases(seq, valid_bases=VALID_BASES)
            self.assertEqual(invalid["num_invalid"], 0)

    def test_invalid(self):
        idx = 2
        seq = TOY_DATASET["sequences"][idx]
        invalid = find_invalid_bases(seq, valid_bases=VALID_BASES)
        self.assertEqual(invalid["num_invalid"], 1)


class TestKMers(unittest.TestCase):
    def setUp(self) -> None:
        self.top_n = 3

    def test_5_mer_aggregate(self):
        # We expect the top 3 5-mers to be AAAAA, BBBBB and CCCCC
        aggregate_k_mers = find_top_k_mers(
            TOY_DATASET, top=self.top_n, k=5, per_sequence=False
        )
        ground_truth = ["AAAAA", "BBBBB", "CCCCC"]
        for i, kmer_tuple in enumerate(aggregate_k_mers):
            self.assertEqual(kmer_tuple[0][0], ground_truth[i])

    def test_5_mer_per_seq(self):
        per_seq_k_mers = find_top_k_mers(
            TOY_DATASET, top=self.top_n, k=5, per_sequence=True
        )

        # We expect the first sequence to only have 1 5-mer as
        # it is just made up of 1 nucleobase
        self.assertEqual(len(per_seq_k_mers[0]), 1)

        # The next 2 sequences must have 2 5-mers, because they
        # are made of the same base except for 1 character
        self.assertEqual(len(per_seq_k_mers[1]), 2)
        self.assertEqual(len(per_seq_k_mers[2]), 2)


if __name__ == "__main__":
    unittest.main()
