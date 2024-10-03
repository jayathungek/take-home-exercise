import sys
import argparse
from collections import Counter
from typing import List, Dict, Callable, Generator, Union, Tuple

from utils import load_json, save_json, plot_dna_dataset, plot_dnt_confmat


def substrings(seq: str, win_len: int, hop_len: int) -> Generator:
    """
    Moves a sliding window over a DNA sequence to generate
    substrings of a given length at a given interval.

    Args:
        seq (str): The sequence to move the window over
        win_len (int): The length of the K-mer
        hop_len (int): Distance to move the window for next K-mer

    Yields:
        Generator: The next K-mer in the sequence
    """
    seq_len = len(seq)
    for i in range(0, seq_len, hop_len):
        substring = seq[i : i + win_len]
        if len(substring) != win_len:
            break  # We have overshot the length of the sequence
        yield substring


def get_gc_stats(sequence: str, valid_nucleobases: Union[List, None] = None) -> Dict:
    """
    Returns some GC-related stats for a single DNA sequence. Calculates:
        - The GC distribution
        - GC skew

    Args:
        sequence (str): The DNA sequence to analyse
        valid_nucleobases (List, optional): List of characters
        representing the valid nucleobases. Any other bases are ignored
        for the purposes of statistical calculations. If None, all bases found in the
        sequence are treated as valid. Defaults to None.

    Returns:
        Dict: Object with all stats for the sequence
    """
    if valid_nucleobases is None:
        total_valid = len(sequence)
    else:
        num_invalid = len([base for base in sequence if base not in valid_nucleobases])
        total_valid = len(sequence) - num_invalid

    num_c = sequence.count("C")
    num_g = sequence.count("G")
    gc_distribution = (num_c + num_g) / total_valid
    gc_skew = (num_g - num_c) / (num_g + num_c)

    gc_stats = {"gc_distribution": gc_distribution, "gc_skew": gc_skew}

    return gc_stats


def get_dinucleotide_freqs(sequence: str) -> Dict:
    """
    Finds the frequency of all 2-mers in a single DNA sequence
    Args:
        sequence (str): The DNA sequence to analyse

    Returns:
        Dict: Object with all stats for the sequence
    """
    valid_dinucleotides = [
        i + j
        for i in SETTINGS["valid_nucleobases"]
        for j in SETTINGS["valid_nucleobases"]
    ]
    dinucleotide_counts = {d: 0 for d in valid_dinucleotides}
    dinucleotides = substrings(sequence, win_len=2, hop_len=2)
    for d in dinucleotides:
        if d in dinucleotide_counts.keys():
            dinucleotide_counts[d] += 1
    total_dinucleotides = sum(dinucleotide_counts.values())
    dinucleotide_frac = {
        d: count / total_dinucleotides for d, count in dinucleotide_counts.items()
    }
    return dinucleotide_frac


def find_invalid_bases(sequence: str) -> Dict:
    """
    Finds the positions of all invalid bases in a single DNA sequences,
    as defined in the settings.
    Args:
        sequence (str): The DNA sequence to analyse

    Returns:
        Dict: Object with all stats for the sequence
    """
    invalid_dict: Dict = {"num_invalid": 0, "invalid_bases": []}
    for pos, base in enumerate(sequence):
        if base not in SETTINGS["valid_nucleobases"]:
            invalid_dict["num_invalid"] += 1
            invalid_dict["invalid_bases"].append({"pos": pos, "base": base})
    return invalid_dict


def find_palindromes(sequence: str, min_substring_len: int, min_bases: int) -> Dict:
    """
    Finds palindromes of a minimum length and nucleobase count
    in a single DNA sequence.

    Args:
        sequence (str): The DNA sequence to analyse
        min_substring_len (int): Minimum length of the palindromic sequence
        min_bases (int): Minimum number of nucleobases that must be present 
        in the palindrome

    Returns:
        Dict: Object with all stats for the sequence. The start
        position in the DNA sequence and the length are reported for
        each palindrome found.
    """
    palindrome_stats: Dict = {"num_palindromes": 0, "palindromes": []}

    def expand_frontier(start: int, end: int):
        """
        Iterate through the sequence, treating each character (or pair
        of characters) as the centre of a possible palindrome. Keep 
        expanding outwards from the centre until a boundary is reached
        on either side of the sequence. Only add to the list of palindromes
        if long enough and contains enough unique bases.
        """
        while start >= 0 and end < len(sequence) and sequence[start] == sequence[end]:
            substring = sequence[start : end + 1]
            if (
                len(substring) > min_substring_len
                and len(Counter(substring)) >= min_bases
            ):
                palindrome_stats["num_palindromes"] += 1
                palindrome_stats["palindromes"].append(
                    {"pos": start, "length": len(substring)}
                )
            start -= 1
            end += 1

    for i in range(len(sequence)):
        expand_frontier(i, i)  # Find odd-numbered palindromes
        expand_frontier(i, i + 1)  # Find even-numbered palindromes

    return palindrome_stats


def find_k_mers(sequence: str, k: int) -> Dict:
    """_summary_
    Finds all substrings of length k (k-mers) of a single DNA
    sequence

    Args:
        sequence (str): The DNA sequence to analyse
        k (int): Length of the substring

    Returns:
        Dict: A mapping of each K-mer to its count
    """
    k_mer_stats: Dict = {}

    for substring in substrings(sequence, win_len=k, hop_len=1):
        if substring not in k_mer_stats.keys():
            k_mer_stats[substring] = 0
        k_mer_stats[substring] += 1

    return k_mer_stats


def find_top_k_mers(
    seq_object: Dict, top: int, k: int, per_sequence: bool = False
) -> List:
    """
    Reports the top N K-mers for all the DNA sequences in a
    given sequence object

    Args:
        seq_obj (Dict): Object containing a list of DNA sequences
        top (int): How many top K-mers to report
        per_sequence (bool, optional): Whether or not to aggregate the results
        of the top K-mers. If True, reports the top N K-mers for each DNA
        sequence individually, otherwise reports the top N for
        the whole dataset. Defaults to False.

    Returns:
        List: If arg `per_sequence` is False, the list contains only 1 element: 
        a list of tuples showing the aggregated results. If True, 
        it contains a list of tuples for each DNA sequence
    """

    def get_top_values(d: Dict) -> List:
        items = list(d.items())
        # Sort by count of dict values
        items.sort(key=lambda kv: kv[1], reverse=True)
        top_n = items[:top]
        return top_n

    all_k_mer_stats = apply_foreach(seq_object, find_k_mers, k=k)
    filtered = []
    if per_sequence:
        for seq_results in all_k_mer_stats:
            k_mer_stats = get_top_values(seq_results)
            filtered.append(k_mer_stats)
    else:
        aggregated_k_mer_stats: Dict = {}
        for seq_results in all_k_mer_stats:
            for key, value in seq_results.items():
                if key not in aggregated_k_mer_stats.keys():
                    aggregated_k_mer_stats[key] = 0
                aggregated_k_mer_stats[key] += value
        k_mer_stats = get_top_values(aggregated_k_mer_stats)
        filtered.append(k_mer_stats)
    return filtered


def apply_foreach(seq_object: Dict, analysis_func: Callable, **kwargs) -> List:
    """
    Utility function that applies a given analysis function to a
    collection of DNA sequences. Can pass in keyword arguments
    as required by the analysis function.

    Args:
        seq_obj (Dict): Object containing a list of DNA sequences
        analysis_func (Callable): A function that applies some analysis to
        a single DNA sequence

    Returns:
        List: A collection of DNA stats relevant to whatever analysis function was
        passed to `analysis_func`
    """
    stats: List = []
    for seq in seq_object["sequences"]:
        seq_stats = analysis_func(seq, **kwargs)
        stats.append(seq_stats)
    return stats


def find_max(stats: List, field: str) -> Tuple[int, Dict]:
    """
    Finds the entry with the maximum value for a given field.
    Throws KeyError if field is not found in the stat dict

    Args:
        stats (List): Collection of stat objects
        field (str): Field of interest within the object

    Returns:
        Tuple[int, Dict]: The position where the max value was found, and the
        maximum value object itself
    """
    max_index = 0
    max_value = float("-inf")
    ret_dict = stats[0]
    for i, stats_dict in enumerate(stats):
        if stats_dict[field] > max_value:
            max_value = stats_dict[field]
            ret_dict = stats[i]
            max_index = i
    return max_index, ret_dict


def find_min(stats: List, field: str) -> Tuple[int, Dict]:
    """
    Finds the entry with the minimum value for a given field.
    Throws KeyError if field is not found in the stat dict

    Args:
        stats (List): Collection of stat objects
        field (str): Field of interest within the object

    Returns:
        Tuple[int, Dict]: The position where the min value was found, 
        and the minimum value object itself
    """
    min_index = 0
    min_value = float("inf")
    ret_dict = stats[0]
    for i, stats_dict in enumerate(stats):
        if stats_dict[field] < min_value:
            min_value = stats_dict[field]
            ret_dict = stats[i]
            min_index = i
    return min_index, ret_dict


def reduce_avg(stats: List) -> Dict:
    """
    Averages all the stats for a given collection of results.

    Args:
        stats (List): Output of the `apply_foreach` function. A list of
        stats pertaining to a collection of DNA sequences

    Returns:
        Dict: The average values of all the DNA sequences for each
        field in the analysis
    """
    aggregated: Dict = {key: 0.0 for key, val in stats[0].items()}
    for stats_dict in stats:
        for key, value in stats_dict.items():
            aggregated[key] += value

    aggregated = {key: val / len(stats) for key, val in aggregated.items()}
    return aggregated


def parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="main")
    parser.add_argument("seq_file", help="Path to JSON file containing DNA sequences")
    parser.add_argument("-s", "--settings-file", default="./settings.json", type=str)
    parser.add_argument("-p", "--per_seq_k_mer", action="store_true")
    parser.add_argument("-o", "--out-dir", default="./results", type=str)
    parsed = parser.parse_args(args)
    return parsed


if __name__ == "__main__":
    run_args = parse_args(sys.argv[1:])
    SETTINGS = load_json(run_args.settings_file)
    seq_object = load_json(run_args.seq_file)

    # Palindromic substrings over 20 bases long that contain at least 3 bases
    print("Finding palindromes...", end='', flush=True)
    all_palindrome_stats = apply_foreach(
        seq_object,
        find_palindromes,
        min_substring_len=SETTINGS["min_basepair_len"],
        min_bases=SETTINGS["min_bases"],
    )
    seq_dict: Dict = {}
    for i, stats in enumerate(all_palindrome_stats):
        if stats["num_palindromes"] > 0:
            palindromes = [
                (
                    p["pos"],
                    seq_object["sequences"][i][p["pos"] : p["pos"] + p["length"]],
                )
                for p in stats["palindromes"]
            ]
            palindrome_dict = {pos: palindrome for pos, palindrome in palindromes}
            seq_dict[i] = palindrome_dict
    save_json(f"{run_args.out_dir}/palindrome_stats.json", seq_dict)
    print("done")

    # Find k-mers both per-sequence and overall
    print(
        f"Finding top K-mer stats {'' if run_args.per_seq_k_mer else '(aggregate)'} ...",
        end='', 
        flush=True
    )
    k_mer_dict: Dict = {}
    for k in SETTINGS["k_values"]:
        per_seq_k_mers = find_top_k_mers(
            seq_object, top=SETTINGS["top_n"], k=k, per_sequence=run_args.per_seq_k_mer
        )
        k_mer_dict[f"{k}_mers"] = per_seq_k_mers
    save_json(
        f"{run_args.out_dir}/k_mer_stats{'' if run_args.per_seq_k_mer else '_aggregate'}.json",
        k_mer_dict,
    )
    print("done")

    # GC content stats
    print("Finding GC stats...", end='', flush=True)
    all_gc_stats = apply_foreach(seq_object, get_gc_stats)
    gc_dict = {
        "avg_gc": reduce_avg(all_gc_stats),
        "max_gc_dist": find_max(all_gc_stats, "gc_distribution"),
        "max_gc_skew": find_max(all_gc_stats, "gc_skew"),
        "min_gc_dist": find_min(all_gc_stats, "gc_distribution"),
        "min_gc_skew": find_min(all_gc_stats, "gc_skew"),
    }
    save_json(f"{run_args.out_dir}/gc_stats.json", gc_dict)
    print("done")

    # Dinucleotide frequencies
    print("Calculating dinucleotide stats...", end='', flush=True)
    all_dinucleotide_stats = apply_foreach(seq_object, get_dinucleotide_freqs)
    dnt_dict = {
        "all_dnt_stats": all_dinucleotide_stats,
        "avg_dnt": reduce_avg(all_dinucleotide_stats),
    }
    save_json(f"{run_args.out_dir}/dnt_stats.json", dnt_dict)
    print("done")

    # Sequences with errors in them, maybe? Assuming only A, C, G, T are valid bases
    print("Searching for invalid sequences...", end='', flush=True)
    all_invalid_stats = apply_foreach(seq_object, find_invalid_bases)
    invalid_dict: Dict = {}
    for i, stat in enumerate(all_invalid_stats):
        if stat["num_invalid"] > 0:
            invalid_bases = [
                (s["pos"], seq_object["sequences"][i][s["pos"]])
                for s in stat["invalid_bases"]
            ]
            invalid_dict[i] = {pos: invalid_base for pos, invalid_base in invalid_bases}
    save_json(f"{run_args.out_dir}/invalid_stats.json", invalid_dict)
    print("done")

    # Generate graphs
    avg_dnt = dnt_dict["avg_dnt"]
    print("Creating dinucleotide frequency matrix...", end='', flush=True)
    plot_dnt_confmat(
        avg_dnt,
        bases=SETTINGS["valid_nucleobases"],
        outfile=f"{run_args.out_dir}/dnt_confmat.png",
    )
    print("done")

    print("Creating DNA visualisation...", end='', flush=True)
    plot_dna_dataset(
        seq_object,
        bases=SETTINGS["valid_nucleobases"] + ["X"],
        outfile=f"{run_args.out_dir}/dna_visualisation.png",
    )
    print("done")
