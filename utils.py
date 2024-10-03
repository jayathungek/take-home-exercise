import json
import colorsys
from typing import List, Dict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm


def get_colours(n: int) -> List:
    """
    Returns a list of RGB values that are evenly spaced in the
    HLS spectrum, such that they are easily distinguishable
    Args:
        n (int): number of colours needed

    Returns:
        List: List of normalised RGB values
    """
    colours = []
    hue_step = 360.0 / n

    for i in range(n):
        hue = i * hue_step
        saturation = 1
        brightness = 0.4

        rgb = colorsys.hls_to_rgb(hue / 360.0, brightness, saturation)
        colours.append([v for v in rgb])

    return colours


def seq_to_img(sequences: List, colours: List, bases: List) -> np.ndarray:
    """
    Transforms a List of strings representing DNA sequences to a
    N x L x 3 array, where N is the number of DNA sequences,
    L is the length of each sequence. Each nucleobase is converted
    to an RGB pixel so that it can be drawn.

    Args:
        sequences (List): List of DNA sequences
        colours (List): Colours to use
        bases (List): Valid nucleobases

    Returns:
        np.ndarray: The NxLx3 array which can be plotted
    """
    colour_lookup = {base: colour for base, colour in zip(bases, colours)}

    rows = []
    for seq in sequences:
        row = []
        for base in seq:
            try:
                row.append(colour_lookup[base])
            except KeyError:
                # invalid base found, use reserved error key
                row.append(colour_lookup["X"])
        rows.append(row)

    return np.array(rows)


def plot_dna_dataset(seq_obj: Dict, bases: List, outfile: str) -> None:
    """
    Represents a list of DNA sequences as an image, with each nucleobase
    given a unique, easily identifiable colour so that large patterns can be found

    Args:
        seq_obj (Dict): Object containing a list of DNA sequences
        bases (List): List of valid bases
        outfile (str): Where to save the image
    """
    sequences = seq_obj["sequences"]
    colours = get_colours(len(bases))
    custom_cmap = ListedColormap(colours)
    dna_image = seq_to_img(sequences, colours, bases)
    norm = BoundaryNorm(np.arange(len(bases) + 1) - 0.5, len(bases))
    plt.imshow(dna_image, cmap=custom_cmap, norm=norm)
    cbar = plt.colorbar(ticks=np.arange(len(bases)), label="Nucleobases")
    cbar.ax.set_yticklabels(bases)
    plt.savefig(outfile, dpi=300, bbox_inches="tight", format="png")
    plt.cla()


def plot_dnt_confmat(dnt_stats: Dict, bases: List, outfile: str) -> None:
    """
    Represents the given dinucleotide frequency stats as a
    confusion matrix
    Args:
        seq_obj (Dict): Stats for each dinucleotide pairing
        bases (List): List of valid bases
        outfile (str): Where to save the image
    """
    num_bases = len(bases)
    values = np.array(list(dnt_stats.values()))
    confmat = values.reshape(num_bases, num_bases)
    fig, ax = plt.subplots()
    for i in range(confmat.shape[0]):
        for j in range(confmat.shape[1]):
            ax.text(
                j, i, f"{confmat[i, j]:.4f}", ha="center", va="center", color="white"
            )

    ax.set_xticks(np.arange(num_bases))
    ax.set_yticks(np.arange(num_bases))
    ax.set_xticklabels(bases)
    ax.set_yticklabels(bases)
    ax.imshow(confmat)
    plt.savefig(outfile, dpi=300, bbox_inches="tight", format="png")
    plt.cla()


def load_json(path: str) -> Dict:
    """
    Convenience function for loading JSON
    Args:
        path (str): Path to JSON file
    Returns:
        Dict: The loaded JSON object
    """
    with open(path) as fh:
        loaded = json.load(fh)
    return loaded


def save_json(path: str, obj: object) -> None:
    """
    Convenience function for saving JSON

    Args:
        path (str): JSON file's save path
        obj (Dict): Object to serialise
    """
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=4)
 