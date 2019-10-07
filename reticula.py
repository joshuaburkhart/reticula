## accept flags and run various parts of the analysis

# coding=utf-8
import sys

from src.lib.Alignment import MultipleAlignment
from src.lib.FileConverter import add_alignments_from_txt
from src.lib.ScoringFunction import entropy
from src.lib.ScoringFunction import sum_of_pairs
from src.lib.ScoringMatrix import blosum_62

__author__ = 'burkhart'
USAGE = 'Usage: python3 hw3.py <alignment file>\n' \
        'Example: python3 hw3.py ./data/input/hw3.txt'
# %% cell 0
if len(sys.argv) != 2:
    print(USAGE)
    exit()


#%% create cell

def print_entropy_and_sop(infile):
    """
    prints minimum entropy, entropy, and sum of pairs
    for a multiple alignment found in infile
    :param infile:
    """
    ma = MultipleAlignment()
    add_alignments_from_txt(infile, ma)

    print("Minimum Entropy: {0}".format(
        round(
            ma.minimize_aligned_char_scores(entropy),
            3)))

    print("Entropy: {0}".format(
        round(
            ma.sum_aligned_char_scores(entropy),
            3)))

    print("Sum of Pairs: {0}".format(
        round(
            ma.sum_aligned_char_scores(sum_of_pairs, blosum_62),
            3)))


print_entropy_and_sop(sys.argv[1])