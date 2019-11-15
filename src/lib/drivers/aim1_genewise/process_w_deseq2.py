# from https://github.com/wckdouglas/diffexpr/blob/master/example/deseq_example.ipynb

import argparse
import json as jp
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.lib.classes import py_DESeq2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", "-v", type=bool, default=False,
                        help="Provide verbose output during execution.")
    parser.add_argument("--datafile", "-d", type=str,
                        help="Provide the full path to a GCT (.gct) file.")
    parser.add_argument("--classfile", "-c", type=str,
                        help="Provide the full path to a DESeq2 class (.cls) file.")
    parser.add_argument("--output", "-o", type=str,
                        help="Provide an output file path.")
    return parser


# load arguments
args = build_parser().parse_args(sys.argv[1:])

# load data & class files
df = pd.read_table(Path(args.datafile))
design_matrix = pd.read_table(Path(args.classfile))

# execute DESeq2
dds = py_DESeq2(count_matrix=df,
                design_matrix=design_matrix,
                design_formula='~ sample',
                gene_column='Name')  # <- This is the DESeq2 "gene ID" column... should be "Name" in GCT
dds.run_deseq()
res = dds.normalized_result()  # TODO: confirm this is the preferred function to generate final DESeq2 output

# store results
with open(Path(args.output), 'w') as outfile:
    jp.dump(res, outfile)

# generate & store test figure
plt.scatter(res.logF2FoldChange,
            -np.log2(res.padj))
plt.savefig(Path(args.imgout))
