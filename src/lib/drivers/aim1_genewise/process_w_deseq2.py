# from https://github.com/wckdouglas/diffexpr/blob/master/example/deseq_example.ipynb

import pickle as pk

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.lib.classes import py_DESeq2

GCT_PATH = "<GCT PATH>"
design_matrix_path = "<DESIGN MATRIX>"
RES_BINARY_PATH = "<RES BINARY PATH>"
PLT_BINARY_PATH = "<PLT BINARY PATH>"
TEST_PLOT_PATH = "<TEST PLOT PATH>"

df = pd.read_table(GCT_PATH)
design_matrix = pd.read_table(design_matrix_path)

dds = py_DESeq2(count_matrix=df,
                design_matrix=design_matrix,
                design_formula='~ sample',
                gene_column='id')  # <- This is the DESeq2 "gene ID" column.
dds.run_deseq()
res = dds.normalized_result()

# save results
pk.dump(res, open(RES_BINARY_PATH, "rb"))

# try generating a test figure
plt.scatter(res.logF2FoldChange,
            -np.log2(res.padj))
plt.savefig(TEST_PLOT_PATH)
