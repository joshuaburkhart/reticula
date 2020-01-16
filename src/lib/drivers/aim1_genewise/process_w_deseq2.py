# from https://gist.githubusercontent.com/wckdouglas/3f8fb27a3d7a1eb24c598aa04f70fb25/raw/a9f019ae71f7ae46d576b2c4603a8123753b29eb/py_deseq.py
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, numpy2ri, Formula

numpy2ri.activate()
pandas2ri.activate()
from rpy2.robjects.packages import importr

deseq = importr('DESeq2')
'''
Adopted from: https://stackoverflow.com/questions/41821100/running-deseq2-through-rpy2
'''

to_dataframe = robjects.r('function(x) data.frame(x)')


class py_DESeq2:
    '''
    DESeq2 object through rpy2

    input:
    count_matrix: should be a pandas dataframe with each column as count, and a id column for gene id
        example:
        id    sampleA    sampleB
        geneA    5    1
        geneB    4    5
        geneC    1    2
    design_matrix: an design matrix in the form of pandas dataframe, see DESeq2 manual, samplenames as rownames
                treatment
    sampleA1        A
    sampleA2        A
    sampleB1        B
    sampleB2        B

    design_formula: see DESeq2 manual, example: "~ treatment""
    gene_column: column name of gene id columns, exmplae "id"
    '''

    def __init__(self, count_matrix, design_matrix, design_formula, gene_column='id'):
        try:
            assert gene_column in count_matrix.columns, 'Wrong gene id column name'
            gene_id = count_matrix[gene_column]
        except AttributeError:
            sys.exit('Wrong Pandas dataframe?')

        self.dds = None
        self.deseq_result = None
        self.resLFC = None
        self.comparison = None
        self.normalized_count_df = None
        self.gene_column = gene_column
        self.gene_id = count_matrix[self.gene_column]
        self.samplenames = count_matrix.columns[count_matrix.columns != self.gene_column]
        self.count_matrix = pandas2ri.py2rpy_pandasdataframe(count_matrix.set_index(self.gene_column))
        self.design_matrix = pandas2ri.py2rpy_pandasdataframe(design_matrix)
        self.design_formula = Formula(design_formula)
        self.dds = deseq.DESeqDataSetFromMatrix(countData=self.count_matrix,
                                                colData=self.design_matrix,
                                                design=self.design_formula)

    def run_deseq(self, **kwargs):
        self.dds = deseq.DESeq(self.dds, **kwargs)

    def get_deseq_result(self, contrast=None, **kwargs):

        self.comparison = deseq.resultsNames(self.dds)
        if contrast:
            if len(contrast) == 3:
                contrast = numpy2ri.numpy2rpy(np.array(contrast))
            else:
                assert len(contrast) == 2, 'Contrast must be length of 3 or 2'
                contrast = robjects.ListVector({None: con for con in contrast})
            print('Using contrast: ', contrast)
            self.deseq_result = deseq.results(self.dds, contrast=contrast, **kwargs)
        else:
            self.deseq_result = deseq.results(self.dds, **kwargs)
        self.deseq_result = to_dataframe(self.deseq_result)
        self.deseq_result = pandas2ri.rpy2py_dataframe(self.deseq_result)  ## back to pandas dataframe
        self.deseq_result[self.gene_column] = self.gene_id.values

    def normalized_count(self):
        normalized_count_matrix = deseq.counts_DESeqDataSet(self.dds, normalized=True)
        normalized_count_matrix = to_dataframe(normalized_count_matrix)
        # switch back to python
        self.normalized_count_df = pandas2ri.rpy2py_dataframe(normalized_count_matrix)
        self.normalized_count_df[self.gene_column] = self.gene_id.values
        return self.normalized_count_df


# from https://github.com/wckdouglas/diffexpr/blob/master/example/deseq_example.ipynb

import argparse
import json as jp
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Process RNA sequencing data using DESeq2.')
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
args = build_parser().parse_args()

print(Path(args.datafile))
print(Path(args.classfile))

# load data & class files
df = pd.read_csv(Path(args.datafile))
print(df.head())
print(df.size)
print(df.shape)
print(df.ndim)

ds = pd.read_csv(Path(args.classfile))
ds.index = ds.SAMPID
print(ds.head())
print(ds.size)
print(ds.shape)
print(ds.ndim)

# execute DESeq2
dds = py_DESeq2(count_matrix=df,
                design_matrix=ds,
                design_formula="~ SMTS",
                gene_column='id')  # <- This is the DESeq2 "gene ID" column... should be "id" in GCT
dds.run_deseq()
dds.get_deseq_result()
res = dds.deseq_result
res.head()
res = dds.normalized_count()  # TODO: confirm this is the preferred function to generate final DESeq2 output

# store results
with open(Path(args.output), 'w') as outfile:
    jp.dump(res, outfile)

# generate & store test figure
plt.scatter(res.logF2FoldChange,
            -np.log2(res.padj))
plt.savefig(Path(args.imgout))
