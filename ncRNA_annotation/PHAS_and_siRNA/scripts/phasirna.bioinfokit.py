import sys
import pandas as pd
import numpy as np
import bioinfokit
from pandas import read_csv
from bioinfokit.analys import norm


# load sRNA expression matrix
path=snakemake.input[0]
df = pd.read_table(path)
df.head(5)


# drop some columns
#df = df.drop(['Locus', 'main'], axis=1)
#df.head(5)

# make the column of locus name as index column
df = df.set_index('Name')
df.head(5)

# now, normalize raw counts using CPM method
nm = norm()
nm.cpm(df=df)

# get CPM normalized dataframe
cpm_df = nm.cpm_norm
cpm_df.head(2)

cpm_df.to_csv(snakemake.output[0], sep="\t")
