import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import sys

# Activate pandas-to-R dataframe conversion
pandas2ri.activate()

# Load DESeq2 package
base = importr("base")
deseq2 = importr("DESeq2")

# Input arguments
counts_file = sys.argv[1]
output_file = sys.argv[2]

# Load count data
counts = pd.read_csv(counts_file, sep="\t", index_col=0)

# Define experimental design
# Update this to match your actual conditions
conditions = pd.Series(["control", "treated", "treated"], index=counts.columns)

# Convert data to R data types
counts_r = pandas2ri.py2rpy(counts)
conditions_r = robjects.FactorVector(conditions)

# Create DESeq2 dataset
dds = deseq2.DESeqDataSetFromMatrix(countData=counts_r,
                                    colData=robjects.DataFrame({"condition": conditions_r}),
                                    design=robjects.Formula("~ condition"))

# Run DESeq2 analysis
dds = deseq2.DESeq(dds)
res = deseq2.results(dds)

# Convert results to pandas DataFrame
res_df = pandas2ri.rpy2py(res)
res_df.index = counts.index

# Save results to CSV
res_df.to_csv(output_file)
