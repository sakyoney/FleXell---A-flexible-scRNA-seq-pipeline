import sys
import pandas as pd

def aggregate_counts(file_list):
     dfs = [pd.read_csv(f, sep='\t', comment='#', index_col=0) for f in file_list]
     aggregated = pd.concat(dfs, axis=1)
     aggregated.columns = [f.split('_counts.txt')[0] for f in file_list]
     return aggregated

if __name__ == "__main__":
     file_list = sys.argv[1:]
     aggregated = aggregate_counts(file_list)
     aggregated.to_csv(sys.stdout, sep='\t')