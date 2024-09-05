import time
import json

import pandas as pd
import numpy as np
import sys
sys.path.append("private-pgm/mechanisms/")
sys.path.append("private-pgm/src/")
from mst import MST
from mbi import Domain, Dataset


DIR = "/Users/golobs/Documents/GradSchool/RHAIL/experiment_files/"

def main():
    files = [
        DIR + "JH_batch1_salmon.merged.gene_counts.tsv",
        DIR + "JHU_batch2_addl_salmon.merged.gene_counts.tsv",
        DIR + "JHU_batch2_salmon.merged.gene_counts.tsv",
        DIR + "JHU_batch3_salmon.merged.gene_counts.tsv",
        DIR + "JHU_mixed_batch_salmon.merged.gene_counts.tsv",
    ]
    full_data = pd.DataFrame()

    for f in files:
        df = pd.read_csv(f, sep='\t')
        df = df.iloc[:500].T
        df = df.rename(columns=df.iloc[0] + "_" + df.iloc[1])
        df = df.iloc[2:]
        full_data = pd.concat([full_data, df], ignore_index=False)

    full_data.to_csv(DIR + "SYNAPSE_gene_counts_full.csv")

def create_mst_data():
    start = time.process_time()
    df = pd.read_csv(DIR + "SYNAPSE_gene_counts_full.csv")
    cols = df.iloc[:, 1:20].columns
    # df = df.rename(columns=dict(zip(cols, list(map(str, range(19))))))
    # cols = list(map(str, range(19)))

    df[cols] = df[cols].astype(np.int64)
    ranges = dict(df[cols].max() - df[cols].min())
    # config = json.loads(json.dumps(ranges, cls=NpEncoder))
    # domain = Domain(config.keys(), config.values())

    domain = Domain(cols, ranges.values())
    data = Dataset(df[cols], domain)
    synth = MST(data, 1.0, 1e-9)


    end = time.process_time()
    print(round(end - start))



class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


if __name__ == '__main__':
    # main()
    create_mst_data()

