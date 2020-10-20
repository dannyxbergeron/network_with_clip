import numpy as np
import pandas as pd

from pybedtools import BedTool as bt

data_file = snakemake.input.sno_host
ref_file = snakemake.input.host_ref
begraph_file = snakemake.input.bedgraph


def load_data():
    df = pd.read_csv(data_file, sep='\t')
    return df


def load_ref():
    df = pd.read_csv(ref_file, sep='\t')
    return df


def load_bedgraph():
    df = pd.read_csv(bedgraph_file, sep='\t')
    return df


def get_intron(df_, ref_df):

    def get_loc(num):
        loc = df.at[i, f'loc{num}']
        id = df.at[i, f'ex_int_id{num}']

        if loc == 'intron' or loc == 'intron_exon':
            row = ref_df.loc[ref_df.ex_id == id].values[0]
            return row[1], row[2]
        elif loc == 'exon_intron':
            ex_int_num = df.at[i, f'ex_int_num{num}']
            transcript_id = ref_df.loc[ref_df.ex_id == id].values[0][6]
            row = ref_df.loc[(ref_df.ex_num == ex_int_num) &
                             (ref_df.trans_id == transcript_id) &
                             (ref_df.ex_id.str.contains('intron'))]
            return row[1], row[2]
        else:
            return np.nan, np.nan

    df = df.copy(deep=True)

    for i in df.index:
        intron_start1, intron_end1 = get_loc(1)
        intron_start2, intron_end2 = get_loc(2)


def bedtools(df1, df2):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=True, sorted=False)
    new_cols = ['chr1', 'start1', 'end1', 'gene_name1', 'DG', 'strand',
                'chr2', 'start2', 'end2', 'gene_id2', 'gene_name2', 'strand2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    return intersect_df


def main():

    df = load_data()
    ref_df = load_ref()
    bedgraph_df = load_bedgraph()

    df = get_intron(df, ref_df)


if __name__ == '__main__':
    main()
