from statistics import mean

import pandas as pd
import numpy as np

from pybedtools import BedTool as bt

data_file = snakemake.input.alt_splice
bg_files = snakemake.input.bg

bed_viz = snakemake.output.bed_viz
out_file = snakemake.output.ext_ratio


def load_df(file, header=True):
    if header:
        df = pd.read_csv(file, sep='\t')
    else:
        df = pd.read_csv(file, sep='\t',
                         names=['chr', 'start', 'end', 'value'])
    return df

def load_beds(files):

    dfs = {}
    for file in files:
        name = file.split('/')[-1].split('.')[0].replace('intron_', '')
        dfs[name] = load_df(file, header=False)

    return dfs


def get_regions(df_):

    THRESH = 5
    OFFSET = 2

    df = df_.copy(deep=True)
    df = df.loc[~pd.isnull(df.int_portion_start2)]
    df_val = df[[
        'chr1', 'int_portion_start2', 'intron_start', 'intron_end',
        'sno_start', 'sno_end', 'strand1', 'DG', 'name1', 'name2'
    ]].values

    master_list = []
    for chr, portion2, int_start, int_end, sno_start, sno_end, strand, DG, name1, name2 in df_val:
        p_start, p_end = [int(x) for x in portion2.split('-')]
        sno_start -= (OFFSET + 1)
        sno_end += OFFSET
        if p_start < sno_start:
            ext_start, ext_end = p_start, sno_start
            snoExt_start, snoExt_end = p_start, sno_end
        else:
            ext_start, ext_end = sno_end, p_end
            snoExt_start, snoExt_end = sno_start, p_end

        ext = [chr, ext_start, ext_end, DG, 'ext', strand]
        master_list.append(ext)

        left = [chr, int_start, snoExt_start, DG, 'left', strand]
        right = [chr, snoExt_end, int_end - 1, DG, 'right', strand]

        if snoExt_start - int_start > THRESH:
            master_list.append(left)
        if int_end - snoExt_end > THRESH:
            master_list.append(right)

    master_df = pd.DataFrame(master_list, columns=['chr', 'start', 'end', 'DG', 'side', 'strand'])
    master_df['chr'] = 'chr' + master_df['chr']
    master_df['start'] = master_df.start.map(int)
    master_df['end'] = master_df.end.map(int)
    master_df.to_csv(bed_viz, sep='\t', index=False, header=False)

    return master_df

def bedtools(df1, df2, sum_bp):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=False, sorted=False)
    new_cols = ['chr', 'start', 'end', 'DG', 'side', 'strand',
                'chr2', 'start2', 'end2', 'score', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})

    intersect_df['sum_score'] = intersect_df['score'] * intersect_df['overlap']
    average = sum(intersect_df['sum_score']) / sum_bp

    return average


def get_values(df_, beds):

    df = df_.copy(deep=True)
    # Just for testing with SNORD2...
    # df = df.loc[df.DG == '273061_L0|273063_L0|488949_L1|488953_L1|692461_P0']

    means_ratio = {}
    for DG in set(df.DG):
        tmp = df.loc[df.DG == DG]

        ext = tmp.loc[df.side == 'ext']
        rest = tmp.loc[~(df.side == 'ext')]

        ext_sum = sum([
            ext.at[x, 'end'] - ext.at[x, 'start']
            for x in ext.index
        ])
        rest_sum = sum([
            rest.at[x, 'end'] - rest.at[x, 'start']
            for x in rest.index
        ])

        ratios = []
        for name, bed in beds.items():

            ext_val = bedtools(ext, bed, ext_sum)
            rest_val = bedtools(rest, bed, rest_sum)

            ratio = ext_val / rest_val
            ratios.append(ratio)

        mean_ratio = mean(ratios)
        means_ratio[DG] = mean_ratio
        print(DG)
        print(mean(ratios))

    return means_ratio

def main():

    df = load_df(data_file)
    beds = load_beds(bg_files)

    df_regions = get_regions(df)

    ratio_dict = get_values(df_regions, beds)

    df['ext_ratio'] = df['DG'].map(ratio_dict)

    print(df)
    df.to_csv(out_file, sep='\t', index=False)



if __name__ == '__main__':
    main()
