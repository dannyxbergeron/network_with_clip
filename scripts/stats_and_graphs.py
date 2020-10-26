from collections import Counter

import numpy as np
import pandas as pd

in_file = snakemake.input.data_file
cons_file = snakemake.input.cons

def load_df():
    df = pd.read_csv(in_file, sep='\t')
    df.drop_duplicates(subset=['DG'], inplace=True)

    good_biotypes = [
        'snoRNA', 'protein_coding', 'lncRNA', 'rRNA', 'tRNA',
        'snRNA', 'scaRNA', 'miRNA',
    ]
    filtered_biotype = []
    for bio in df['biotype2'].values:
        if bio in good_biotypes:
            filtered_biotype.append(bio)
        else:
            filtered_biotype.append('other')

    df['filtered_biotype2'] = filtered_biotype
    return df


def stats(df, cons_df):

    # ----------------- Interacting biotypes --------------------
    print('============= all interactions count ================')
    print('Initial nb of interactions after merging and filtering: {}'.format(len(df)), end='\n\n')
    counter = Counter(list(df.filtered_biotype2))
    print('-- Biotypes interacting with snoRNA after removing intergenic and filtering --')
    print(counter.items(), end='\n\n')
    # ----------------- Interacting biotypes --------------------

    # ----------------- gene_host_interactions --------------------
    print('============= host interactions count ================')
    df = df.loc[df.host_interaction]
    print('Nb of snoRNA-host interactions: {}'.format(len(df)), end='\n\n')
    counter = Counter(list(df.filtered_biotype2))
    print(counter.items(), end='\n\n')

    # 10 were lost because of their nonsens thing (either not the good
    # overlapping gene or just no protein_coding transcript...)
    print('Number of filtered sno-protein_coding interactions: {}'.format(len(cons_df)))
    # print(set(df.loc[df.filtered_biotype2 == 'protein_coding'].DG)- set(cons_df.DG))

    intra_df = cons_df.loc[cons_df.interaction_type == 'intra']
    print('Number of intra interactions: {}'.format(len(intra_df)))

    sno_host = intra_df.drop_duplicates(subset=['merged_name'])
    print('Number of sno_host involved: {}'.format(len(sno_host)), end='\n\n')
    # print(list(sno_host.single_id2))

    intra_df = intra_df.loc[intra_df.other_mean_cons >= 0.2]
    print('Number interations showing some kind of conservation (>= 0.2): {}'.format(len(intra_df)))

    intra_df_filt = intra_df.drop_duplicates(subset=['merged_name'])
    print('Number of sno_host involved: {}'.format(len(intra_df_filt)), end='\n\n')

    # intra_df['2bit'] =''
    # for i in intra_df.index:
    #     chr = intra_df.at[i, 'chr1']
    #     start = min(intra_df.at[i, 'start1'], intra_df.at[i, 'start2'])
    #     end = max(intra_df.at[i, 'end1'], intra_df.at[i, 'end2'])
    #     name = intra_df.at[i, 'merged_name']
    #
    #     intra_df.at[i, '2bit'] = f'twoBitToFa ../chrom/chr{chr}.2bit:chr{chr}:{start}-{end} {name}'
    #     print(name, f'twoBitToFa ../chrom/chr{chr}.2bit:chr{chr}:{start}-{end} {name}.fa')
    #
    #     sno_start = intra_df.at[i, 'sno_start']
    #     sno_end = intra_df.at[i, 'sno_end']
    #     print(name, f'snoRNA ----->twoBitToFa ../chrom/chr{chr}.2bit:chr{chr}:{sno_start}-{sno_end} {name}.fa')

    intra_df.other_len_cons = intra_df.other_len_cons.map(int)
    intra_df.other_len_cons = intra_df.other_len_cons.map(int)
    print(intra_df[['name1', 'name2', 'E', 'other_len_cons', 'other_mean_cons',
                    'sno_tpm', 'target_tpm']])

    print('\n----------------------------------------------------------------------')
    print(cons_df.columns)

    # ----------------- gene_host_interactions --------------------


def main():

    df = load_df()
    cons_df = pd.read_csv(cons_file, sep='\t')

    stats(df, cons_df)




if __name__ == '__main__':
    main()
