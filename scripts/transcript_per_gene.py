from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu as mwu

import matplotlib.pyplot as plt
import seaborn as sns

parsed_file = snakemake.input.parsed
snodb_file = snakemake.input.snodb
sno_host_file = snakemake.input.sno_host

def load_df(file):
    df = pd.read_csv(file, sep='\t')
    return df


def analyse(gtf_df, snodb_df, sno_host_df):

    gtf_df = gtf_df.loc[(gtf_df.feature.isin(['transcript']))
                        & (gtf_df.gene_biotype == 'protein_coding')]
    print('Number of total transctipts', len(gtf_df))
    print('Number of total protein_coding genes', len(set(gtf_df.gene_id)))
    print('----------------------------------------------------------')

    host_list = set(snodb_df['host gene id'])
    host_list = [x for x in host_list if not pd.isna(x)]

    host_df = gtf_df.loc[gtf_df.gene_id.isin(host_list)]
    non_host_df = gtf_df.loc[~(gtf_df.gene_id.isin(host_list))]
    # sno_host_df = sno_host_df.loc[sno_host_df.interaction_type == 'intra']
    sno_host_list = sno_host_df.single_id2
    network_host_df = host_df.loc[host_df.gene_id.isin(sno_host_list)]

    # host_df = host_df.loc[~(host_df.gene_id.isin(network_host_df.gene_id))] # CHANGED
    # snoRNA not having the same strand as the host ?? TODO !!!!
    host_len_trans = len(host_df)
    host_len_genes = len(set(host_df.gene_id))

    non_host_len_trans = len(non_host_df)
    non_host_len_genes = len(set(non_host_df.gene_id))

    network_len_trans = len(network_host_df)
    network_len_genes = len(set(network_host_df.gene_id))


    print('Total transcript for other genes ', non_host_len_trans)
    print('Total genes for other genes ', non_host_len_genes)
    print('----------------------------------------------------------')
    print('Total transcript for host_genes ', host_len_trans)
    print('Total genes for host_genes ', host_len_genes)
    print('----------------------------------------------------------')
    print('Total transcript for network genes', network_len_trans)
    print('Total genes for other network host genes', network_len_genes)
    print('=======================================\n')

    print('Host genes transcript per gene: {:.2f}'.format(host_len_trans/host_len_genes))
    print('Other genes transcript per gene: {:.2f}'.format(non_host_len_trans/non_host_len_genes))
    print('Host that the snoRNA interacts with: {:.2f}'.format(network_len_trans/network_len_genes))



    all_host = list(Counter(list(host_df.gene_id)).values())
    all_non_host = list(Counter(list(non_host_df.gene_id)).values())
    network_host = list(Counter(list(network_host_df.gene_id)).values())

    nb_genes = [non_host_len_genes, host_len_genes, network_len_genes]

    return all_non_host, all_host, network_host, nb_genes

def graph(data, nb_genes):

    MAX_VAL = 60

    # ============= STATS ==================
    print('\n=========== STATS - mann-whitney u test ==========')
    all_stats, all_pval = mwu(data[1], data[0], alternative='two-sided')
    print(f'For all_host vs all_non_host p_value: {all_pval}')
    net_stats, net_pval = mwu(data[2], data[1], alternative='two-sided')
    print(f'For network host interacting vs all_host p_value: {net_pval}')
    print('===================================================\n')
    # ============= STATS ==================


    groups = [
        f'all_non_host ({nb_genes[0]} genes)',
        f'all_snoRNA_hosts ({nb_genes[1]} genes)',
        f'network_snoRNA_host ({nb_genes[2]} genes)'
    ]

    for i in range(len(data)):
        data[i] = [x if x <=MAX_VAL else MAX_VAL for x in data[i]]

    fig, ax = plt.subplots()
    fig.canvas.draw()

    for label, data in zip(groups[:2], data[:2]):
        ax = sns.distplot(data, hist = False, kde = True,
                     kde_kws = {'shade': True, 'linewidth': 1},
                     label=label, ax=ax)

    tick_labels = [
        int(tick_label)
        for tick_label in ax.get_xticks().tolist()
    ]

    tick_labels[-2] = str(tick_labels[-2]) + '+'
    ax.set_xticklabels(tick_labels)

    plt.title('Transcripts/gene for snoRNA host genes')
    plt.xlabel('Number of transcripts per gene')
    plt.legend()
    plt.show()


def main():

    gtf_df = load_df(parsed_file)
    snodb_df = load_df(snodb_file)
    sno_host_df = load_df(sno_host_file)

    all, all_non, network, nb_genes = analyse(gtf_df, snodb_df, sno_host_df)

    graph([all, all_non, network], nb_genes)


if __name__ == '__main__':
    main()
