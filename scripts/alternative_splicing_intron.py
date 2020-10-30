import numpy as np
import pandas as pd
import scipy.stats as stats

import matplotlib.pyplot as plt

from pybedtools import BedTool as bt

parsed_file = snakemake.input.parsed
snodb_file = snakemake.input.snodb
sno_host_file = snakemake.input.cons

out_file = snakemake.output.alt_splice

THRES = 300


def load_df(file):
    df = pd.read_csv(file, sep='\t')
    print(file)
    print(df.columns)
    print('-----------------------------')
    return df


def get_sno_and_host(gtf_df, snodb_df, sno_host_df):

    sno_df = gtf_df.loc[(gtf_df.gene_biotype == 'snoRNA') & (gtf_df.feature == 'gene')]
    prot_cod_df = gtf_df.loc[gtf_df.gene_biotype == 'protein_coding']

    snoRNA_ids = set(sno_df.gene_id)
    protein_coding_set = set(prot_cod_df.gene_id)

    snodb_df = snodb_df.loc[(snodb_df.gene_id_annot2020.isin(snoRNA_ids)) &
                            (snodb_df['host gene id'].isin(protein_coding_set))]

    snodb_dict = dict(zip(snodb_df.gene_id_annot2020, snodb_df['host gene id']))
    # CHANGED !!!!!!!!!!!!!!!!!!!
    snodb_dict.update(dict(zip(sno_host_df.single_id1, sno_host_df.single_id2)))

    colnames = ['seqname', 'start', 'end', 'gene_id', 'gene_name', 'strand']
    sno_df = sno_df.loc[sno_df.gene_id.isin(snodb_dict.keys())][colnames]
    prot_cod_df = prot_cod_df.loc[prot_cod_df.gene_id.isin(snodb_dict.values())]

    return snodb_dict, prot_cod_df, sno_df


def create_introns(df):
    master_list = []
    df_values = df.values
    for i, row in enumerate(df_values):
        chr, start, end, exon_number, exon_id, strand = row
        master_list.append(row)
        if i != len(df) - 1:
            if exon_number < df_values[i+1][3]:
                intron_id = f'intron_{exon_id}'
                if strand == '+':
                    int_start = end
                    int_end = df_values[i + 1][1]
                else:
                    int_start = df_values[i + 1][2]
                    int_end = start
                intron_row = [chr, int_start, int_end,
                              exon_number, intron_id, strand]
                master_list.append(intron_row)
    return pd.DataFrame(master_list, columns=df.columns)


def bedtools(df1, df2):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=True, sorted=False)
    new_cols = ['seqname1', 'start1', 'end1', 'gene_id', 'gene_name', 'strand1',
                'seqname2', 'start2', 'end2', 'exon_number', 'exon_id', 'strand2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    return intersect_df


def get_modulation(int_start, int_end, full_intron_exon_df, bt1_):

    bt1 = bt1_.copy(deep=True)
    sno_start = bt1['start'].values[0]
    sno_end = bt1['end'].values[0]

    bt1['start'] = int_start + 1
    bt1['end'] = int_end - 1

    intersect_df = bedtools(bt1, full_intron_exon_df)

    intersect_df = intersect_df.loc[intersect_df.exon_id.str.contains('intron')]
    min_distance = 1000000
    for start, end in intersect_df[['start2', 'end2']].values:
        if (start == int_start and end != int_end):
            min_distance = min(min_distance, abs(end - sno_end), abs(int_end - sno_end))
        elif (end == int_end and start != int_start):
            min_distance = min(min_distance, abs(start - sno_start),  abs(int_start - sno_start))
    if min_distance != 1000000:
        return True, min_distance
    return False, np.nan


def get_sno_intron(snodb_host_dict, prot_cod_df, sno_df_):

    sno_df = sno_df_.copy(deep=True)
    ref_df = prot_cod_df.copy(deep=True)

    intron_start = []
    intron_end = []
    to_remove = []
    splicing_modulation = []
    distances = []
    for idx in sno_df.index:
        sno_id = sno_df.at[idx, 'gene_id']
        sno_name = sno_df.at[idx, 'gene_name']
        host_id = snodb_host_dict[sno_id]
        host_df = ref_df.loc[ref_df.gene_id == host_id]
        sno_data_start = sno_df.at[idx, 'start']
        sno_data_end = sno_df.at[idx, 'end']

        tmp = host_df.loc[(host_df.feature == 'transcript') &
                          (host_df.transcript_biotype == 'protein_coding')]

        if len(tmp) == 0:
            # no protein_coding transcript for this gene...
            to_remove.append(sno_id)
            continue

        tmp = tmp.sort_values('transcript_name')

        tmp_values = tmp[['start', 'end', 'transcript_id', 'transcript_name']].values
        for start, end, transcript_id, transcript_name in tmp_values:
            if sno_data_start > start and sno_data_end < end:
                host_transcript_name = transcript_name
                host_transcript_id = transcript_id
                break
        else:
            # snoRNA not in a protein_coding transcript...
            to_remove.append(sno_id)
            # print(tmp[['gene_name']].values[0][0])
            # print('---------------------------------')
            continue

        exon_df = host_df.loc[(host_df.transcript_id == host_transcript_id)
                              & (host_df.feature == 'exon')]
        exon_df = exon_df[['seqname', 'start', 'end', 'exon_number', 'exon_id', 'strand']]
        exon_intron_df = create_introns(exon_df)

        bt1 = sno_df.loc[sno_df.index == idx]
        bt2 = exon_intron_df

        intersect_df = bedtools(bt1, bt2)

        # snoRNA not the same strand as the host...
        if len(intersect_df) < 1:
            to_remove.append(sno_id)
            continue

        # snoRNA in an exon
        if 'intron' not in intersect_df['exon_id'].values[0]:
            to_remove.append(sno_id)
            continue

        # snoRNA partly in an exon
        if len(intersect_df) > 1:
            to_remove.append(sno_id)
            continue

        int_start = intersect_df.start2.values[0]
        int_end = intersect_df.end2.values[0]

        intron_start.append(int_start)
        intron_end.append(int_end)

        full_exon_df = host_df.loc[host_df.feature == 'exon']
        full_exon_df = full_exon_df[['seqname', 'start', 'end', 'exon_number', 'exon_id', 'strand']]
        full_intron_exon_df = create_introns(full_exon_df)

        s_module, distance = get_modulation(int_start, int_end, full_intron_exon_df, bt1)
        distances.append(distance)
        splicing_modulation.append(s_module)

    to_remove_set = set(to_remove)
    filt_sno_df = sno_df.loc[~sno_df.gene_id.isin(to_remove_set)].copy(deep=True)
    filt_sno_df['host'] = filt_sno_df.gene_id.map(snodb_host_dict)
    filt_sno_df['intron_start'] = intron_start
    filt_sno_df['intron_end'] = intron_end
    filt_sno_df['splicing_hypothesis'] = splicing_modulation
    filt_sno_df['distance'] = distances

    print('==================================')
    print('len original: {}, len filtered: {}'.format(len(sno_df), len(filt_sno_df)))
    original = sno_df.copy(deep=True)
    original['host'] = original.gene_id.map(snodb_host_dict)
    print('Original nb of hosts:{}, filtered nb: {}'.format(len(set(original.host)),
                                                            len(set(filt_sno_df.host))))
    print('==================================')

    return filt_sno_df


def get_stats(sno_df, sno_host_df):

    def get_ratio(df, tresh):
        tmp = df.copy(deep=True)
        tmp_pos = tmp.loc[(tmp.distance <= tresh) & (tmp.splicing_hypothesis)]
        all_true = list(tmp_pos["splicing_hypothesis"]).count(True)
        all = len(tmp)
        ratio = all_true / all
        return all, all_true, ratio

    def fisher(net_true, nb_net, other_true, nb_other):
        net_neg = nb_net - net_true
        other_neg = nb_other - other_true

        obs = np.array([[net_true, net_neg], [other_true, other_neg]])
        return stats.fisher_exact(obs)

    sno_df_in_data = sno_df.loc[sno_df.gene_id.isin(sno_host_df.single_id1)]
    sno_df_not_in_data = sno_df.loc[~sno_df.gene_id.isin(sno_host_df.single_id1)]

    master_list = []
    threshold = [75, 100, 150, 300, 500, 1000, 10000]
    print(sno_df)
    for thres in threshold:

        nb_all, all_true, all_ratio = get_ratio(sno_df, thres)
        nb_net, net_true, net_ratio = get_ratio(sno_df_in_data, thres)
        nb_other, other_true, other_ratio = get_ratio(sno_df_not_in_data, thres)

        print(f'------------------ {thres} -------------------')
        print(f'All snoRNA ratio: {all_ratio:.2f} ({all_true}/{nb_all})')
        print(f'snoRNA inteacting with their host ratio: {net_ratio:.2f} ({net_true}/{nb_net})')
        print(f'Other snoRNA ratio: {other_ratio:.2f} ({other_true}/{nb_other})\n')
        print('------> Fisher exact test:', end=' ')
        odds, p_val = fisher(net_true, nb_net, other_true, nb_other)
        print(f'Odds: {odds:.2f}, pValue: {p_val:.5f}')

        master_list.append([net_true, nb_net, other_true, nb_other, thres])

    return pd.DataFrame(master_list,
                        columns=['net_true', 'nb_net', 'other_true', 'nb_other', 'thres'])


def graph(df):

    # Data
    r = list(range(len(df)))

    fig, axes = plt.subplots(1, 2, figsize=(12, 8))
    barWidth = 0.85
    names = [str(x) for x in df.thres]
    colnames = [('net_true', 'nb_net'), ('other_true', 'nb_other')]
    colors = ['blue', 'red']
    titles = ['snoRNA interacting with their host intron\nwith possible alternative splicing',
              'Other snoRNA with a possible alternative splicing']
    for (pos, nb), ax, color, title in zip(colnames, axes, colors, titles):

        pos = df[pos] / df[nb] * 100
        rest = 100 - pos

        ax.bar(r, pos, color=color, edgecolor='white', width=barWidth, label='Positive')
        ax.bar(r, rest, bottom=pos, color='grey',
               edgecolor='white', width=barWidth, label='rest')
        ax.set_title(title)

        # Custom x axis
        ax.set_xticks(r)
        ax.set_xticklabels(names)

        ax.set_xlabel("Distance of the splicing event from the snoRNA (pb)")
        ax.set_ylabel("% of snoRNA")

        ax.legend()

    # Show graphic
    plt.legend()
    plt.show()


def graph_cumsum(df_, net_sno):

    df_copy = df_.copy(deep=True)
    df_copy.sort_values('distance', inplace=True)
    net_df = df_copy.loc[df_copy.gene_id.isin(net_sno)]
    other_df = df_copy.loc[~df_copy.gene_id.isin(net_sno)]

    print(len(net_df), len(other_df))

    fig, axes = plt.subplots(figsize=(12, 8))
    dfs = [net_df, other_df]
    colors = ['blue', 'red']
    labels = ['snoRNA-host interacting', 'snoRNA-host not interacting']

    for df, color, label in zip(dfs, colors, labels):
        size = len(df)
        value = 100 / size

        x = [0]
        y = [0]
        for dist in df.distance.values:
            if dist == -1:
                continue
            if dist > 2000:
                break

            y.append(y[-1])
            x.append(dist)
            x.append(dist)
            y.append(y[-1] + value)

        plt.plot(x, y, color=color, label=label, linewidth=3)

    plt.grid(b=True, which='major', color='lightgray', linestyle='-')
    plt.title('snoRNA splicing modulation potential')
    plt.xlabel('Distance from the closest splicing event')
    plt.ylabel('Cumulative % of snoRNA')
    plt.legend()
    plt.show()


def main():

    gtf_df = load_df(parsed_file)
    snodb_df = load_df(snodb_file)
    sno_host_df = load_df(sno_host_file)
    sno_host_df = sno_host_df.loc[sno_host_df.interaction_type == 'intra']

    snodb_host_dict, prot_cod_df, sno_df = get_sno_and_host(gtf_df, snodb_df, sno_host_df)

    sno_df = get_sno_intron(snodb_host_dict, prot_cod_df, sno_df)

    # Get stats to make a graph
    stats_df = get_stats(sno_df, sno_host_df)
    graph(stats_df)

    graph_cumsum(sno_df, list(sno_host_df.single_id1))

    print('===============================================================')
    name_id_dict = dict(zip(prot_cod_df.gene_id, prot_cod_df.gene_name))
    sno_df['host_name'] = sno_df['host'].map(name_id_dict)
    # sno_df = sno_df.loc[sno_df.distance < 500]
    sno_df['in_net'] = np.where(sno_df.gene_id.isin(sno_host_df.single_id1),
                                True,
                                False)
    print(sno_df[['seqname', 'start', 'end', 'gene_name',
                  'host_name', 'splicing_hypothesis', 'distance', 'in_net']])

    splicing_dict = dict(zip(sno_df.gene_id, sno_df.splicing_hypothesis))
    distance_dict = dict(zip(sno_df.gene_id, sno_df.distance))
    int_start_dict = dict(zip(sno_df.gene_id, sno_df.intron_start))
    int_end_dict = dict(zip(sno_df.gene_id, sno_df.intron_end))

    sno_host_df['splicing_hypothesis'] = sno_host_df.single_id1.map(splicing_dict)
    sno_host_df['splice_dist'] = sno_host_df.single_id1.map(distance_dict)
    sno_host_df['intron_start'] = sno_host_df.single_id1.map(int_start_dict)
    sno_host_df['intron_end'] = sno_host_df.single_id1.map(int_end_dict)

    sno_host_df.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
