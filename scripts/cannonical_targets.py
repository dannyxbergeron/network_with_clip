import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

network_file = snakemake.input.merged_double
sno_host_loc_file = snakemake.input.sno_host_loc
sno_host_file = snakemake.input.sno_host
snodb_file = snakemake.input.snodb
bio_function_file = snakemake.input.bio_function

out_can = snakemake.output.svg_can_targets
out_bio_func = snakemake.output.svg_bio_functions
out_box_type = snakemake.output.svg_box_type

def load_df(file):
    df = pd.read_csv(file, sep='\t')
    return df

def get_ids(df, network_df, loc_df):

    intra = set(df.loc[df.interaction_type == 'intra'].single_id1)
    diff_loc = set(df.loc[~(df.interaction_type == 'intra')].single_id1)

    others = set(network_df.single_id1)
    others = others - intra - diff_loc
    others = others.intersection(set(loc_df.gene_id))

    print('--> Length of intra, diff_loc and others')
    print(len(intra), len(diff_loc), len(others))
    print('-------------------------\n')

    return intra, diff_loc, others

def analyse_can_targets(intra, diff_loc, others, snodb_df):

    def get_percent(ids, can_sno):
        ids_set = set(ids)
        can_sno_set = set(can_sno)
        total = len(ids_set)
        intersect = ids_set.intersection(can_sno_set)
        return len(intersect) / total * 100


    snodb_df = snodb_df[[
        'gene_id_annot2020', 'rrna', 'snrna'
    ]].copy(deep=True)
    snodb_df.dropna(thresh=2, inplace=True)
    cannonical_snoRNA = list(snodb_df.gene_id_annot2020)

    intra_percent = get_percent(intra, cannonical_snoRNA)
    diff_loc_percent = get_percent(diff_loc, cannonical_snoRNA)
    others_percent = get_percent(others, cannonical_snoRNA)

    print(f'Intra % of cannonical: {intra_percent:.1f}%')
    print(f'Diff loc % of cannonical: {diff_loc_percent:.1f}%')
    print(f'Others % of cannonical: {others_percent:.1f}%')

    return ['Same intron', 'Not same intron', 'Others'], [intra_percent, diff_loc_percent, others_percent]

def graph_cannonical(names, data):

    total = [[100] for x in data]
    data = [[x] for x in data]

    fig, ax = plt.subplots(figsize=(8,8))

    bar_tot = sns.barplot(data=total, color='#fc8d62')
    bar_can = sns.barplot(data=data, color='#66c2a5')

    tick_labels = [
        label
        for label, tick_label in zip(names, ax.get_xticks().tolist())
    ]
    ax.set_xticklabels(tick_labels)

    plt.title('Number of cannonical snoRNA in the different groups', size=15)
    plt.xlabel('Groups', size=12)
    plt.ylabel('Percentage of snoRNAs', size=12)

    top_bar = mpatches.Patch(color='#fc8d62', label='Orphans')
    bottom_bar = mpatches.Patch(color='#66c2a5', label='Cannonical snoRNAs')
    plt.legend(handles=[top_bar, bottom_bar])

    # plt.show()
    plt.savefig(out_can, format='svg')
    plt.close()

def prepare_data(intra, diff_loc, others, loc_df, bio_func_df):

    def put_in_df(ids, name):
        tmp = pd.DataFrame(ids, columns=['gene_id'])
        tmp['group'] = name
        return tmp

    intra_df = put_in_df(intra, 'intra')
    diff_loc_df = put_in_df(diff_loc, 'diff_loc')
    others_df = put_in_df(others, 'others')
    master_df_original = pd.concat([intra_df, diff_loc_df, others_df])

    master_df = master_df_original.copy(deep=True)
    master_df['host_id'] = master_df.gene_id.map(dict(zip(loc_df.gene_id,
                                                                   loc_df.host_id)))
    master_df['host_function'] = master_df.host_id.map(dict(zip(bio_func_df.host_id,
                                                                bio_func_df.host_function)))
    master_df['host_function'] = master_df['host_function'].fillna('Not investigated')

    master_df_gb = master_df[['group', 'host_function', 'gene_id']]\
                        .groupby(['group', 'host_function'])\
                        .count()\
                        .reset_index()
    total_pergroup = master_df_gb[['group', 'gene_id']].groupby('group').sum().reset_index()
    master_df_gb['sum_group'] = master_df_gb.group.map(dict(zip(total_pergroup.group,
                                                                total_pergroup.gene_id)))
    master_df_gb['percentage_func'] = master_df_gb['gene_id'] / master_df_gb['sum_group'] * 100
    master_df_gb = master_df_gb.round({'percentage_func': 1})

    return master_df_gb, master_df_original

def graph_bio_functions(data_df):

    groups = ['intra', 'diff_loc', 'others']
    new_group_names = ['Same intron', 'Not same intron', 'Others']
    # host_fct = sorted(list(set(data_df.host_function)), reverse=True)
    host_fct = [
        'ribosomal protein',
        'Ribosome biogenesis & translation',
        'RNA binding, processing, splicing',
        'poorly characterized',
        'Other',
        'Not investigated'
    ]
    colors = [
        '#80b1d3',
        '#fdb462',
        '#fb8072',
        '#8dd3c7',
        '#bebada',
        '#d9d9d9',
        '#ffffb3',
    ]


    print(host_fct)

    fig, ax = plt.subplots(1, figsize=(10, 8))
    x = np.arange(len(groups))  # the label locations
    width = 0.9  # the width of the bars

    bottom = np.zeros(len(groups))
    for func, color in zip(host_fct, colors):
        tmp = data_df.loc[data_df.host_function == func]
        tmp_dict = dict(zip(tmp.group, tmp.percentage_func))
        tmp_data = [
            tmp_dict[group]
            if group in tmp_dict
            else 0
            for group in groups
        ]
        label = func.capitalize() if not func.startswith('RNA') else func
        ax.bar(x, tmp_data, width, label=label, bottom=bottom, color=color)
        bottom += np.array(tmp_data)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Percentage of host in fuction categories', size=12)
    ax.set_xlabel('Groups', size=12)
    ax.set_title('Proportions of function of snoRNA host gene according to their locations', size=15)
    ax.set_xticks(x)
    ax.set_xticklabels(new_group_names)

    # ax.legend()
    plt.legend(bbox_to_anchor=([1, 1, 0, 0]))

    fig.subplots_adjust(right=0.7)
    plt.savefig(out_bio_func, format='svg')
    # plt.show()
    plt.close()

def process_box_type(merged_df_):
    merged_df = merged_df_.copy(deep=True)
    merged_df.dropna(subset=['box_type'], inplace=True) # 2 scaRNA
    data_df = merged_df.groupby(['group', 'box_type']).count().reset_index()
    total = data_df.groupby('group').sum().reset_index()
    data_df['total_pergroup'] = data_df.group.map(dict(zip(total.group, total.gene_id)))
    data_df['percent_box_type'] = data_df.gene_id / data_df.total_pergroup * 100
    data_df = data_df.round({'percent_box_type': 1})

    return data_df

def graph_box_type(data_df_):

    order_dict = {'intra': 0, 'diff_loc': 1, 'others': 2}
    data_df = data_df_.copy(deep=True)
    data_df['order'] = data_df.group.map(order_dict)
    data_df.sort_values(['order'], inplace=True)
    labels = ['Same intron', 'Not same intron', 'Others']

    print(data_df)

    total = data_df.groupby('group').sum().reset_index()
    cd = data_df.loc[data_df.box_type == 'C/D']

    fig, ax = plt.subplots(figsize=(8,8))

    bar_tot = sns.barplot(x='group', y='percent_box_type', data=total, color='#7570b3')
    bar_can = sns.barplot(x='group', y='percent_box_type', data=cd, color='#d95f02')

    tick_labels = [
        label
        for label, tick_label in zip(labels, ax.get_xticks().tolist())
    ]
    ax.set_xticklabels(tick_labels)

    plt.title('Proportion of snoRNA box type in different groups', size=15)
    plt.xlabel('Groups', size=12)
    plt.ylabel('Percentage of snoRNAs', size=12)

    top_bar = mpatches.Patch(color='#7570b3', label='H/ACA')
    bottom_bar = mpatches.Patch(color='#d95f02', label='C/D')
    plt.legend(handles=[top_bar, bottom_bar])

    # plt.show()
    plt.savefig(out_box_type, format='svg')
    # plt.close()

def main():

    df = load_df(sno_host_file)
    loc_df = load_df(sno_host_loc_file)
    snodb_df = load_df(snodb_file)
    network_df = load_df(network_file)

    print('df', df.columns)
    print('loc_df', loc_df.columns)
    print('snodb_df', snodb_df.columns)
    print('network_df', network_df.columns)
    print('==================================\n')

    intra, diff_loc, others = get_ids(df, network_df, loc_df)

    names, data = analyse_can_targets(intra, diff_loc, others, snodb_df)

    graph_cannonical(names, data)

    # Bio function graph
    print('\n============================== Bio function analysis ==============================\n')
    bio_func_df = load_df(bio_function_file)
    data_df, merged_df = prepare_data(intra, diff_loc, others, loc_df, bio_func_df)

    graph_bio_functions(data_df)

    # Box type analysis
    print('\n============================== Box type analysis ==============================\n')
    merged_df['box_type'] = merged_df.gene_id.map(dict(zip(snodb_df.gene_id_annot2020,
                                                           snodb_df['box type'])))
    data_df = process_box_type(merged_df)
    graph_box_type(data_df)




if __name__ == '__main__':
    main()
