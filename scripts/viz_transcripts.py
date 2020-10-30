from collections import defaultdict
import os

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

plt.rcParams['svg.fonttype'] = 'none'

data_file = snakemake.input.data_file

parsed_gtf = snakemake.input.parsed_gtf
multiple_file = snakemake.input.multiple_list
name_id = snakemake.input.gene_id
sno_host = snakemake.input.sno_host

transcript_tpm = snakemake.input.trans_tpm

COLORS = ['#6a3d9a', '#ff7f00', '#1f78b4', '#33a02c', '#b15928',
          '#a6cee3', '#b2df8a', '#fdbf6f', '#cab2d6', '#fb9a99',
          '#e31a1c', '#ffff99', 'black', 'grey', '#b3de69',
          '#ccebc5', '#fccde5', '#fb9a99']


def validate_name(gene_name):

    def end(gene_name):
        print('Error - {} contains more than one corresponding id'.format(
            gene_name))
        exit()

    multiple = []
    with open(multiple_file, 'r') as f:
        multiple = f.read().splitlines()

    if gene_name in multiple:
        end(gene_name)

    corresp_df = pd.read_csv(name_id, sep='\t')
    corresp_dict = dict(zip(corresp_df.gene_name, corresp_df.gene_id))

    if gene_name in corresp_dict:
        return corresp_dict[gene_name]
    else:
        end(gene_name)


def get_snos(gene_id):

    df = pd.read_csv(sno_host, sep='\t')
    df = df.loc[df.host_id == gene_id]
    return list(df.sno_id)


def load_ref(gene_name):

    df = pd.read_csv(parsed_gtf, sep='\t', dtype={'transcript_support_level': float})
    df = df.fillna(value={'transcript_support_level': 0})  # CHANGED !!!!

    gene_id = validate_name(gene_name)

    snos_list = get_snos(gene_id)

    snos = df.loc[(df.gene_id.isin(snos_list)) & (df['feature'] == 'gene')]
    snos = snos[['gene_name', 'start', 'end']]

    df = df.loc[df.gene_id == gene_id]

    df = df[['seqname', 'feature', 'start', 'end', 'strand', 'exon_id',
             'exon_number', 'gene_biotype', 'gene_name', 'gene_id',
             'transcript_id', 'transcript_name', 'transcript_biotype', 'transcript_support_level']]
    df = df.loc[(df['feature'] == 'gene') | (df['transcript_support_level'] <= 5)]
    df = df.loc[df.feature.isin(['gene', 'exon', 'transcript'])]

    strand = df.values[0][4]

    return df, snos, strand


def get_data(gene_name):

    df = pd.read_csv(data_file, sep='\t')
    df = df.loc[df.name2 == gene_name]
    # df = df.loc[df['E'] < 0] #CHANGED
    # df = df[['start2', 'end2', 'name1', 'E']] # CHANGED!!!!
    df.reset_index(inplace=True, drop=True)

    return df


def add_sno_data(sno_data, ax, gene_start, gene_end, yloc, BOX_HEIGHT, off):

    def getHeight(energie, b_h):
        h = energie / -15 * b_h
        if h < 0.2 or str(energie) == "nan":
            # return 0.3
            return 0.5  # CHANGED
        elif h > 1:
            return 1
        return h

    corr = yloc * 0.05
    legloc = yloc + 0.25 * BOX_HEIGHT
    new_off = off / 2
    for i in sno_data.index:
        start = sno_data.at[i, 'start2'] - gene_start
        width = sno_data.at[i, 'end2'] - sno_data.at[i, 'start2']
        height = getHeight(sno_data.at[i, 'E'], BOX_HEIGHT)
        sno_name = sno_data.at[i, 'name1']

        # adjust = gene_end * 0.00075  # just to be able to see the boxes
        adjust = 0
        rect = Rectangle((start-adjust, yloc-(height / 2)),
                         width + (2 * adjust),
                         height)

        rect_leg = Rectangle((gene_end + new_off, legloc - corr/2),
                             new_off,
                             BOX_HEIGHT / 8 + corr)
        ax.text(gene_end + off*1.2, legloc + BOX_HEIGHT * 0.01, sno_name)

        if sno_data.at[i, 'host_interaction']:
            tar_start = sno_data.at[i, 'start1'] - gene_start
            tar_width = sno_data.at[i, 'end1'] - sno_data.at[i, 'start1']
            tar_height = height / 2

            tar_rect = Rectangle((tar_start-adjust, yloc - height/2),
                                 tar_width + (2 * adjust),
                                 height / 2)
            pc = PatchCollection([tar_rect], facecolor=COLORS[i], alpha=.5,
                                 edgecolor=None, zorder=300)
            ax.add_collection(pc)

        pc = PatchCollection([rect, rect_leg], facecolor=COLORS[i], alpha=1,
                             edgecolor=None, zorder=200)

        ax.add_collection(pc)
        legloc -= 0.25 * BOX_HEIGHT + corr


def draw_sno_in_host(ax, start_, t_start, t_end,
                     sno_in_host, yloc, BOX_HEIGHT, high):
    snos = []
    for i in sno_in_host.index:
        w = sno_in_host.at[i, 'end'] - sno_in_host.at[i, 'start']
        s = sno_in_host.at[i, 'start'] - start_
        name = sno_in_host.at[i, 'gene_name'].replace('SNOR', '')

        if s >= t_start and w + s <= t_end:
            rect = Rectangle((s, yloc-BOX_HEIGHT/4), w, BOX_HEIGHT/2)
            snos.append(rect)
            pc = PatchCollection(snos, facecolor='#fd8d3c', alpha=1,
                                 edgecolor=None, zorder=100)
            ax.add_collection(pc)
            ax.text(s,
                    yloc-(BOX_HEIGHT * (0.2714 + high * 0.0142857)),  # (high / 30)
                    name,
                    fontname='Comic sans ms',
                    fontsize=6)


def graph_transcripts(ax, transcripts, BOX_HEIGHT, STEPS):

    yloc = 1

    df = pd.read_csv(transcript_tpm, sep='\t')

    transcripts_dict = {t: i for (i, t) in enumerate(transcripts)}
    df = df.loc[df.transcript.isin(transcripts)]
    df['order'] = df.transcript.map(transcripts_dict)
    df.sort_values(by='order', ascending=True, inplace=True)
    df = df.drop(columns=['order'])

    df['max'] = df.max(axis=1)
    df['avg'] = df.mean(axis=1)
    df['std'] = df.std(axis=1)
    print(df[['transcript', 'max', 'avg', 'std']])

    rectangles = []
    max_val = max(df['avg'])
    for i in df.index:
        s = 0
        w = df.at[i, 'avg']
        name = df.at[i, 'transcript']

        rect = Rectangle((s, yloc-BOX_HEIGHT/2), w, BOX_HEIGHT)
        rectangles.append(rect)
        ax.hlines(yloc, 0, max_val*1.05, color='white', zorder=1)

        ax.text(0, yloc + 0.4 * STEPS, name,
                fontsize=8)

        yloc += STEPS

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(rectangles, facecolor='#33393B')

    # Add collection to axes
    ax.add_collection(pc)

    yloc += STEPS
    ax.set_ylim(0, yloc)
    ax.get_yaxis().set_visible(False)
    ax.set_title('TPM')


def prepare_fig(ref_df, sno_data, sno_in_host, TPM_VALS, strand):

    FONTSIZE = 10
    BOX_HEIGHT = 1
    STEPS = 1.5

    ref_df.reset_index(inplace=True)
    gene = ref_df.loc[ref_df.feature == 'gene']
    start_ = gene.at[0, 'start']
    start = 0
    end = gene.at[0, 'end'] - start_
    gene_name = gene.at[0, 'gene_name']

    # Create figure and axes
    if TPM_VALS:
        fig = plt.figure(constrained_layout=True, figsize=(20, 10))
        gs = fig.add_gridspec(1, 4)
        ax = fig.add_subplot(gs[0, 1:])
        ax2 = fig.add_subplot(gs[0, :1])
    else:
        fig, ax = plt.subplots(figsize=(16, 8))

    offset = 0.05 * end
    ax.set_xlim(start - offset, end + 5*offset)

    tmp = ref_df.loc[ref_df.feature == 'transcript']
    tmp.sort_values(by=['transcript_support_level', 'transcript_id'],
                    ascending=False,
                    inplace=True)
    # print(tmp)
    transcripts = tmp.transcript_id
    # for t in transcripts:
    #     print(t)

    yloc = 1
    exons = []
    high = len(transcripts) + 1
    for t in transcripts:
        print(t)
        tmp = ref_df.loc[ref_df.transcript_id == t]
        tmp.reset_index(inplace=True, drop=True)

        t_id = tmp.at[0, 'transcript_id']
        t_name = tmp.at[0, 'transcript_name']
        biot = tmp.at[0, 'transcript_biotype']
        tsl = tmp.at[0, 'transcript_support_level']
        t_start = tmp.at[0, 'start'] - start_
        t_end = tmp.at[0, 'end'] - start_

        exon_list = [(tmp.at[i, 'start'] - start_,
                      tmp.at[i, 'end'] - tmp.at[i, 'start']) for i in tmp.index[1:]]

        ax.hlines(yloc, t_start, t_end, color='#a6bddb')
        ax.text(0, yloc + 0.4 * STEPS,
                '{} ({}) tsl{} {}'.format(t_id, biot, int(tsl), t_name),
                fontsize=FONTSIZE)

        for e in exon_list:
            s = e[0]
            w = e[1]

            rect = Rectangle((s, yloc-BOX_HEIGHT/2), w, BOX_HEIGHT)
            exons.append(rect)

        # DRAW THE SNORNAs OF THE HOST !!
        draw_sno_in_host(ax, start_, t_start, t_end,
                         sno_in_host, yloc, BOX_HEIGHT, high)

        # Create patch collection with specified colour/alpha
        pc = PatchCollection(exons, facecolor='#2b8cbe', alpha=.95,
                             edgecolor=None, zorder=200)

        # Add collection to axes
        ax.add_collection(pc)

        yloc += STEPS
        exons.clear()

    ax.hlines(yloc, start, end, color='grey')
    ax.text(0, yloc + 0.4 * STEPS, 'snoRNA locations', fontsize=FONTSIZE)

    add_sno_data(sno_data, ax, start_, end, yloc, BOX_HEIGHT, offset)

    yloc += STEPS
    ax.set_ylim(0, yloc)
    ax.get_yaxis().set_visible(False)
    ax.set_title(f'{gene_name} (strand {strand})')

    if TPM_VALS:
        # graph the transcript part
        graph_transcripts(ax2, transcripts, BOX_HEIGHT, STEPS)

    plt.tight_layout()
    # plt.savefig('/data/labmeetings/DKC1/transcripts_DKC1__.svg',
    #             format='svg', transparent=True)
    plt.show()


def main():

    GOI = 'RPS8'
    TPM_VALS = True
    ref_df, sno_in_host, strand = load_ref(GOI)

    sno_data = get_data(GOI)

    prepare_fig(ref_df, sno_data, sno_in_host, TPM_VALS, strand)


if __name__ == "__main__":
    main()
