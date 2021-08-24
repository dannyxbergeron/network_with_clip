from io import StringIO
from os.path import join

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import subprocess
import shlex

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']


input_file = snakemake.input.coords
fasta_dir = snakemake.params.fasta_dir
svg = snakemake.output.svg


def intaRNA(df):

    def get_pred(cmd):
        cmd = shlex.split(cmd)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout = process.communicate()[0].decode('utf-8')
        data = StringIO(stdout)
        tmp = pd.read_csv(data, sep=',')
        if tmp.empty:
            row = [['dummy1', 0, 0, 'dummy2', 0, 0, np.nan, np.nan, np.nan]]
            return pd.DataFrame(row, columns=tmp.columns)
        return tmp

    cols = ['id1', 'start1', 'end1', 'id2', 'start2', 'end2', 'subseqDP', 'hybridDP', 'E', 'DG', 'type']
    master_df = pd.DataFrame([], columns=cols)
    for dg in df.DG.values:
        dg_path = join(fasta_dir, dg)
        for type in ['target', 'simulated']:
            cmd = f'IntaRNA -q {dg_path}-snoRNA.fa -t {dg_path}-{type}.fa --outMode C --outSep=,'
            tmp = get_pred(cmd)
            tmp['DG'] = dg
            tmp['type'] = type
            master_df = pd.concat([master_df, tmp])

    return master_df

def stats(df):

    def gb(tmp):
        type = tmp.type.values[0]
        null = tmp.loc[tmp.E.isnull()]
        print(f'Number of {type} null E values: {len(null)}')
        not_null = tmp.loc[~(tmp.E.isnull())]
        print(f'Number of {type} non-null E values: {len(not_null)}')

    target = df.loc[df.type == 'target']
    simulated = df.loc[df.type == 'simulated']

    print('\n=================== STATS =====================')
    gb(target)
    gb(simulated)
    print('-------------------- END ---------------------\n')


def graph(df_):

    df = df_.copy(deep=True)
    # df = df.fillna(0)
    df = df.dropna()

    target = df.loc[df.type == 'target'].E
    simulated = df.loc[df.type == 'simulated'].E
    data = [target, simulated]
    types = [
        'Target',
        'Simulated',
    ]
    colors = ['#4daf4a', '#e41a1c', '#377eb8']

    fig, ax = plt.subplots(figsize=(10,8))
    fig.canvas.draw()

    for label, data, color in zip(types, data, colors):
        sns.kdeplot(data=data, shade=True, linewidth=1, alpha=.3,
                    label=label, ax=ax, bw=1,
                    color=color)

    for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)


    plt.title('IntaRNA prediction\n(NaN removed, 60 for target and 63 for simulated)', fontsize=25)
    plt.xlabel('Molecular free energy (mfe)', fontsize=20)
    plt.ylabel('Distribution', fontsize=20)
    plt.legend(fontsize=16)
    # plt.savefig(out_file, format='svg')
    plt.show()



def main():

    df = pd.read_csv(input_file, sep='\t')

    intaRNA_df = intaRNA(df)

    stats(intaRNA_df)


    graph(intaRNA_df)



if __name__ == '__main__':
    main()
