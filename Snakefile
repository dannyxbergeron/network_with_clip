from os.path import join

configfile: "config.json"

rule all:
    input:
        single_id = join(config['path']['processed'],
                         config['file']['single_id']),
        merged_windows = join(config['path']['processed'],
                              config['file']['merged_windows']),
        merged_double = join(config['path']['processed'],
                                config['file']['merged_double']),
        tok = 'data/tmp/tok'


rule get_single_id:
    """ keep the more relevant DG and single_id in case of more that one """
    input:
        initial_file = join(config['path']['raw'],
                            config['file']['initial']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype']),
        man_curated = join(config['path']['raw'],
                           config['file']['man_curated']),
    output:
        single_id = join(config['path']['processed'],
                        config['file']['single_id'])
    params:
        multi = join(config['path']['tmp'],
                        config['file']['single_id'])
    conda:
        "envs/python.yaml"
    script:
        "scripts/get_single_ID_and_single_DG.py"


rule check_overlaps:
    """ Merge the windows together """
    input:
        single_id = join(config['path']['processed'],
                         config['file']['single_id']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype'])
    output:
        merged_windows = join(config['path']['processed'],
                              config['file']['merged_windows'])
    conda:
        "envs/python.yaml"
    script:
        "scripts/check_overlaps.py"

rule double_and_filter:
    """ Double the entries for snoRNA-snoRNA for further analysis and filter
        for interactions < 8 nt and remove intergenic """
    input:
        merged_windows = join(config['path']['processed'],
                              config['file']['merged_windows']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype']),
        snoDB = join(config['path']['ref'],
                     config['file']['snoDB'])
    output:
        merged_double = join(config['path']['processed'],
                                config['file']['merged_double'])
    params:
        min_length = 9
    conda:
        "envs/python.yaml"
    script:
        "scripts/analyse_merged.py"


subworkflow RNAplex_analyser:
    workdir:
        "../RNAplex_analyser"
    configfile:
        "../RNAplex_analyser/config.json"

rule intaRNA:
    input:
        RNAplex_analyser("../network_with_clip/data/processed/merged_P-L-S_double_sno_intaRNA.tsv")
    output:
        tok = 'data/tmp/tok'
    shell:
        "echo 'intaRNA done !' && touch {output.tok}"

include: "rules/process_ref.smk"
