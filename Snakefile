from os.path import join

configfile: "config.json"

rule all:
    input:
        full_merge = join(config['path']['processed'],
                          config['file']['full_network_clip']),
        interactions = join(config['path']['cytoscape'],
                            config['network_file']['interactions']),


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
        merged_double = join(config['path']['processed'],
                             config['file']['merged_double']),
        intaRNA = RNAplex_analyser("../network_with_clip/data/processed/merged_P-L-S_double_sno_intaRNA.tsv")
    output:
        intaRNA_tok = 'data/tmp/intaRNA_tok',
        merged_double_inta = join(config['path']['processed'],
                                  config['file']['merged_double_inta']),
    shell:
        "echo 'intaRNA done !' && touch {output.intaRNA_tok}"

subworkflow clip_analysis:
    workdir:
        "../Tommy_stuff/dan_new_analysis"

rule process_clip:
    input:
        merge_clip = clip_analysis("data_processed/all_merged.bed")
    output:
        merge_clip = join(config['path']['clip'],
                          config['file']['clip_data'])
    shell:
        "cp {input.merge_clip} {output.merge_clip}"


rule merge_network_and_clip:
    input:
        merged_double_inta = join(config['path']['processed'],
                                  config['file']['merged_double_inta']),
        merged_clip = join(config['path']['clip'],
                           config['file']['clip_data'])
    output:
        full_merge = join(config['path']['processed'],
                          config['file']['full_network_clip'])
    conda:
        "envs/python.yaml"
    script:
        "scripts/network_analysis.py"

rule build_network:
    input:
        full_merge = join(config['path']['processed'],
                          config['file']['full_network_clip']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype']),
    output:
        interactions = join(config['path']['cytoscape'],
                            config['network_file']['interactions']),
        nodes = join(config['path']['cytoscape'],
                     config['network_file']['nodes']),
        mapped_clip_only = join(config['path']['processed'],
                                config['file']['mapped_clip_only']),
    conda:
        "envs/python.yaml"
    script:
        "scripts/create_network_V2_with_ENCODE_data.py"



# Include files
include: "rules/process_ref.smk"
