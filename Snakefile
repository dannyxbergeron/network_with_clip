from os.path import join

configfile: "config.json"

rule all:
    input:
        # full_merge = join(config['path']['processed'],
        #                   config['file']['full_network_clip']),
        # interactions = join(config['path']['cytoscape'],
        #                     config['network_file']['interactions']),
        # edges = join(config['path']['cytoscape'],
        #              config['network_file']['edges']),
        beds = expand(join(config['path']['beds'], 'sorted_clean_{bed}.bedgraph'),
                      bed=config['bedgraphs']),
        intron_tok = expand(join(config['path']['tissues_cell_bg'],
                                        'intron_{bed}.bedgraph'),
                            bed=config['bedgraphs']),
        ext_ratio = join(config['path']['sno_host_data'],
                         config['file']['ext_ratio']),
        dist_num_interactions = '/data/articles/SNORD2_article/svgs/dist_num_interactions.svg',
        intra_sno_hosts = '/data/articles/SNORD2_article/svgs/distance_target_from_snoRNA.svg',
        upstream_vs_downstream = '/data/articles/SNORD2_article/svgs/upstream_vs_downstream.svg',
        cannonical_targets = '/data/articles/SNORD2_article/svgs/cannonical_targets.svg',
        bp_distance = join(config['path']['branchpoint'], config['file']['branchpoint']),
        best_branch_points = join(config['path']['branchpoint'],
                                  config['file']['best_branchpoint']),
        svg_branch_point = '/data/articles/SNORD2_article/svgs/branchpoint_binding.svg',
        sno_intron_coordinates = join(config['path']['sno_intron'],
                                      config['file']['sno_intron_coord']),
        merged_mfe = join(config['path']['sno_intron'], config['file']['merged_mfe']),


rule merge_raw_and_filter:
    input:
        raw_files = expand(join(config['path']['gab_raw'],
                                '{SRR}_sno_interaction_snobothside.coco.csv'),
                           SRR=config['raw_files'].keys())
    output:
        initial_file = join(config['path']['raw'],
                            config['file']['initial'])
    conda:
        "envs/python.yaml"
    script:
        "scripts/merge_raw_files.py"


rule get_single_id:
    """ Get single_id in case of more that one """
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
        "scripts/get_single_ID.py"


rule get_single_dg:
    """ Resolve the multiple match per dg """
    input:
        single_id = join(config['path']['processed'],
                         config['file']['single_id'])
    output:
        single_id_and_dg = join(config['path']['processed'],
                                config['file']['single_id_and_dg'])
    params:
        min_length = 8,
        host_max_offset = 10,
    conda:
        "envs/python.yaml"
    script:
        "scripts/get_single_dg_coco.py"


rule check_overlaps:
    """ Merge the windows together """
    input:
        single_id = join(config['path']['processed'],
                         config['file']['single_id_and_dg']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype'])
    output:
        merged_windows = join(config['path']['processed'],
                              config['file']['merged_windows'])
    conda:
        "envs/python.yaml"
    script:
        "scripts/check_overlaps.py"


rule add_offsets:
    """ Add the offsets (in the host gene or the target host/sno) for the
        sno and the targets """
    input:
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype']),
        merged_windows = join(config['path']['processed'],
                              config['file']['merged_windows'])
    output:
        merged_with_offset = join(config['path']['processed'],
                                  config['file']['merged_offset'])
    conda:
        "envs/python.yaml"
    script:
        "scripts/get_offsets.py"


rule double_and_filter:
    """ Double the entries for snoRNA-snoRNA for further analysis and filter
        for interactions < 8 nt and remove intergenic """
    input:
        merged_windows = join(config['path']['processed'],
                              config['file']['merged_offset']),
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
        intaRNA = RNAplex_analyser(
            "../network_with_clip/data/processed/merged_P-L-S_double_sno_intaRNA.tsv")
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
        snodb = join(config['path']['ref'],
                     config['file']['snoDB']),
    output:
        interactions = join(config['path']['cytoscape'],
                            config['network_file']['interactions']),
        nodes = join(config['path']['cytoscape'],
                     config['network_file']['nodes']),
        edges = join(config['path']['cytoscape'],
                     config['network_file']['edges']),
        mapped_clip_only = join(config['path']['processed'],
                                config['file']['mapped_clip_only']),
    conda:
        "envs/python.yaml"
    script:
        "scripts/create_network_V2_with_ENCODE_data.py"


# Include files
include: "rules/process_ref.smk"

# Include the vizualisation subflow
include: "rules/graphs_and_analysis.smk"

# Include the host_interaction part
include: "rules/host_interactions.smk"

# Include rules for the branch point prediction
include: "rules/branch_point_prediction.smk"

# Include rules for getting mfe for folding
include: "rules/mfe_intron_folding.smk"
