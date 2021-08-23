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
        svg = '/data/articles/SNORD2_article/svgs/SNORA12_interactions.svg',
        alternative_splicing_intron = '/data/articles/SNORD2_article/svgs/alternative_splicing_intron.svg',



# Include raw processing
include: "rules/process_raw.smk"

# Include building of the network
include: "rules/network.smk"

# Include files
include: "rules/process_ref.smk"

# Include the host_interaction part
include: "rules/host_interactions.smk"

# Include the vizualisation subflow
include: "rules/graphs_and_analysis.smk"

# Include rules for the branch point prediction
include: "rules/branch_point_prediction.smk"

# Include rules for getting mfe for folding
include: "rules/mfe_intron_folding.smk"
