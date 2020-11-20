rule viz:
    input:
        data_file = join(config['path']['processed'],
                         config['file']['merged_double_inta']),
        parsed_gtf = join(config['path']['ref'],
                          config['file']['light_parsed_gtf']),
        multiple_list = join(config['path']['ref'],
                             config['file']['multiple_gene_name']),
        gene_id = join(config['path']['ref'],
                       config['file']['name_id']),
        sno_host = join(config['path']['ref'],
                        config['file']['sno_host']),
        trans_tpm = join(config['path']['count_matrices'],
                         config['file']['kallisto_trans_tpm']),
    output:
        'viz.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/viz_transcripts.py"


rule clip_candidates:
    input:
        full_merged = join(config['path']['processed'],
                           config['file']['full_network_clip'])
    output:
        tmp_candidates = join(config['path']['tmp'],
                              config['file']['clip_candidates']),
        bed_file = join(config['path']['tmp'],
                        config['file']['clip_candidates_bed'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/clip_candidates.py"


rule get_stats:
    input:
        data_file = join(config['path']['processed'],
                         config['file']['merged_double_inta']),
        ext_ratio = join(config['path']['sno_host_data'],
                         config['file']['ext_ratio']),
    output:
        'stats.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/stats_and_graphs.py"


rule pearson_corr:
    input:
        tpm = join(config['path']['count_matrices'],
                   config['file']['coco_tpm']),
        parsed_gtf = join(config['path']['ref'],
                          config['file']['light_parsed_gtf']),
    output:
        pearson_corr = join(config['path']['count_matrices'],
                            config['file']['pearson'])
    threads:
        6
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/multiprocessing_pearson.py"


rule correlation_sno_host:
    input:
        pearson_corr = join(config['path']['count_matrices'],
                            config['file']['pearson']),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        cons = join(config['path']['sno_host_data'],
                    config['file']['cons'])
    output:
        'correlation.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/correlation_sno_host.py"


rule gene_ontology:
    input:
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        cons = join(config['path']['sno_host_data'],
                    config['file']['cons']),
        bio_function = join(config['path']['ref'],
                            config['file']['bio_function']),
    output:
        'gene_ontology.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/gene_ontology.py"


rule graph_merged:
    input:
        data_file = join(config['path']['processed'],
                         config['file']['merged_double_inta']),
    output:
        'graph_merged.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/graph_merged.py"


rule sankey:
    input:
        sno_host = join(config['path']['ref'],
                        config['file']['prot_coding_sno_host']),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        snodb = join(config['path']['ref'],
                     config['file']['snoDB']),
        gene_bed = join(config['path']['ref'],
                        config['file']['gene_bed_biotype']),
        net_sno_host = join(config['path']['sno_host_data'],
                            config['file']['sno_host_data']),
    output:
        svg = '/data/labmeetings/host_interactions/sankey.svg'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/sankey.py"

rule pearson_correlation:
    input:
        tpm = join(config['path']['count_matrices'],
                   config['file']['coco_tpm']),
    output:
        'pearson_corr.tok'
    conda:
        "../envs/python_scikit.yaml"
    script:
        "../scripts/pearson_correlation.py"


rule get_significant_NMD:
    input:
        transId_geneId_geneName = join(config['path']['ref'],
                                       config['file']['transId_geneId_geneName']),
        NMD_tpm = join(config['path']['count_matrices'],
                       config['file']['NMD_tpm_transcript']),
    output:
        sign_nmd_trans = join(config['path']['NMD'],
                              config['file']['sign_trans_nmd']),
    params:
        deseq_dir = join(config['path']['NMD'], 'DESeq2'),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get_all_significants_transcripts.py"


rule NMD_graph:
    input:
        NMD_tpm = join(config['path']['count_matrices'],
                       config['file']['NMD_tpm_transcript']),
        transId_geneId = join(config['path']['ref'],
                              config['file']['transId_geneId_geneName']),
        sign_nmd_trans = join(config['path']['NMD'],
                              config['file']['sign_trans_nmd']),
    output:
        'NMD.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/NMD.py"


rule EIF4A3_analysis:
    input:
        EIF4A3_CLIP = join(config['path']['tmp'],
                           config['file']['EIF4A3_CLIP']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype']),
    output:
        'EIF4A3_analysis.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/EIF4A3_analysis.py"
