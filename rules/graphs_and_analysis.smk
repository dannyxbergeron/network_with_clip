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
