rule host_interacting:
    input:
        data_file = join(config['path']['processed'],
                         config['file']['merged_double_inta']),
        parsed_gtf = join(config['path']['ref'],
                          config['file']['light_parsed_gtf'])
    output:
        sno_host = join(config['path']['sno_host_data'],
                        config['file']['sno_host_data']),
        host_ref = join(config['path']['sno_host_data'],
                        config['file']['host_ref']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/host_interaction_analysis.py"


# rule create_bedgraph_cons:
#     input:
#         phastconst = '/data/heavy_stuff/phastconst/phastCons100way.BedGraph.gz',
#         host_ref = join(config['path']['sno_host_data'],
#                         config['file']['host_ref']),
#     output:
#         bedgraph = join(config['path']['sno_host_data'],
#                         config['file']['host_bedgraph']),
#     conda:
#         "../envs/python.yaml"
#     shell:
#         "zcat {input.phastconst} | scripts/createConstBedgraph {input.host_ref} > "
#         "{output.bedgraph}"


rule conservation:
    input:
        sno_host = join(config['path']['sno_host_data'],
                        config['file']['sno_host_data']),
        host_ref = join(config['path']['sno_host_data'],
                        config['file']['host_ref']),
        bedgraph = join(config['path']['sno_host_data'],
                        config['file']['host_bedgraph']),
        tpm = join(config['path']['count_matrices'],
                   config['file']['coco_tpm'])
    output:
        cons = join(config['path']['sno_host_data'],
                    config['file']['cons'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/conservation.py"


rule transcript_per_gene:
    input:
        parsed = join(config['path']['ref'],
                      config['file']['parsed_gtf']),
        snodb = join(config['path']['ref'],
                     config['file']['snoDB']),
        sno_host = join(config['path']['sno_host_data'],
                        config['file']['sno_host_data']),
    output:
        'transcript_per_gene.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/transcript_per_gene.py"

rule alternative_splicing_intron:
    input:
        parsed = join(config['path']['ref'],
                      config['file']['light_parsed_gtf']),
        snodb = join(config['path']['ref'],
                     config['file']['snoDB']),
        sno_host = join(config['path']['sno_host_data'],
                        config['file']['sno_host_data']),
    output:
        'alternative_splicing_intron.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/alternative_splicing_intron.py"
