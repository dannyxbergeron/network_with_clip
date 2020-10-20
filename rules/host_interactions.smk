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


rule create_bedgraph_cons:
    input:
        phastconst = '/data/heavy_stuff/phastconst/phastCons100way.BedGraph.gz',
        host_ref = join(config['path']['sno_host_data'],
                        config['file']['host_ref']),
    output:
        bedgraph = join(config['path']['sno_host_data'],
                        config['file']['host_bedgraph']),
    conda:
        "../envs/python.yaml"
    shell:
        "zcat {input.phastconst} | scripts/createConstBedgraph {input.host_ref} > "
        "test"
        # "{output.bedgraph}"


rule conservation:
    input:
        sno_host = join(config['path']['sno_host_data'],
                        config['file']['sno_host_data']),
        host_ref = join(config['path']['sno_host_data'],
                        config['file']['host_ref']),
        bedgraph = join(config['path']['sno_host_data'],
                        config['file']['host_bedgraph']),
    output:
        'cons.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/conservation.py"
