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
        cons = join(config['path']['sno_host_data'],
                    config['file']['cons'])
    output:
        alt_splice = join(config['path']['sno_host_data'],
                          config['file']['alt_splice']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/alternative_splicing_intron.py"

rule prepare_bedgraph_search:
    input:
        alt_splice = join(config['path']['sno_host_data'],
                          config['file']['alt_splice']),
    output:
        introns = join(config['path']['tmp'],
                       config['file']['introns_sno_in_host']),
    shell:
        "colTab -f {input.alt_splice} -c chr1,intron_start,intron_end "
        "| awk '{{if(NF == 3 && NR > 1){{print \"chr\" $0}}}}' "
        "| sort -k1,1 -k2,2n -k3,3 "
        "| uniq "
        "| awk 'BEGIN{{FS=\"\t\";OFS=\"\t\"}}{{print$1, $2-100, $3+100}}' "
        "> {output.introns}"

rule get_intron_bed:
    input:
        bed = join(config['path']['beds'], 'sorted_clean_{bed}.bedgraph'),
        introns = join(config['path']['tmp'],
                       config['file']['introns_sno_in_host']),
    output:
        bg = join(config['path']['tissues_cell_bg'],
                  'intron_{bed}.bedgraph'),
    conda:
        "../envs/python.yaml"
    shell:
        "set +o pipefail; cat {input.bed} | scripts/getIntronBG {input.introns} "
        "> {output.bg}"

rule reads_in_extensions:
    input:
        alt_splice = join(config['path']['sno_host_data'],
                          config['file']['alt_splice']),
        bg = expand(join(config['path']['tissues_cell_bg'], 'intron_{bed}.bedgraph'),
                    bed=config['bedgraphs'])
    output:
        bed_viz = join(config['path']['tmp'],
                       config['file']['bedgraph_viz']),
        ext_ratio = join(config['path']['sno_host_data'],
                         config['file']['ext_ratio']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/reads_in_extensions.py"
