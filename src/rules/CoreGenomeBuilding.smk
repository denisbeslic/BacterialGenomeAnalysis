rule roary:
    input:
        expand(os.path.join(config['Results'],"Prokka/{sample}/{sample}.gff"), sample = list(samples.index)),
        multiqc = os.path.join(config['Results'],"multiqc/multiqc.html")
    params:
        krakendb = config['kraken'],
        outputdir = os.path.join(config['Results'],"Roary")
    output:
        a = os.path.join(config['Results'],"Roary.done")
    log:
        "logs/Roary/test.log"
    threads: 8
    conda:
        "../envs/roary.yaml"
    shell:
        "(roary -f {params.outputdir} -p {threads} -e -n -v {input} -qc -k {params.krakendb} && touch {output.a}) 2> {log}"
