rule trimmomatic_pe:
    input:
        r1=lambda wildcards: samples.at[wildcards.sample, 'fq1'] \
                if wildcards.sample in samples.index else '',
        r2=lambda wildcards: samples.at[wildcards.sample, 'fq2'] \
                if wildcards.sample in samples.index else ''
    output:
        r1=os.path.join(config['Results'],"trimmed/{sample}_1.fastq"),
        r2=os.path.join(config['Results'],"trimmed/{sample}_2.fastq"),
        # reads where trimming entirely removed the mate
        r1_unpaired=os.path.join(config['Results'],"trimmed/{sample}_1.unpaired.fastq"),
        r2_unpaired=os.path.join(config['Results'],"trimmed/{sample}_2.unpaired.fastq")
    params:
        trimmer=config['trimmomatic']
    threads:
        4
    log:
        "logs/trimmomatic/{sample}.log"
    #conda:
    #    "../envs/trimmomatic.yaml"
    wrapper:
        "0.35.1/bio/trimmomatic/pe"
