rule fastQC1:
    input:
        trimmed1=os.path.join(config['Results'],"trimmed/{sample}_1.fastq")
    output:
        html=os.path.join(config['Results'],"fastQC/after_trim/{sample}_1.html"),
        zip=os.path.join(config['Results'],"fastQC/after_trim/{sample}_1.zip")
    threads: 4
    wrapper:
        "0.35.0/bio/fastqc"


rule fastQC2:
    input:
        trimmed2=os.path.join(config['Results'],"trimmed/{sample}_2.fastq"),
    output:
        html=os.path.join(config['Results'],"fastQC/after_trim/{sample}_2.html"),
        zip=os.path.join(config['Results'],"fastQC/after_trim/{sample}_2.zip")
    threads: 4
    wrapper:
        "0.35.0/bio/fastqc"


rule QUAST:
    input:
        expand(os.path.join(config['Results'],"Prokka/{sample}/{sample}.gff"),  sample =list(samples.index)),
        contigs = os.path.join(config['Results'],"SPAdes/{sample}/contigs.fasta"),
        refAnnotated = os.path.join(config['Results'],"Prokka/{sample}/{sample}.gff")
    params:
        outdir = os.path.join(config['Results'],"QUAST/{sample}"),
        ref = config['reference']
    output:
        os.path.join(config['Results'],"QUAST/{sample}/report.txt")
    log:
        "logs/QUAST/{sample}.log"
    threads: 8
    conda:
        "../envs/fastQC.yaml"
    shell:
        "(quast.py --threads {threads} -o {params.outdir} {input.contigs} -r {params.ref} --features CDS:{input.refAnnotated}) &> {log}"


rule multiQC:
    input:
        expand(os.path.join(config['Results'],"fastQC/after_trim/{sample}_1.html"), sample =list(samples.index)),
        expand(os.path.join(config['Results'],"fastQC/after_trim/{sample}_2.html"), sample =list(samples.index)),
        expand(os.path.join(config['Results'],"Prokka/{sample}/{sample}.gff"), sample =list(samples.index)),
        expand(os.path.join(config['Results'],"QUAST/{sample}/report.txt"), sample =list(samples.index))
    output:
        os.path.join(config['Results'],"multiqc/multiqc.html")
    wrapper:
        "0.35.0/bio/multiqc"
