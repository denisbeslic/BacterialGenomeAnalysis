rule spades:
    input:
        fastqc=expand(os.path.join(config['Results'],"fastQC/after_trim/{sample}_1.html"), sample =list(samples.index)),
        fastqc2=expand(os.path.join(config['Results'],"fastQC/after_trim/{sample}_2.html"),  sample =list(samples.index)),
        r1=os.path.join(config['Results'],"trimmed/{sample}_1.fastq"),
        r2=os.path.join(config['Results'],"trimmed/{sample}_2.fastq")
    params:
        os.path.join(config['Results'],"SPAdes/{sample}/")
    output:
        os.path.join(config['Results'],"SPAdes/{sample}/contigs.fasta")
    log:
        "logs/SPAdes/{sample}.log"
    threads: 2
    conda:
        "../envs/spades.yaml"
    shell:
        "(spades.py --careful --cov-cutoff auto -1 {input.r1} -2 {input.r2} -k 55,77 -o {params}) &> {log}"
