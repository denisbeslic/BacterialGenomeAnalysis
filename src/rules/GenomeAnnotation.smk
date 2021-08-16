rule prokka:
    input:
        indir = os.path.join(config['Results'],"SPAdes/{sample}/contigs.fasta")
    params:
        out = os.path.join(config['Results'],"Prokka"),
        genus = config["genus"],
        species = config["species"],
        prefix = "{sample}/{sample}"
    output:
        os.path.join(config['Results'],"Prokka/{sample}/{sample}.gff")
    log:
        "logs/Prokka/{sample}.log"
    conda:
        "../envs/prokka.yaml"
    shell:
        "(prokka --outdir {params.out} --force --prefix {params.prefix} --usegenus --genus {params.genus} --species {params.species} {input}) 2> {log}"
