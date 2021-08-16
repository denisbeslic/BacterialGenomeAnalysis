rule raxml:
    input:
        r =  os.path.join(config['Results'],"Roary.done")
    params:
        dir=os.path.join(config['Results'],"raxml/coregenome/"),
        filename="coregenome",
        align = os.path.join(config['Results'],"Roary/core_gene_alignment.aln")
    output:
        finaloutput = os.path.join(config['Results'],"raxml/coregenome/RAxML_bestTree.coregenome")
    log:
        "logs/RAxML/coregenome.log"
    threads: 8
    conda:
        "../envs/raxml.yaml"
    shell:
        "(raxmlHPC -m GTRGAMMA -p 12345 -s {params.align} -w {params.dir} -# 20 -n {params.filename}) &> {log}"

rule roary_visualise:
    input:
        os.path.join(config['Results'],"raxml/coregenome/RAxML_bestTree.coregenome")
    params:
        os.path.join(config['Results'],"Roary/gene_presence_absence.csv")
    output:
        os.path.join(config['Results'],"roary_visualise.done")
    log:
        "logs/roary_visualise/coregenome.log"
    conda:
        "../envs/raxml.yaml"
    shell:
        "(python rules/roary_plots.py {input} {params} && touch {output}) &> {log}"

#needs Biopython
