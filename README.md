# Bacterial Genome Analysis
Automated user-friendly pipeline for whole bacterial genome analysis based on snakemake

![grafik](https://user-images.githubusercontent.com/30258228/168625956-7e6bf88e-89c0-4c88-b23d-17ba7f35a8c9.png)

## Aim of the pipeline

The aim of this project is to create an automated analysis workflow that is able to perform
comparative analysis on bacterial strains. Ultimately, the pipeline should be both easy to
maintain and update, and easy to (re-) use. It is ideally tailored to utilize powerful
multi-processor computing servers that are commonly employed by research groups and
institutions. That being said, given the emphasis on accessibility and re-usability our
pipeline should in theory be executable on most current mid- to high- end desktops and
laptops, albeit more time consuming and taxing on the respective machines.

To implement our pipeline we used Snakemake, a Python-based workflow
management tool, which automates the execution of many command line applications and
ensure reproducibility. The key idea in Snakemake is that the workflow is determined
automatically from top (the files you want) to bottom (the files you have), by applying very
general rules with the wildcards you give to Snakemake. Rules can either use shell
commands, plain Python code or external Python or R scripts to create necessary output
and inputs. Snakemake workflows can be easily executed on workstations, clusters, the
grid, and in the cloud with little to no user modification. By specifying the number of threads
to use per rule and the number of cores to use in execution, Snakemake can run multiple
jobs in parallel. Snakemake can automatically deploy required software dependencies of a
workflow using Conda .

We were able to create a workflow that takes raw whole-genome fastq reads
sequenced from a number of bacterial strains, obtained from sequencing platforms Illumina
MiSeq and/or Hiseq machines, in paired-end mode. Our pipeline then performs a de novo
assembly and annotates the assembled genomes for all individual bacterial samples in order
to subsequently obtain the core-genes across the samples. Finally - given this core genome - a phylogenetic tree is built with maximum likelihood estimation methods

## Used Tools

### Trimmomatic
Input : Raw Reads in Fastq format

Output : Reads with trimmed adaptors in Fastq format

To filter out the adapters of raw reads we used Trimmomatic, a commonly used tool for
adaptor trimming. We decided to perform quality trimming by utilizing a key module of our
chosen Genome assembler - SPAdes - which is considered very robust. For the
implementation we used the trimmomatic wrapper and defined the trimming options to only
cut off the adapter in the config file.
```
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10:8:True
```

### SPAdes
Input : Reads with trimmed adaptors in Fastq format

Output : FASTA file containing assembled contigs

After trimming the adapters (and performing a quality control check: see Quality Control
section) from the raw reads we performed a de novo genome assembly with SPAdes. De
novo genome assembly is a method for merging overlapping reads into contiguous
sequences (contigs or scaffolds) without mapping the reads to a reference genome. Most of
the de novo assembler tools use de Bruijn graphs to produce an assembly.
SPAdes is a de Bruijn graph based assembler; it is designed and intended for small
genomes and can take sequencing data from Illumina and IonTorrent platforms. It breaks
reads into fixed-size k-mers, builds a graph, cleans the graph, then finds paths through the
graph. These paths end up as contigs and scaffolds. A recent study compared different
genome assemblers on E.Coli samples. SPAdes was developed primarily for bacterial data.
What makes SPAdes an especially attractive choice is that it is essentially multiple tools
in one; and importantly, it does not sacrifice robustness or accuracy whilst maintaining a
high level of versatility. It is this integrated approach however that makes it much
simpler to incorporate in pipelines. It must be noted that it can be quite RAM intensive.
SPAdes comes in several separate modules (for Illumina reads):
- BayesHammer – read error correction tool for Illumina reads, which works well on
both single-cell and standard data sets.
- SPAdes – iterative short-read genome assembly module; values of K(-mers) are
selected automatically based on the read length and data set type.
- MismatchCorrector – a tool which improves mismatch and short indel rates inresulting contigs and scaffolds; this module uses the BWA tool.

We were using following command to call SPAdes:
```
spades.py --careful --cov-cutoff auto -1 {input.r1} -2 {input.r2} -k
33,55,77,99,121 -o {params}
```

Additionally, the entire assembly of fastq reads does not have to be restarted if you wish to
change certain parameters (such as length of k-mers to use [ --restart-from
<check_point>] ), or in the case of a crash [ --continue] , as the assembler may be
continued from the last available check-point. Check-points are made after: error correction
module completion, iteration for each specified k(-mer) value of assembly module, and
mismatch module completion

### Prokka 
Input : FASTA file containing assembled contigs

Output : annotated genome in GFF3 format, containing both sequences and annotations

Once we obtained the set of contigs in, the next step is to annotate the genomes. Genome
annotation is the process of finding and labeling interesting features on the genome. This
can include signal-peptides, ribosomal and transfer RNAs in the genome.
Prokka relies on external feature prediction tools to identify the coordinates of
genomic features within contigs. These tools are Prodigal for Coding Sequence, RNAmmer
for rRNA, Aragorn for tRNA, SignalP for Signal leader peptides and Infernal for non-coding
RNA. The ideal input for Prokka are finished sequences without gaps, but it is expected that
the typical input will be a set of scaffold sequences produced by de novo assembly with
SPAdes.
We used the following command to call Prokka in our pipeline:
```
prokka --outdir {params.out} --force --prefix {params.prefix}
--usegenus --genus {params.genus} --species {params.species} {input}
```

### Roary
Input : annotated genome in GFF3 format, containing both sequences and annotations

Output : multi-FASTA alignment of all of the core genes called core_gene_alignment.aln

Generating a pan genome describes the process of merging all genes from various strains
within a species. A Pan-genome includes genes present in all strains and genes that are
unique to each strain. This then can be used to identify the set of core genes called core
genome, which represents genes that are present only in all strains of a species. With this
information it is possible to gain a better picture of the conserved genes and in our case to
get a better understanding of the phylogeny of the different strains.
We decided to use Roary a tool that rapidly builds large-scale pan-genomes,
identifying the core and accessory genes. Roary takes as input the GFF3 files and clusters
the predicted proteins to allow for the full genomic variation of the input set to be explored.
The basic method is to convert input files to protein sequences, filter and precluster the
proteins, perform an all against all comparison using BLASTP, and cluster with MCL. Finally,
the pre-clustering results from CD-HIT are merged together with the results of MCL. Using
conserved gene neighborhood information, homologous groups containing paralogs are split into groups of true orthologs. A graph is constructed of the relationships of the clusters
based on the order of occurrence in the input sequences, allowing for the clusters to be
ordered and thus providing context for each gene.
Roary creates multiple output files, including a spreadsheet detailing the presence
and absence of each gene in each isolate, number of isolates a gene is found in, frequency
of genes per isolate, functional annotation, QC information and sorting information. A
multiple sequence file of all the nucleotide sequences for each cluster is also created and
aligned using PRANK to create a multiple sequence alignment file. For quality control Roary
allows the option to use Kraken database, which will produce a report listing all the top
species for each sample. Kraken's accuracy is comparable with Megablast, with slightly
lower sensitivity and very high precision.
```
roary -f {params.outputdir} -p {threads} -e -n -v {input.gff} -qc -k
{params.krakendb} && touch {output.a}
```

### RAxML7
Input : core_gene_alignment.aln

Output : Raxml.bestTree

RAxML (Randomized Axelerated Maximum Likelihood) is a popular program for
phylogenetic analysis of large datasets under maximum likelihood. It is an efficient maximum
likelihood tree search algorithm which returns trees so that under the assumed statistical
model the observed data is most probable.
We built the Phylogenetic trees by calling this command line:
```
raxmlHPC -m GTRGAMMA -p 12345 -s {params.align} -# 20 -n
{params.filename} -w {params.dir} -T {threads}
```

### Visualisation
GraPhLan requires several annotation files the user has to specify based on his
dataset and it is only useful when the user has a large pan genome or pan genomes from
different species, which was not the main goal of our implementation.
We included Marco Galradinis roary_plots.py script in the rules directory and can therefore
simply be called with the following command:
```
python rules/roary_plots.py {input} {params}
```
It provides 3 figures: the tree compared to a matrix with the presence and absence of
core and accessory genes, a pie chart of the breakdown of genes and the number of isolate
they are present in and a graph with the frequency of genes versus the number of genomes.
Additionally, we produced total_vs_unique_gene graphs by essentially using the code
incorporated in Roary. Finally, to improve the look of the phylogeny trees, we recommend
using the online tool iTOL12 for a more detailed display, annotation and management of the
tree.

### Quality Control

#### FastQC
The quality of the raw reads will be checked with FastQC. This tool is widely used in the field
and is considered a robust, efficient, and versatile quality control software for a varied range
of genetic data. It outputs a quality report which can be viewed to give an indication on how
good the respective reads are. We implemented FastQC by using wrappers after the
trimming of the adapters.
#### QUAST
To test the quality of our assembled genomes afterwards we used QUAST - a quality
assessment tool for evaluating and comparing genome assemblies. It is possible to assess
the quality of assemblies with or without a reference genome. The shell command used is:
quast.py --threads {threads} -o {params.outdir} {input} -r
{params.ref} --features CDS:{params.refAnnotated}
QUAST compares the de novo sample from SPAdes by taking its contigs.fasta file with the
reference genome provided in config.yaml by the user. It is possible to provide the annotated
genomes from prokka as reference_annotations. Based on the quality check with QUAST it
is possible to redo the genome assembling with SPAdes changing some parameters.
#### MultiQC
MultiQC searches a given directory for analysis logs and compiles a HTML report. It is a
general use tool, perfect for summarising the output from numerous bioinformatics tools. We
used a snakemake wrapper to utilize it in the workflow. It takes all FastQC, QUAST and even
Prokka report files and outputs the quality report as a HTML file.




