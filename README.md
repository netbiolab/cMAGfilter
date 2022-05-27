# cMAGfilter
Filtering faulty circularized genome using conspecific MAGs.

### Workflow introduction
Reconstruction of the circularized genomes is the ultimate goal of prokaryotic genome assembly. The development of accurate long-read sequencing technology enables the assembly of circularized genomes from highly complex metagenomic samples. However, prokaryotic genomes tend to have many repetitive sequences, and those often result in faulty assembly closing, thereby generating circularized genomes with significant gaps. cMAGfilter filters out the circularized metagenome-assembled genomes (cMAGs) with such gaps using their conspecific MAGs. For a given cMAG and its conspecific MAGs, it first searches core contigs, the contigs shared by most of the conspecific MAGs, from conspecific MAGs. Next, it calculates the core contig retrieval rate from the cMAG and filters out the cMAG using the information.

![](images/introductory.png)

### Requirements
cMAGfilter requires Python>=3.6 and [mummer4](https://mummer4.github.io/) package.
You can install mummer from its [tarball](https://github.com/mummer4/mummer/releases) or from [bioconda](https://bioconda.github.io/recipes/mummer4/README.html?highlight=mummer4#package-package%20&#x27;mummer4&#x27;).
Please locate the mummer4 package softwares in PATH or specify the location with -nuc parameter.

### Setting up cMAGfilter
'''
python3 setup.py install --user
'''
