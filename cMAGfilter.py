#!/usr/bin/env python3

     #   #   #    #####  ### # #  #   ###  #   #
     ## ##  # #  #       #     #  #  #   # # ## 
 ### # # # ##### # #### #### # # ### ##### ## 
##   #   # #   # #    #  #   # #  #  #     # 
 ### #   ##     # ####   #   # #  ##  ###  #

# Input
# (a) Circular contig
#    : Circular contig FASTA file
# (b) Conspecific genomes directory
#    : Directory containing conspecific genome files in FASTA format (# of conspecific genomes >= 5) 
# (c) Output directory
#    : Directory where the outputs will be stored (The output directory must not already exist)

# Output
# (a) all_by_all_conspecific_genomes 
#    : All by all alignment results of conspecific genomes
# (b) conspecific_genomes.contig_report.tsv 
#    : Report of contig processing from nucmer results of conspecific genomes
# (c) core_contig_list.txt 
#    : List of core contigs and corresponding conspecific genomes 
#      (Core_contig [TAB] Source_genome)
# (d) core_contigs.fna 
#    : FASTA file of core contigs
# (e) [circular_name]_core_contigs 
#    : Nucmer results of aligning core contigs to the circular contig
# (f) [circular_name]_core_contigs_alignment.contig_stat.tsv 
#    : Stats from the point of view of the core contigs of the result of aligning the core contigs to the circular contig
#      (Core_contig [TAB] Alignment_length [TAB] Contig_length [TAB] Coverage)
# (g) [circular_name]_core_contigs_alignment.summary.tsv
#    : Stats/Summary from the point of view of the circular contig of the result of aligning the core contigs to the circular contig
#      (Circular_contig [TAB] Cir_Aln_Len [TAB] Cir_Len [TAB] Cir_Aln_Len/Cir_Len [TAB] Total_Aln_Len 
#       [TAB] Total_Ctg_Len [TAB] Total_Aln_Len/Total_Ctg_Len [TAB] Aln_Ctg_Cnt [TAB] Ctg_Cnt [TAB] Aln_Ctg_Cnt/Ctg_Cnt=Core_contig_retrieval_rate)
# (h) log

# Dependency
# (a) mummer4 (nucmer, delta-filter, and show-coords programs will be used)

# Steps
# (a) Identify core contigs shared by most conspecific genomes
# (b) Align core contigs to the query circular contig and calculate the retreival rate of the core contig

from cmagfilter.cmagfilter_run import run

if __name__ == "__main__":
    run()
