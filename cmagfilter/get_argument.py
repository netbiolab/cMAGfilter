import argparse

def get_argument():
    parser = argparse.ArgumentParser(description='cMAGfilter : Calculate the retrieval rate of the core contig for the circular contig')
    ## Positional arguments
    parser.add_argument('circular_contig',help='Circular contig FASTA file')
    parser.add_argument('conspecific_genomes_dir',help='Directory containing conspecific genome files in FASTA format (# of conspecific genomes >= 5)')
    parser.add_argument('output_dir',help='Directory where the outputs will be stored (The output directory must not already exist)')
    ## Optional arguments
    parser.add_argument('-x','--extension',default='fna',metavar='EXTENSION',help='Extension of contig/genome FASTA files (Default : %(default)s)')
    parser.add_argument('-cei','--contig_existence_identity_threshold',default=95,type=float,metavar='CEI_THRESHOLD[%]',help='Contig is present in other conspecific genome when alignment identity > N%% (Default : %(default)s)')
    parser.add_argument('-cec','--contig_existence_coverage_threshold',default=50,type=float,metavar='CEC_THRESHOLD[%]',help='Contig is present in other conspecific genome when the N%% of contig sequences are covered by other genome (Default : %(default)s)')
    parser.add_argument('-ccl','--core_contig_length_threshold',default=5000,type=int,metavar='CCL_THRESHOLD[bp]',help='Contig longer than this threshold could be a core contig (Default : %(default)s)')
    parser.add_argument('-ccc','--core_contig_coverage_threshold',default=80,type=float,metavar='CCC_THRESHOLD[%]',help='Contig is core contig when the contig is present in N%% of other conspecific genomes (Default : %(default)s)')
    parser.add_argument('-t','--threads',default=1,type=int,metavar='THREADS',help='Maximum number of threads to use when running nucmer (Default : %(default)s)')
    parser.add_argument('-nuc','--nucmer_path',default='PATH',type=str,metavar='NUCMER_PATH',help='nucmmer4 executable path (Default : %(default)s)')
    parser.add_argument('-F','--force_writing', default=False, action="store_true", help='force_writing (Default : %(default)s)')
    args = parser.parse_args()

    return (args.circular_contig, 
                args.conspecific_genomes_dir, 
                args.output_dir, 
                args.extension, 
                args.contig_existence_identity_threshold, 
                args.contig_existence_coverage_threshold, 
                args.core_contig_length_threshold, 
                args.core_contig_coverage_threshold,
                args.threads, 
                args.nucmer_path, 
                args.force_writing)