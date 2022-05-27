from cmagfilter.get_argument import get_argument
from cmagfilter.integrity import integrity
from cmagfilter.data_process import all_by_all_alignment
from cmagfilter.data_process import define_core_contig
from cmagfilter.data_process import align_back

def run():
    # Get arguments
    (circular_contig, conspecific_genomes_dir, output_dir, 
        extension, contig_existence_identity_threshold,
        contig_existence_coverage_threshold, core_contig_length_threshold,
        core_contig_coverage_threshold, threads, nucmer_path, force_writing) = get_argument()

    # Program integrity check
    (nucmer_executable, delta_filter_executable, show_coords_executable) = integrity(nucmer_path, circular_contig, conspecific_genomes_dir,
                                                                                        extension, force_writing, output_dir)

    #EXECUTE!
    global_log = open(output_dir + "/log.txt", "w")
       
    contig2length, contig2genome_hit, total_genomes = all_by_all_alignment(conspecific_genomes_dir, contig_existence_identity_threshold,
                                                                                contig_existence_coverage_threshold, extension, nucmer_executable,
                                                                                delta_filter_executable, show_coords_executable, threads, global_log, output_dir)

    define_core_contig(contig2length, contig2genome_hit, total_genomes,
        extension, core_contig_length_threshold, core_contig_coverage_threshold,
        conspecific_genomes_dir, global_log, output_dir)

    align_back(circular_contig, extension, nucmer_executable,
        show_coords_executable, threads, global_log, output_dir)

    global_log.close()