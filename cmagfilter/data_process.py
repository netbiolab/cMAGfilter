import os
from itertools import combinations
from .logging import print2

def to_genome_name(genome_path, extension):
    ext = "." + extension
    genome_file_name = genome_path.split("/")[-1]
    ext_index = genome_file_name.index(ext)
    return genome_file_name[:ext_index]

def print2(string, f_log):
    print(string)
    print(string,file=f_log,flush=True)

def prolong(previous_list, new_start, new_end):
    if new_start > new_end:
        exit()
       
    result_list = list()
    affected_list_idx = list()
    for i in range(len(previous_list)):
        start, end = previous_list[i]
        if (end >= new_start) and (start <= new_end):
            affected_list_idx.append(i)
            
    if len(affected_list_idx) >= 1:  
        affected_min_idx = affected_list_idx[0]
        affected_max_idx = affected_list_idx[-1]
        if  new_start <= previous_list[affected_min_idx][0]:
            concat_start = new_start
        else:
            concat_start = previous_list[affected_min_idx][0]
        
        if new_end >= previous_list[affected_max_idx][1]:
            concat_end = new_end
        else:
            concat_end = previous_list[affected_max_idx][1]
            
        added = False
        for i in range(len(previous_list)):
            if (i < affected_min_idx) or (i > affected_max_idx) :
                result_list.append(previous_list[i])
            else:
                if not added:
                    added = True
                    result_list.append([concat_start, concat_end])               
    else:
        added = False
        prev_end = -1
        for i in range(len(previous_list)):
            start, end = previous_list[i]
            if (new_start > prev_end) and (new_start < start):
                added = True
                result_list.append([new_start, new_end])
                result_list.append(previous_list[i])
            else:
                result_list.append(previous_list[i])
            prev_end = end
        
        if not added:
            result_list.append([new_start, new_end])
    
    return result_list

def alignment_length(lol):
    aln_len = 0
    for start, end in lol:
        unit = end - start + 1
        aln_len += unit

    return aln_len


def all_by_all_alignment(conspecific_genomes_dir, contig_existence_identity_threshold, contig_existence_coverage_threshold, extension, nucmer_executable, delta_filter_executable, show_coords_executable, threads, global_log, output_dir):
    # Step (a) : Identify core contigs shared by most conspecific genomes
    print2('--------------------------------------------------------------------------------------------------------------------', global_log)
    print2('Step (a) : Identify core contigs shared by most conspecific genomes', global_log)
    print2('', global_log)
    files = os.listdir(conspecific_genomes_dir)
    conspecific_genomes = [conspecific_genomes_dir+"/"+fil for fil in files if fil.endswith(extension)]
    conspecific_genomes.sort()
    conspecific_genomes_all_pairs = list(combinations(conspecific_genomes,2))
    conspecific_genomes_all_pairs.sort()
    pair_count = len(conspecific_genomes_all_pairs)
    print2("Total pair count : {}".format(pair_count), global_log)
    print2("Running nucmer for all by all alignment of conspecific genomes...", global_log)

    nucmer_completed = 0
    for pair in conspecific_genomes_all_pairs:
        genome1 = pair[0]
        genome2 = pair[1]

        genome1_name = to_genome_name(genome1, extension)
        genome2_name = to_genome_name(genome2, extension)
        
        os.makedirs(output_dir+"/all_by_all_alignment_results/"+genome1_name,exist_ok=True)
        prefix = output_dir+"/all_by_all_alignment_results/"+genome1_name+'/'+genome1_name+'___'+genome2_name
        delta = prefix+'.delta'
        delta_filter = prefix+'.delta.filter'
        coord = prefix+'.coords'

        nucmer_command = nucmer_executable+' --mum -t '+str(threads)+' -p '+prefix+' '+genome1+' '+genome2
        nucmer_response = os.system(nucmer_command)
        assert (nucmer_response == 0), "nucmer failed (all by all alignment step) - '"+nucmer_command+"'"

        delta_filter_command = delta_filter_executable+' -r -q '+delta+' > '+delta_filter
        delta_filter_response = os.system(delta_filter_command)
        assert (delta_filter_response == 0), "delta-filter failed (all by all alignment step) = '"+delta_filter_command+"'"

        show_coords_command = show_coords_executable+' '+delta_filter+' > '+coord
        show_coords_response = os.system(show_coords_command)
        assert (show_coords_response == 0), "show-coords failed (all by all alignment step) = '"+show_coords_command+"'"

        nucmer_completed += 1
        if nucmer_completed % 10 == 0 or nucmer_completed == pair_count:
            print2("Completed : {} / {} pairs".format(nucmer_completed,pair_count), global_log)



    print2('', global_log)
    print2("Computing length of contigs...", global_log)
    print2('', global_log)

    contig2length = dict()
    for genome in conspecific_genomes:
        genome_name = to_genome_name(genome, extension)
        f_genome = open(genome)
        for line in f_genome:
            line = line.strip()
            if line.startswith('>'):
                seq_name = line[1:].split()[0]
                #to avoid duplicated contig name, concat contig name with source genome
                contig_id = seq_name+'\t'+genome_name
                contig2length[contig_id] = 0
            else:
                contig2length[contig_id] += len(line)
        f_genome.close()



    contig2genome_hit = dict()
    total_genomes = set()
    print2("Processing nucmer results of conspecific genomes...", global_log)

    processing_completed = 0
    for pair in conspecific_genomes_all_pairs:
        genome1 = pair[0]
        genome2 = pair[1]

        genome1_name = to_genome_name(genome1, extension)
        genome2_name = to_genome_name(genome2, extension)

        total_genomes.add(genome1_name)
        total_genomes.add(genome2_name)

        coord_fname = output_dir+"/all_by_all_alignment_results/"+genome1_name+'/'+genome1_name+'___'+genome2_name+".coords"

        start = False
        contig2hit_len = dict()
        with open(coord_fname) as inf:
            for line in inf:
                line = line.strip()
                if not start:
                    if line.startswith('==='):
                        start = True
                        continue
                else:
                    elements = line.split()
                    len1, len2, identity, contig1, contig2 = int(elements[6]),int(elements[7]),float(elements[9]),elements[11],elements[12]

                    #to avoid duplicated contig name, concat contig name with source genome
                    contig1 = contig1+'\t'+genome1_name
                    contig2 = contig2+'\t'+genome2_name

                    #always hit to self
                    try:
                        contig2genome_hit[contig1].add(genome1_name)
                    except KeyError:
                        contig2genome_hit[contig1] = {genome1_name}
                    try:
                        contig2genome_hit[contig2].add(genome2_name)
                    except KeyError:
                        contig2genome_hit[contig2] = {genome2_name}

                    #if identity > threshold --> add aligned length
                    if identity > contig_existence_identity_threshold:
                        try:
                            contig2hit_len[contig1] += len1
                        except KeyError:
                            contig2hit_len[contig1] = len1
                        try:
                            contig2hit_len[contig2] += len2
                        except KeyError:
                            contig2hit_len[contig2] = len2

                        if contig2hit_len[contig1] > (contig2length[contig1] * contig_existence_coverage_threshold / 100):
                            contig2genome_hit[contig1].add(genome2_name)
                        if contig2hit_len[contig2] > (contig2length[contig2] * contig_existence_coverage_threshold / 100):
                            contig2genome_hit[contig2].add(genome1_name)
        
        processing_completed += 1
        if processing_completed % 10 == 0 or processing_completed == pair_count:
            print2("Completed : {} / {} pairs".format(processing_completed,pair_count), global_log)
    
    return contig2length, contig2genome_hit, total_genomes


def define_core_contig(contig2length, contig2genome_hit, total_genomes, extension, core_contig_length_threshold, core_contig_coverage_threshold, conspecific_genomes_dir, global_log, output_dir):
    print2('', global_log)
    print2("Identifying core contigs...", global_log)
    detailed_report_fp = open(output_dir+"/conspecific_genomes.contig_report.tsv","w")

    total_genomes = list(total_genomes)
    total_genomes.sort()
    genome_count = len(total_genomes)
    over_coverage_contigs = list()
    contig_count = len(contig2genome_hit)
    genome2core_contigs = dict()

    print ("#CONTIGS\tSOURCE_GENOME\tLENGTH\tCOVERAGE\tANNOT\t"+"\t".join(total_genomes), file=detailed_report_fp)
    for contig in contig2genome_hit:
        genome_hit = contig2genome_hit[contig]
        hit_status = list()
        hit_count = 0
        for genome in total_genomes:
            if genome in genome_hit:
                hit_count += 1
                hit_status.append('O')
            else:
                hit_status.append('X')

        genome_coverage = hit_count / genome_count * 100
        annot = "not_shared"
        if genome_coverage >= core_contig_coverage_threshold:
            if contig2length[contig] >= core_contig_length_threshold:
                over_coverage_contigs.append(contig)
                c,g = contig.split('\t')
                try:
                    genome2core_contigs[g].add(c)
                except KeyError:
                    genome2core_contigs[g] = {c}
                annot = "core"
            else:
                annot = "shared_but_short"

        print(contig+'\t'+str(contig2length[contig])+'\t'+str(genome_coverage)+'\t'+annot+'\t'+'\t'.join(hit_status), file=detailed_report_fp)

    pass_rate = len(over_coverage_contigs) / contig_count
    print ('#Among '+str(contig_count)+' contigs, '+str(len(over_coverage_contigs))+' contigs are core contigs, '+str(pass_rate*100)+'%', file=detailed_report_fp)
    print2('Among '+str(contig_count)+' contigs, '+str(len(over_coverage_contigs))+' contigs are core contigs, '+str(pass_rate*100)+'%', global_log)
    detailed_report_fp.close()
    print2("Write complete : conspecific_genomes.contig_report.tsv", global_log)


    printed_contigs = set()
    contig_seq_out = open(output_dir+"/core_contigs.fna","w")

    for neighbor in genome2core_contigs:
        core_contigs = genome2core_contigs[neighbor]
        genome_seq = conspecific_genomes_dir+"/"+neighbor+"."+extension
        print_flag = False
        with open(genome_seq) as fin:
            for line in fin:
                line = line.strip()
                if line.startswith('>'):
                    contig = line[1:].split()[0]
                    if contig in core_contigs:
                        print_flag = True
                        print (line+'___'+neighbor, file=contig_seq_out)
                        printed_contigs.add(contig+'\t'+neighbor)
                    else:
                        print_flag = False
                else:
                    if print_flag:
                        print(line, file=contig_seq_out)

    contig_seq_out.close()
    print2("Write complete : core_contigs.fna", global_log)
    assert len(printed_contigs) == len(over_coverage_contigs),"Only some of the core contigs have been written."


def align_back(circular_contig, extension, nucmer_executable, show_coords_executable, threads, global_log, output_dir):
    # Step (b) : Align core contigs to the query circular contig and calculate the retreival rate of the core contig
    print2('--------------------------------------------------------------------------------------------------------------------', global_log)
    print2('Step (b) : Align core contigs to the query circular contig and calculate the retreival rate of the core contig', global_log)
    print2('', global_log)
    circular_name = to_genome_name(circular_contig, extension)

    print2('Aligning core contigs to the query circular contig using nucmer...', global_log)
    os.makedirs(output_dir+"/"+circular_name+"_align_back_results", exist_ok=True)
    prefix = output_dir+"/"+circular_name+"_align_back_results/"+circular_name+"_align_back"


    nucmer_command = nucmer_executable+' --maxmatch -t '+str(threads)+' -p '+prefix+" "+circular_contig+" "+output_dir+"/core_contigs.fna"
    nucmer_response = os.system(nucmer_command)
    assert (nucmer_response == 0), "nucmer failed (align back step) - '"+nucmer_command+"'"

    show_coords_command = show_coords_executable+' '+prefix+'.delta > '+prefix+'.coords'
    show_coords_response = os.system(show_coords_command)
    assert (show_coords_response == 0), "show-coords failed (align back step) = '"+show_coords_command+"'"

    print2('Nucmer finished', global_log)
    print2('', global_log)
    print2('Calculating core contig retrieval rate...', global_log)

    circular_len = 0
    with open(circular_contig) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                continue
            else:
                circular_len += len(line)

    core_contig2length = dict()
    with open(output_dir+"/core_contigs.fna") as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_contig = line[1:]
                core_contig2length[current_contig] = 0
            else:
                core_contig2length[current_contig] += len(line)

    start = False
    alignment2contig = dict()
    with open(output_dir+"/"+circular_name+"_align_back_results/"+circular_name+"_align_back.coords") as f:
        for line in f:
            line = line.strip()
            if line.startswith('======='):
                start=True
                continue
            if start:
                s1, e1, bar, s2, e2, bar, len1, len2, bar, identity, bar, circular, contig = line.split()
                s1, e1 = sorted([int(s1), int(e1)])
                s2, e2 = sorted([int(s2), int(e2)])
                try:
                    alignment4circular = prolong(alignment4circular, s1, e1)
                except NameError:
                    alignment4circular = [[s1,e1]]

                try:
                    alignment2contig[contig] = prolong(alignment2contig[contig], s2, e2)
                except KeyError:
                    alignment2contig[contig] = [[s2,e2]]


    by_corectg_outf = open(output_dir+"/"+circular_name+"_core_contigs_alignment.core_contig_stat.tsv","w")
    print('#CORE_CONTIG\tALN_LEN\tCTG_LEN\tCOVERAGE',file=by_corectg_outf)
    summary_outf = open(output_dir+"/"+circular_name+"_core_contigs_alignment.summary.tsv","w")
    print('#CIRCULAR\tCIR_ALN_LEN\tCIR_LEN\tCIR_ALN_LEN/CIR_LEN\tTOT_ALN_LEN\tTOT_CTG_LEN\tTOT_ALN_LEN/TOT_CTG_LEN\tALN_CTG_CNT\tCTG_CNT\tALN_CTG_CNT/CTG_CNT(CORE_CONTIG_RETRIEVAL_RATE)',file=summary_outf)

    total_aln_len = 0
    total_ctg_len = 0
    aln_ctg_cnt = 0
    for contig in core_contig2length:
        try:
            aln_len = alignment_length(alignment2contig[contig])
            total_aln_len += aln_len
            aln_ctg_cnt += 1
        except KeyError:
            aln_len = 0

        ctg_len = core_contig2length[contig]
        total_ctg_len += ctg_len
        coverage = aln_len / ctg_len * 100

        print (contig+'\t'+str(aln_len)+'\t'+str(ctg_len)+'\t'+str(coverage), file=by_corectg_outf)

    circular_aln_len = alignment_length(alignment4circular)
    summary = [circular_aln_len, circular_len, circular_aln_len/circular_len, total_aln_len, total_ctg_len, total_aln_len/total_ctg_len, aln_ctg_cnt, len(core_contig2length), aln_ctg_cnt/len(core_contig2length)]

    print (circular_name+'\t'+'\t'.join(map(str, summary)), file=summary_outf)
    print2('Core contig retrieval rate : {}'.format(aln_ctg_cnt/len(core_contig2length)), global_log)
    by_corectg_outf.close()
    summary_outf.close()
    print2('Write complete : {}_core_contigs_alignment.core_contig_stat.tsv'.format(circular_name), global_log)
    print2('Write complete : {}_core_contigs_alignment.summary.tsv'.format(circular_name), global_log)
    print2('--------------------------------------------------------------------------------------------------------------------', global_log)
    print2('Finished', global_log)


