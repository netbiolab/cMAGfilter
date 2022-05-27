import os
from . import data_process

def check_executable(fpath):
    return os.path.exists(fpath) and os.access(fpath, os.X_OK) and os.path.isfile(fpath)

def check_requirements(path, requirements):
    ref = [True] * len(requirements)
    chk = [False] * len(requirements)

    for i in range(len(requirements)):
        requirement = requirements[i]
        if path == "PATH":
            for path in os.environ['PATH'].split(os.pathsep):
                if check_executable(path+"/"+requirement):
                    chk[i] = True
                    break
        else:
            chk[i] = check_executable(path+"/"+requirement)

    return chk == ref

def integrity(nucmer_path, circular_contig, conspecific_genomes_dir, extension, force_writing, output_dir):
    requirements = ['nucmer', 'delta-filter', 'show-coords']
    requirement_fulfilled = check_requirements(nucmer_path, requirements)
    assert requirement_fulfilled, 'CANNOT find nucmer4 (nucmer, delta-filter, show-coords) package from '+nucmer_path+', Please install nucmer'
    if nucmer_path == "PATH":
        nucmer_executable = 'nucmer'
        delta_filter_executable = 'delta-filter'
        show_coords_executable = 'show-coords'
    else:
        nucmer_executable = nucmer_path +'/'+'nucmer'
        delta_filter_executable = nucmer_path +'/'+'delta-filter'
        show_coords_executable = nucmer_path +'/'+'show-coords'

    assert os.path.isfile(circular_contig),"Circular contig file does not exist"
    assert os.path.isdir(conspecific_genomes_dir),"Conspecific genomes directory does not exist"
    assert circular_contig.endswith(extension),"Circular contig does not have extension '{}'".format(extension)
    
    files = os.listdir(conspecific_genomes_dir)
    conspecific_genomes = [conspecific_genomes_dir+"/"+fil for fil in files if fil.endswith(extension)] 
    assert len(conspecific_genomes) >= 5, "The number of conspecifc genomes must be at least 5"
    
    if not force_writing:
        assert (not os.path.isdir(output_dir)), "Output directory already exists"
        os.makedirs(output_dir)
    else:
        os.makedirs(output_dir, exist_ok=True)

    print('Integrity check finished\n')
    print('Calculate the retrieval rate of the core contig for the circular contig\n')
    print('Circular contig : {}'.format(circular_contig))
    print('The number of conspecific genomes : {}'.format(len(conspecific_genomes)))
    print('')

    return nucmer_executable, delta_filter_executable, show_coords_executable