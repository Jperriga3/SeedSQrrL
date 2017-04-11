import subprocess
import glob
import gzip
import os
import re
import shutil
from Bio import SeqIO
from collections import defaultdict

_read_min_length = 100

_manifest_fmt = """
project = {sample}
job = genome,denovo,accurate
parameters = -GE:not=0 -AS:nop=0 -NW=mrnl=0 SOLEXA_SETTINGS -CO:msr=no

readgroup = {sample}_raw
data = {forward} {reverse}
technology = solexa
strain = {sample}
"""


def mitobim_notebook(seed_file):
    sample = os.path.splitext(os.path.basename(seed_file))[0]
    out_dirname = sample
    # If pattern matches SAMPLEID_GENE, extract sample name only
    if re.match('[A-Z]\d+_[A-Z]{2}_[A-Z]\d+$', sample) is None:
        sample = sample.rsplit('_', 1)[0]

    try:
        sample_dir = glob.glob('RawReads*/*{}'.format(sample))[0]

    except IndexError:
        print('No sample dir found for {}'.format(sample))
        return
    try:
        sample_dir = glob.glob('{}/redo/*{}'.format(sample_dir, sample))[0]
        print('Found redo dir for {}, using it'.format(sample))
    except IndexError:
        pass

    reads = glob.glob('{}/{}_*.gz'.format(sample_dir, sample))
    try:
        forward = [read for read in reads if 'R1' in read][0]
        reverse = [read for read in reads if 'R2' in read][0]
    except IndexError:
        print('Forward or reverse reads not found for {}'.format(sample))
        return

    try:
        files_to_link = {'forward': os.path.abspath(forward),
                         'reverse': os.path.abspath(reverse),
                         'seed': os.path.abspath(seed_file)}
        make_dirs_and_links(out_dirname, files_to_link=files_to_link)
    except OSError:
        print('Output dirs already exist for {}, skipping it. Delete output dirs to force re-run'
              .format(seed_file))
        return

    mira_trim(out_dirname, forward, reverse)

    extract_good_pairs_and_singletons(to_process='./{}/2-read-trimming/mixed-trimmed.fastq'.format(out_dirname),
                                      out_dir='./{}/2-read-trimming/'.format(out_dirname))

    flash_merge(out_dirname)

    run_mitobim(out_dirname)

    rename_contigs('./{samp}/{samp}_mt_candidate.fasta'.format(samp=out_dirname),
                   seed_file)


def make_dirs_and_links(sample, files_to_link):
    os.mkdir(sample)
    os.chdir(sample)

    os.mkdir('1-raw')
    os.chdir('1-raw')
    os.symlink(files_to_link['forward'], 'pe_raw_1.fastq.gz')
    os.symlink(files_to_link['reverse'], 'pe_raw_2.fastq.gz')
    os.symlink(files_to_link['seed'], 'seed.fasta')
    os.chdir('..')

    os.mkdir('2-read-trimming')
    os.chdir('2-read-trimming')
    os.symlink('pe-trimmed-1.fastq.gz', '{}_trimmed-minlength-{}-pe-1.fastq.gz'
               .format(sample, _read_min_length))
    os.symlink('pe-trimmed-2.fastq.gz', '{}_trimmed-minlength-{}-pe-2.fastq.gz'
               .format(sample, _read_min_length))
    os.symlink('se-trimmed.fastq.gz', '{}_trimmed-minlength-{}-se.fastq.gz'
               .format(sample, _read_min_length))
    os.chdir('..')

    os.mkdir('3-read-merging')
    os.chdir('3-read-merging')
    os.symlink('out.extendedFrags.fastq.gz', '{}_trimmed-minlength-{}-merged.fastq.gz'
               .format(sample, _read_min_length))
    os.chdir('..')

    os.mkdir('4-MITObim')
    os.symlink('4-MITObim/mt-candidate.fasta', '{}_mt_candidate.fasta'.format(sample))

    os.chdir('..')


def mira_trim(sample, forward, reverse):
    forward = os.path.abspath(forward)
    reverse = os.path.abspath(reverse)
    os.chdir(os.path.join(sample, '2-read-trimming'))

    # Write the manifest file
    with open('manifest.conf', 'w') as fh:
        fh.write(_manifest_fmt.format(sample=sample, forward=forward, reverse=reverse))

    with open('trim.log', 'w') as err_fh:
        # Run mira trimming
        subprocess.check_call(['mira', 'manifest.conf'], stderr=err_fh)
        # Extract trimmed reads
        err_fh.write('\n\n################\n\nextracting trimmed reads\n\n')
        subprocess.check_call(['miraconvert', '-f', 'maf', '-t', 'fastq', '-C', '-X',
                               str(_read_min_length),
                               '{sample}_assembly/{sample}_d_chkpt/readpool.maf'.format(sample=sample),
                               'mixed-trimmed'],
                              stderr=err_fh)

    os.chdir('..')
    os.chdir('..')


def flash_merge(sample):
    os.chdir(os.path.join(sample, '3-read-merging'))
    with open('flash.log', 'w') as err_fh:
        subprocess.check_call(['flash', '../2-read-trimming/pe-trimmed-1.fastq.gz',
                               '../2-read-trimming/pe-trimmed-2.fastq.gz', '-M', '150', '-z', '-d', './'],
                              stderr=err_fh)
    os.chdir('..')
    os.chdir('..')


def run_mitobim(sample):
    os.chdir(os.path.join(sample, '4-MITObim'))
    with open('MITObim.log', 'w') as err_fh:
        subprocess.check_call(['mitobim.pl', '--quick', '../1-raw/seed.fasta', '-readpool',
                               '../3-read-merging/out.extendedFrags.fastq.gz', '-sample', 'sample',
                               '-ref', 'seed', '-end', '50', '-kbait', '15', '--clean'],
                              stderr=err_fh)
    # Now symlink results
    iterations = glob.glob('*iteration*')
    iterations = human_sorted(iterations)
    last_iter = iterations[-1].replace('iteration', '')
    os.symlink(os.path.join(iterations[-1], 'temp_baitfile.fasta'), 'mt-candidate.fasta')
    os.symlink(os.path.join(iterations[-1], 'sample-readpool-it{}.fastq'.format(last_iter)),
               'mt-candidate-readpool.fastq')

    os.chdir('..')
    os.chdir('..')


def human_sorted(l):
    """ Sort the given iterable in the way that humans expect."""

    def convert(text):
        return int(text) if text.isdigit() else text

    def alphanum_key(key):
        return [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(l, key=alphanum_key)


def open_fastq(filename, mode):
    """
    Function that opens a fastq file handle
    """
    if mode not in ('r', 'w'):
        raise IOError('Mode must be r or w')

    if filename.endswith('.gz'):
        fh = gzip.open(filename, mode + 'b')
    else:
        fh = open(filename, mode)

    return fh


def extract_good_pairs_and_singletons(to_process, out_dir):
    """
    The function parses a fastq file (gzipped supported)
    and separates paired end from singleton reads
    """
    id_dict = defaultdict(list)

    to_process_fh = open_fastq(to_process, mode='r')
    pe_1_fh = open_fastq(out_dir + '/pe-trimmed-1.fastq', mode='w')
    pe_2_fh = open_fastq(out_dir + '/pe-trimmed-2.fastq', mode='w')
    se_fh = open_fastq(out_dir + '/se-trimmed.fastq', mode='w')

    for read in SeqIO.parse(to_process_fh, 'fastq'):
        read_id = read.id[:-2]
        id_dict[read_id].append(read)
        if len(id_dict[read_id]) == 2:
            # Sort by full read id (including the .X part at the end)
            id_dict[read_id] = sorted(id_dict[read_id], key=lambda x: x.id)
            SeqIO.write(id_dict[read_id][0], pe_1_fh, 'fastq')
            SeqIO.write(id_dict[read_id][1], pe_2_fh, 'fastq')
            del id_dict[read_id]

    for read_id in id_dict.keys():
        SeqIO.write(id_dict[read_id][0], se_fh, 'fastq')

    to_process_fh.close()
    pe_1_fh.close()
    pe_2_fh.close()
    se_fh.close()
    # Shell out to parallelize gzip
    subprocess.check_call('echo {dir}/pe-trimmed-1.fastq '
                          '{dir}/pe-trimmed-2.fastq {dir}/se-trimmed.fastq '
                          '| xargs -n 1 -P 3 gzip -9'.format(dir=out_dir), shell=True)

    del id_dict


def rename_contigs(contigsfile, seedfile):
    with open(seedfile) as fh_seeds:
        seed_names = [line for line in fh_seeds.readlines() if line.startswith('>')]
    with open(contigsfile, mode='r') as contig_fh:
        contig_lines = contig_fh.readlines()
    out_lines = []
    for seed_name, contig_name, contig in zip(seed_names, contig_lines[0::2],
                                              contig_lines[1::2]):
        print "SEED_NAME IS " + str(seed_name)
        print "CONTIG_NAME IS" + str(contig_name)
        assert (seed_name.split()[0] in contig_name)
        out_lines += [seed_name.replace('\n', '_mitobim_assembly\n'), contig]
    with open(contigsfile, mode='w') as contig_fh:
        contig_fh.writelines(out_lines)


def split_seeds(seed_file):
    with open(seed_file) as fh:
        lines = fh.readlines()
    for seed_name, seed_seq in zip(lines[0::2], lines[1::2]):
        gene_name = seed_name.split('_')[0].replace('>', '')
        new_file = seed_file.replace('.seeds', gene_name + '.seeds')
        with open(new_file, 'w') as fh:
            fh.writelines([seed_name, seed_seq])
    shutil.move(seed_file, seed_file + '.bak')


if __name__ == '__main__':
    all_seeds = glob.glob('seeds/*.seeds')
    for seed in all_seeds:
        mitobim_notebook(seed)
