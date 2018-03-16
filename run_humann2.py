from __future__ import print_function

import os
import time
import glob
import shutil
import argparse
import subprocess


def join_humann2(outdir):
    genefamilies_out = os.path.join(outdir, 'humann2_genefamilies.tsv')
    cmd = 'humann2_join_tables ' \
          '--input {} ' \
          '--output {} ' \
          '--file_name genefamilies_relab'.format(outdir, genefamilies_out)
    p = subprocess.Popen(cmd, shell=True)
    p.wait()

    pathcoverage_out = os.path.join(outdir, 'humann2_pathcoverage.tsv')
    cmd = 'humann2_join_tables ' \
          '--input {} ' \
          '--output {} ' \
          '--file_name pathcoverage'.format(outdir, pathcoverage_out)
    p = subprocess.Popen(cmd, shell=True)
    p.wait()

    pathabundance_out = os.path.join(outdir, 'humann2_pathabundance.tsv')
    cmd = 'humann2_join_tables ' \
          '--input {} ' \
          '--output {} ' \
          '--file_name pathabundance_relab'.format(outdir, pathabundance_out)
    p = subprocess.Popen(cmd, shell=True)
    p.wait()

    output_files = glob.glob(os.path.join(outdir, 'humann2*.tsv'))
    return output_files


def normalize_humann2(genefamilies, pathabundances):
    """
    Deciding between copies per million or relative abundance schemes is dependent on the statistical tests
    you're using downstream. I'll default to relative abundance here.
    :return:
    """
    normalized_genefamilies = []
    for genefamily in genefamilies:
        output_filename = genefamily.replace('genefamilies.tsv', 'genefamilies_relab.tsv')
        cmd = 'humann2_renorm_table --input {} --output {} --units relab'.format(genefamily, output_filename)
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        normalized_genefamilies.append(output_filename)

    normalized_pathabundances = []
    for pathabundance in pathabundances:
        output_filename = pathabundance.replace('pathabundance.tsv', 'pathabundance_relab.tsv')
        cmd = 'humann2_renorm_table --input {} --output {} --units relab'.format(pathabundance, output_filename)
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        normalized_pathabundances.append(output_filename)

    return normalized_genefamilies, normalized_pathabundances


def run_humann2(fastq_file, metaphlan_dir):
    print('\nRunning MetaPhlAn2...')

    outdir = os.path.join(os.path.dirname(fastq_file), 'humann2_output')

    try:
        os.mkdir(outdir)
    except OSError:
        shutil.rmtree(outdir)
        os.mkdir(outdir)

    cmd = 'humann2 --input {} ' \
          '--output {} ' \
          '--memory-use maximum ' \
          '--threads 10 ' \
          '--metaphlan {} ' \
          '--translated-query-coverage-threshold 80.0 ' \
          '--identity-threshold 40.0'.format(fastq_file, outdir, metaphlan_dir)  # Run MetaPhlAn2 2.6.0
    print(cmd)

    p = subprocess.Popen(cmd, shell=True)
    p.wait()

    genefamilies = glob.glob(os.path.join(outdir, '*genefamilies.tsv'))
    pathabundances = glob.glob(os.path.join(outdir, '*pathabundance.tsv'))
    pathcoverage = glob.glob(os.path.join(outdir, '*pathcoverage.tsv'))

    return genefamilies, pathabundances, pathcoverage, outdir


def run_bbduk(fastq_file, min_len=50):
    outname = fastq_file.replace('.fastq', '_filtered.fastq')
    cmd = "bbduk.sh -Xmx1g in={in1} out={out1} qtrim=w " \
                "trimq=10 k=25 minlength={ml} ref=adapters" \
        .format(in1=fastq_file,
                out1=outname,
                ml=min_len)
    p = subprocess.Popen(cmd, shell=True)
    p.wait()
    return outname


class Humann2(object):
    def __init__(self, args):
        print('\033[92m' + '\033[1m' + '\nHUMAnN2' + '\033[0m')

        # Arguments
        self.args = args
        self.fastq_filenames = args.fastq_filenames
        self.metaphlan_dir = args.metaphlan_dir

        # Metadata
        self.num_reads = len(self.fastq_filenames)
        self.workdir = os.path.dirname(self.fastq_filenames[0])

        print('INPUT FILE(S): {}'.format(self.fastq_filenames))
        print('METAPHLAN DIRECTORY: {}'.format(self.metaphlan_dir))
        print('WORK DIRECTORY: {}'.format(self.workdir))

        if self.num_reads > 1:
            print('Paired end reads not yet supported')
            quit()
        else:
            self.fastq_r1 = self.fastq_filenames[0]

        # Step 1: Quality filter sequences
        filtered_fastq = run_bbduk(fastq_file=self.fastq_r1)

        # Step 2: Create gene families, pathway abundance, and pathway coverage output files
        genefamilies, pathabundances, pathcoverage, outdir = run_humann2(filtered_fastq, self.metaphlan_dir)

        # Step 3: Normalize output files
        genefamilies_relab, pathabundances_relab = normalize_humann2(genefamilies=genefamilies,
                                                                     pathabundances=pathabundances)
        print('GENEFAMILIES_RELAB: {}'.format(str(genefamilies_relab)))
        print('PATHABUNDANCES_RELAB: {}'.format(str(pathabundances_relab)))

        # Step 4: Join the output files
        output_files = join_humann2(outdir=outdir)
        print('FINAL OUTPUT FILES: {}'.format(str(output_files)))


if __name__ == '__main__':
    start = time.time()

    parser = argparse.ArgumentParser(description='HUMAnN2 wrapper')
    parser.add_argument('-f', '--fastq_filenames',
                        required=True,
                        help='Path to the FastQ file(s) you would like to perform HUMAnN2 analysis on. ',
                        nargs='*')
    parser.add_argument('-m', '--metaphlan_dir',
                        required=True,
                        help='Path to a folder containing MetaPhlAn2 (v2.6.0)')
    arguments = parser.parse_args()

    x = Humann2(arguments)

    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)

    print('\033[92m' + '\033[1m' + '\nFinished HUMAnN2 functions in %d:%02d:%02d ' % (h, m, s) + '\033[0m')
