from __future__ import print_function

import os
import time
import glob
import shutil
import argparse
import subprocess

class Humann2(object):

    def run_humann2(self, fastq_file):
        print('\nRunning MetaPhlAn2...')

        # Run metaphlan 2.6.0 (required version)
        metaphlan_dir = '/home/dussaultf/Applications/metaphlan2_old'

        outdir = os.path.join(os.path.dirname(fastq_file), 'humann2_output')

        try:
            os.mkdir(outdir)
        except:
            shutil.rmtree(outdir)
            os.mkdir(outdir)

        cmd = 'humann2 --input {} ' \
              '--output {} ' \
              '--memory-use maximum ' \
              '--threads 10 ' \
              '--metaphlan {}'.format(fastq_file, outdir, metaphlan_dir)
        print(cmd)

        p = subprocess.Popen(cmd, shell=True)
        p.wait()

        genefamilies = glob.glob(os.path.join(outdir,'*genefamilies.tsv'))
        pathabundance = glob.glob(os.path.join(outdir,'*pathabundance.tsv'))
        pathcoverage = glob.glob(os.path.join(outdir,'*pathcoverage.tsv'))

        return genefamilies, pathabundance, pathcoverage, outdir


    def normalize_humann2(self, genefamilies):
        """
        Deciding between copies per million or relative abundance schemes is dependent on the statistical tests
        you're using downstream. I'll default to relative abundance here.
        :return:
        """
        normalized_genefamilies = []
        for genefamily in genefamilies:
            output_filename = genefamily.replace('genefamilies.tsv', 'genefamilies_relab.tsv')
            cmd = 'humann2_renorm_table --input {} --output {} --units relab'.format(genefamily, output_filename)
            print(cmd)
            p = subprocess.Popen(cmd, shell=True)
            p.wait()
            normalized_genefamilies.append(output_filename)
        return normalized_genefamilies

    def join_humann2(self, outdir):
        genefamilies_out = os.path.join(self.workdir, 'humann2_genefamilies.tsv')
        cmd = 'humann2_join_tables ' \
              '--input {} ' \
              '--output {} ' \
              '--file_name genefamilies_relab'.format(outdir, genefamilies_out)
        p = subprocess.Popen(cmd, shell=True)
        p.wait()

        pathcoverage_out = os.path.join(self.workdir, 'humann2_pathcoverage.tsv')
        cmd = 'humann2_join_tables ' \
              '--input {} ' \
              '--output {} ' \
              '--file_name pathcoverage'.format(outdir, pathcoverage_out)
        p = subprocess.Popen(cmd, shell=True)
        p.wait()

        pathabundance_out = os.path.join(self.workdir, 'humann2_pathabundance.tsv')
        cmd = 'humann2_join_tables ' \
              '--input {} ' \
              '--output {} ' \
              '--file_name pathabundance_relab'.format(outdir, pathabundance_out)
        p = subprocess.Popen(cmd, shell=True)
        p.wait()

        output_files = glob.glob(os.path.join(outdir,'humann2*.tsv'))
        return output_files

    def __init__(self, args):
        print('\033[92m' + '\033[1m' + '\nHumann2' + '\033[0m')

        # Arguments
        self.args = args
        self.fastq_filenames = args.fastq_filenames

        # Metadata
        self.num_reads = len(self.fastq_filenames)
        self.workdir = os.path.dirname(self.fastq_filenames[0])

        print('INPUT FILE(S): {}'.format(self.fastq_filenames))

        if self.num_reads > 1:
            print('Paired end reads not yet supported')
            quit()
        else:
            self.fastq_r1 = self.fastq_filenames[0]

        # Step 1: Create gene families, pathway abundance, and pathway coverage output files
        genefamilies, pathabundance, pathcoverage, outdir = self.run_humann2(self.fastq_r1)

        # Step 2: Normalize output files
        genefamilies_relab = self.normalize_humann2(genefamilies=genefamilies)
        print('GENEFAMILIES_RELAB: {}'.format(str(genefamilies_relab)))

        # Step 3: Join the output files
        output_files = self.join_humann2(outdir=outdir)
        print('FINAL OUTPUT FILES: {}'.format(str(output_files)))


if __name__ == '__main__':
    start = time.time()

    parser = argparse.ArgumentParser(description='Humann2 wrapper')
    parser.add_argument('fastq_filenames',
                        help='Path to the FastQ file(s) you would like to perform Humann2 analysis on. ',
                        nargs='*')
    arguments = parser.parse_args()

    x = Humann2(arguments)

    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)

    print('\033[92m' + '\033[1m' + '\nFinished Humann2 functions in %d:%02d:%02d ' % (h, m, s) + '\033[0m')