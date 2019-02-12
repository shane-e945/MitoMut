import pysam

import gc
import os
import subprocess

from collections import Counter
from itertools import islice, tee

from multiprocessing import Process

class FileParser:
    def __init__(self, file_path, header, reference, write_directory, quality):
        """Class extracts deletion candidates from BAM/SAM file for analysis"""

        # Setting all basic class variables
        self.header = header
        self.file_path = file_path
        self.reference = reference
        self.quality = quality

        # Trash file for stdout of subprocess command
        self.null = open(os.devnull, 'w')

        # Write path with file_path appended for potential concurrent runs
        self.file_prefix = os.path.join(write_directory, os.path.basename(file_path))


        # Finds and then filters all potential deletion reads based on cigar
        self.find_candidates()
        self.filter_candidates('{}_bad_cigar_one.bam'.format(self.file_prefix))
        self.filter_candidates('{}_bad_cigar_two.bam'.format(self.file_prefix))

        # Combines strand one and two reads
        self.combine_candidates('{}_combined_one.fsa'.format(self.file_prefix),
                                '{}_unmapped_one.bam'.format(self.file_prefix),
                                '{}_bad_cigar_one.bam.new'.format(self.file_prefix))

        self.combine_candidates('{}_combined_two.fsa'.format(self.file_prefix),
                                '{}_unmapped_two.bam'.format(self.file_prefix),
                                '{}_bad_cigar_two.bam.new'.format(self.file_prefix))

        # Runs files through BLAT for alignment
        print('Analyzing collected reads\n')
        self.blat_file('{}_combined_one.fsa'.format(self.file_prefix),
                       '{}_one.axt'.format(self.file_prefix))
        self.blat_file('{}_combined_two.fsa'.format(self.file_prefix),
                       '{}_two.axt'.format(self.file_prefix))

        self.null.close()

    def find_candidates(self):
        """Finds all reads that potentially represent deletions"""

        # Collects improperly aligned mitochondrial reads
        print('Collecting improperly aligned reads')
        subprocess.call('samtools view -q {} -b -f 64 {} -h {} > {}_bad_cigar_two.bam'.format(
            self.quality, self.file_path, self.header, self.file_prefix), shell=True)
        subprocess.call('samtools view -q {} -b -f 128 {} -h {} > {}_bad_cigar_one.bam'.format(
            self.quality, self.file_path, self.header, self.file_prefix), shell=True)

        # Collects unmapped reads
        print('Collecting unmapped reads')
        one = Process(target=subprocess.call, args=['samtools view -b -f 68 {} > {}_unmapped_one.bam'.format(
            self.file_path, self.file_prefix)], kwargs={'shell': True})
        two = Process(target=subprocess.call, args=['samtools view -b -f 132 {} > {}_unmapped_two.bam'.format(
            self.file_path, self.file_prefix)], kwargs={'shell': True})

        one.start()
        two.start()

        one.join()
        two.join()

    @staticmethod
    def filter_candidates(bam_file):
        """Filters all potential candidates with bad cigars"""

        subprocess.call('samtools index {}'.format(bam_file),
                        shell=True)

        def check_cigar(read):
            return read.cigarstring != '{}M'.format(len(read.seq))

        bam_one = pysam.AlignmentFile(bam_file, 'rb')
        bam_two = pysam.AlignmentFile('{}.new'.format(bam_file), 'wb',
                                      template=bam_one)

        for read in bam_one.fetch():
            if check_cigar(read):
                bam_two.write(read)

        bam_one.close()
        bam_two.close()

        os.remove(bam_file)
        os.remove('{}.bai'.format(bam_file))

    def combine_candidates(self, new_file, file_one, file_two):
        """Combines two bam files together"""

        subprocess.call('samtools cat -o {}_tmp.bam {} {}'.format(self.file_prefix, file_one, file_two),
                        shell=True, stdout=self.null)
        subprocess.call('samtools fasta {}_tmp.bam > {}'.format(self.file_prefix, new_file),
                        shell=True, stdout=self.null, stderr=self.null)

        os.remove('{}_tmp.bam'.format(self.file_prefix))
        os.remove(file_one)
        os.remove(file_two)

    def blat_file(self, file, new_file):
        """Runs a file through BLAT to get an axt file"""

        subprocess.call('blat {} {} -out=axt {}'.format(self.reference, file, new_file), shell=True,
                        stdout=self.null)
        os.remove(file)


class ReadParser:
    def __init__(self, file_path, genome_length,
                 required_reads, required_size, write_directory):

        """Uses files created from FileParser to find potential deletions"""

        # Setting all basic class variables
        self.file_path = file_path
        self.genome_length = genome_length
        self.required_size = required_size
        self.required_reads = required_reads

        # Write path with file_path appended for potential concurrent runs
        self.file_prefix = os.path.join(write_directory, os.path.basename(file_path))

        self.reads_one = self.axt_candidates('{}_one.axt'.format(self.file_prefix))
        self.reads_two = self.axt_candidates('{}_two.axt'.format(self.file_prefix))

        os.remove('{}_one.axt'.format(self.file_prefix))
        os.remove('{}_two.axt'.format(self.file_prefix))

        self.reads_one = self.filter_candidates(self.reads_one)
        self.reads_two = self.filter_candidates(self.reads_two)

    @staticmethod
    def axt_candidates(axt_file):
        """Returns a list of potential deletions from a given .axt file"""

        # Opens the axt file and uses a generators to parse line by line
        with open(axt_file) as reads:
            # Creates generators of the split lines indexing over every 4th line
            lines_one, lines_two = tee((line.rstrip('\n').split() for line
                                        in islice(reads, 0, None, 4)))

            # Uses one of the generators of the file to count read id's
            id_counter = Counter(line[4] for line in lines_one)

            # Uses the second file generator to filter based upon mate presence
            cleaned_one, cleaned_two = tee(line for line in lines_two if id_counter[line[4]] == 2)

        read_dict = dict()
        for counter, line in enumerate(cleaned_one):
            # Ensures dual id comes from the same strand
            if line[4] not in read_dict:
                read_dict[line[4]] = line[7]

            elif read_dict[line[4]] != line[7]:
                read_dict[line[4]] = None

        lines = (line for line in cleaned_two if read_dict[line[4]] is not None)

        # Combines matching reads into dictionary entries for easy logging
        read_dict_two = dict()
        for read in lines:
            if read[4] not in read_dict_two:
                read_dict_two[read[4]] = [read[4], read[2], read[3], read[7]]
            else:
                read_dict_two[read[4]].extend((read[2], read[3], read[7]))

        lines = [read_dict_two[read_key] for read_key in read_dict_two if len(read_dict_two[read_key]) == 7]

        # Ensures potentially large data structures are removed immediately
        gc.collect()

        return lines

    def filter_candidates(self, reads):
        """Takes circular nature of MT into account by testing multiple indices"""

        reads = [(read[2], read[4]) for read in reads if min(self.mt_length(int(read[2]), int(read[4])))
                 >= self.required_size]

        # Returning a set of all counted deletions
        count = Counter(reads)
        return {(count[read], read[0], read[1]) for read in reads if count[read] >=
                self.required_reads}

    def mt_length(self, start, end):
        """Calculates the true length of a given mitochondrial deletion"""

        one = (start - end + self.genome_length) % self.genome_length
        two = (end - start + self.genome_length) % self.genome_length

        return one, two


# Data Holder to make deletion writing code more readable
class Read:
    def __init__(self):
        self.total_count = 0
        self.one_count = 0
        self.two_count = 0

        self.start = 0
        self.end = 0


class FileWriter:
    def __init__(self, file_path, header, write_directory, genome_length, reads_one, reads_two):
        """Writes the result deletions to one tab-delimited text file"""

        self.file_prefix = os.path.join(write_directory, os.path.basename(file_path))
        self.reads = self.collect_results(reads_one, reads_two)

        self.genome_length = genome_length
        self.header = header

        # Finding coverage to determine heteroplasmy level
        with pysam.AlignmentFile(file_path, 'rb') as samfile:
            self.coverage = [0 for i in range(genome_length+1)]

            for pileupcolumn in samfile.pileup(self.header, 0, self.genome_length):
                self.coverage[pileupcolumn.pos+1] = pileupcolumn.n

        self.write_results()

    @staticmethod
    def collect_results(reads_one, reads_two):
        """Combines strand one and two deletions into one file"""
        reads = dict()

        for read in reads_one:
            new_read = Read()

            new_read.total_count = read[0]
            new_read.one_count = read[0]
            new_read.start = read[1]
            new_read.end = read[2]

            reads[(read[1], read[2])] = new_read

        for read in reads_two:
            if (read[1], read[2]) not in reads:
                new_read = Read()

                new_read.total_count = int(read[0])
                new_read.two_count = int(read[0])
                new_read.start = read[1]
                new_read.end = read[2]

                reads[(read[1], read[2])] = new_read

            else:
                old_read = reads[(read[1], read[2])]

                old_read.two_count = int(read[0])
                old_read.total_count = old_read.one_count + old_read.two_count

        return reads

    def write_results(self):
        """Writes the reads to resulting text file"""

        def tab_line(*args):
            return '{}\n'.format('\t'.join(map(str, args)))

        with open('{}_results.txt'.format(self.file_prefix), 'w') as results:
            results.write(tab_line('Total Reads', 'S1 Reads', 'S2 Reads', 'Start', 'End',
                                   'Heteroplasmy Level'))

            # Extracts information from read data class
            fourth_length = self.genome_length // 4
            lower_cutoff = self.genome_length - self.genome_length/15
            upper_cutoff = self.genome_length/15

            for read in sorted(self.reads.values(), key=lambda r: r.total_count, reverse=True):
                start = int(read.start)
                end = int(read.end)
                if start < fourth_length  and end > fourth_length * 3:
                    read.start, read.end = read.end, read.start

                else:
                    read.start, read.end = min(start, end), max(start, end)
                
                start = int(read.start)
                end = int(read.end)

                if start > lower_cutoff and end < upper_cutoff:
                    continue

                if start < end:
                    avg_coverage = sum(self.coverage[start:end]) / (end-start+1)
                else:
                    length = self.genome_length - end + start + 1
                    avg_coverage = (sum(self.coverage[:start+1]) + sum(self.coverage[end:]))/length

                heteroplasmy = float(read.total_count) / (int(read.total_count) + avg_coverage)
                results.write(tab_line(read.total_count, read.one_count, read.two_count, read.start, read.end, heteroplasmy))


if __name__ == '__main__':
    from time import time
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('BAM', help='The BAM file that MitoMut will analyze', type=str)

    # Command line interface options w/ given defaults
    parser.add_argument('-c', '--header', help='Sets the mitochondrial DNA header', type=str, default='MT')
    parser.add_argument('-f', '--fasta', help='Sets the reference genome', type=str, default='mt.fasta')
    parser.add_argument('-q', '--qual', help='Specifies the required quality scores to be considered', type=int,
                        default=30)
    parser.add_argument('-s', '--support', help='Number of reads required to support a deletion', type=int, default=5)
    parser.add_argument('-e', '--extract', help='Extracts the mitochondrial section of the reference',
                        action='store_true')
    parser.add_argument('-si', '--size', help='Specifies the minimum number of missing bp to be considered a deletion',
                        type=int, default=5)
    parser.add_argument('-l', '--length', help='The assumed length of the mitochondrial genome', type=int,
                        default=16569)
    parser.add_argument('-d', '--directory', help='Sets the directory used to write result files', type=str, default='')

    args = parser.parse_args()

    if args.extract:
        new_fasta = '{}.mt.fasta'.format(args.fasta)
        subprocess.call('samtools faidx {} {} > {}'.format(args.fasta, args.header,
                                                           new_fasta), shell=True)
        args.fasta = new_fasta

    # Runs file through MitoMut process and times it
    print('BAM/SAM File   : {}\n'.format(args.BAM))

    print('Write Directory: {}'.format(args.directory))
    print('Reference File : {}'.format(args.fasta))
    print('Minimum Quality: {}'.format(args.qual))
    print('Minimum Reads  : {}'.format(args.support))
    print('Minimum Length : {}\n'.format(args.size))

    start = time()

    FileParser(args.BAM, args.header, args.fasta, args.directory, args.qual)
    reads = ReadParser(args.BAM, args.length, args.support, args.size, args.directory)
    FileWriter(args.BAM, args.header, args.directory, args.length, reads.reads_one, reads.reads_two)

    mins, secs = divmod(time() - start, 60)
    print('Elapsed Time   : {} mins, {} secs'.format(int(mins), int(secs)))
