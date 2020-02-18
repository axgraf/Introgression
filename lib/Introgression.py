# Created by Dr. Alexander Graf
# Mail: graf@genzentrum.lmu.de

import sys
import os
import time
import multiprocessing as mp
import queue
import pysam
from lib.VCF import SNP, VCF


class Introgression:

    def __init__(self, filtered_vcf_file_gzipped, window_size, step_size, binned_output_folder, threads):
        self.vcf_file = filtered_vcf_file_gzipped
        self.window_size = window_size  # window size to use
        self.step_size = step_size  # number of bases to jump for next calculation
        self.binned_output_folder = binned_output_folder  # output of window calculation; one file per chromosome
        self.threads = threads
        self.__setup()
        self.__tabix_vcf()

    def __setup(self):
        if not os.path.exists(self.binned_output_folder):
            os.mkdir(self.binned_output_folder)

    def __tabix_vcf(self):
        if os.path.basename(self.vcf_file).endswith('bgz'):
            if os.path.exists(self.vcf_file+".tbi"):
                vcf_stats = os.stat(self.vcf_file)
                tbi_stats = os.stat(self.vcf_file+".tbi")
                if tbi_stats.st_mtime < vcf_stats.st_mtime:
                    print("Index file is older than vcf file - recreate index")
                    pysam.tabix_index(self.vcf_file, preset="vcf", keep_original=True, force=True)
            else:
                pysam.tabix_index(self.vcf_file, preset="vcf", keep_original=True, force=True)

        else:
            sys.stderr.write("The File '{}' is not a bzip file".format(self.vcf_file))


    def sliding_window_generator(self, chrom_size):
        numOfChunks = ((chrom_size - self.window_size) / self.step_size) + 1
        for i in range(0, int(numOfChunks) * self.step_size, self.step_size):
            yield (i, i + self.window_size)

    def __create_empty_samples_genotype_dict(self, sample_names):
        alleles_window_sample = {}
        for sample in sample_names:
            alleles_window_sample.setdefault(sample, [0, 0])  # first item is snow sheep; second is reference
        return alleles_window_sample

    def determine_alleles(self, vcf_arr, sample_names):
        alleles_window_sample = self.__create_empty_samples_genotype_dict(sample_names)  # in case no variants were found
        for snp in vcf_arr:
            for sample_name, sample in snp.samples_dict.items():
                genotype = sample.GT.replace('/', '|').split('|')
                species_interest_alleles = 0
                reference_alleles = 0
                if genotype[0] == '1':
                    species_interest_alleles += 1
                else:
                    reference_alleles += 1
                if genotype[1] == '1':
                    species_interest_alleles += 1
                else:
                    reference_alleles += 1
                alleles_window_sample[sample_name][0] += species_interest_alleles
                alleles_window_sample[sample_name][1] += reference_alleles
        return alleles_window_sample

    def __write_window_to_file(self, alleles_window_sample_dict, window_writer, chrom, start_window, end_window, write_header):
        if write_header:
            window_writer.write("Chromosome\twindow_start\twindow_end\tgenotype\t")
            for idx, sample_name in enumerate(alleles_window_sample_dict.keys()):
                delim = "\t"
                if len(alleles_window_sample_dict) - 1 == idx:
                    delim = "\n"
                window_writer.write(sample_name + delim)
        is_write_row_info = True
        second_line = ""
        for idx, (sample_name, genotype_arr) in enumerate(alleles_window_sample_dict.items()):
            if is_write_row_info:
                window_writer.write("{}\t{}\t{}\t{}\t".format(chrom, start_window, end_window, "species_specific"))
                second_line += "{}\t{}\t{}\t{}\t".format(chrom, start_window, end_window, "reference")
                is_write_row_info = False
            delim = "\t"
            if len(alleles_window_sample_dict) - 1 == idx:
                delim = "\n"
            window_writer.write("{}{}".format(genotype_arr[0], delim))
            second_line += "{}{}".format(genotype_arr[1], delim)
        window_writer.write(second_line)

    def __process_chromosome(self, chromosome_queue, tabix_reader: pysam.TabixFile):
        vcf = VCF()
        samples = vcf.get_sample_names(self.vcf_file)
        while True:
            try:
                chromosome, size = chromosome_queue.get()
            except queue.Empty:
                time.sleep(0.1)
                continue
            if chromosome is None:
                break
            write_header = True

            window_writer = open(os.path.join(self.binned_output_folder, chromosome + "_window.csv"), 'w')
            chunks = self.sliding_window_generator(size)
            print("\nScreening: {}".format(chromosome))
            for start_pos, end_pos in chunks:
                records = tabix_reader.fetch(chromosome, start_pos, end_pos, multiple_iterators=True)
                vcf_arr = [SNP(line, samples) for line in list(records)]
                alleles_window_sample_dict = self.determine_alleles(vcf_arr, samples)
                self.__write_window_to_file(alleles_window_sample_dict, window_writer, chromosome, start_pos, end_pos, write_header)
                write_header = False
            window_writer.close()

    def calculate_introgression(self):
        vcf = VCF()
        tabix_reader = pysam.TabixFile(self.vcf_file, mode='r', encoding='utf-8')
        chromosome_queue = vcf.get_chromosome_queue(self.vcf_file)
        processes = []
        for i in range(0, self.threads):
            introgression_process = mp.Process(target=self.__process_chromosome, args=(chromosome_queue, tabix_reader,))
            introgression_process.daemon = True
            introgression_process.start()
            processes.append(introgression_process)

        for p in processes:
            p.join()

        tabix_reader.close()

