# Created by Dr. Alexander Graf
# Mail: graf@genzentrum.lmu.de

import sys
import os
import gzip
import numpy as np
from itertools import islice
from lib.VCF import SNP, VCF


class Introgression:

    def __init__(self, filtered_vcf_file_gzipped, snow_sheep_sample_names, window_size, step_size, binned_output_folder):
        self.vcf_file = filtered_vcf_file_gzipped
        self.snow_sheep_sample_names = snow_sheep_sample_names
        self.window_size = window_size  # window size to use
        self.step_size = step_size  # number of bases to jump for next calculation
        self.binned_output_folder = binned_output_folder  # output of window calculation; one file per chromosome
        self.__setup()

    def __setup(self):
        if not os.path.exists(self.binned_output_folder):
            os.mkdir(self.binned_output_folder)

    def sliding_window_generator(self, chrom_size):
        numOfChunks = ((chrom_size - self.window_size) / self.step_size) + 1
        for i in range(0, int(numOfChunks) * self.step_size, self.step_size):
            yield (i, i + self.window_size)

#    def __get_vcfs_in_range(self, vcf_arr, min_pos, max_pos, offset):
#        vcf_in_range_arr = []
#        last_idx = 0
#        for idx, vcf in enumerate(islice(vcf_arr, offset, None)):
#            if min_pos <= vcf.position <= max_pos:
#                vcf_in_range_arr.append(vcf)
#                last_idx = idx
#            if vcf.position > max_pos:  # vcf file is sorted by position, thus values exceeding max_pos can be discarded
#                break
#        return last_idx + offset, vcf_in_range_arr

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
                snow_sheep_alleles = 0
                reference_alleles = 0
                if genotype[0] == '1':
                    snow_sheep_alleles += 1
                else:
                    reference_alleles += 1
                if genotype[1] == '1':
                    snow_sheep_alleles += 1
                else:
                    reference_alleles += 1
                alleles_window_sample[sample_name][0] += snow_sheep_alleles
                alleles_window_sample[sample_name][1] += reference_alleles
        return alleles_window_sample


    def __write_window_to_file(self, alleles_window_sample_dict, window_writer, chrom, start_window, end_window, write_header):
        if write_header:
            window_writer.write("Chromosome\twindow_start\twindow_end\tgenotype\t")
            for idx, sample_name in enumerate(alleles_window_sample_dict.keys()):
                delim = "\t"
                if len(alleles_window_sample_dict) - 1 == idx:
                    delim = "\n"
#                window_writer.write(sample_name + "_snow_sheep\t")
                window_writer.write(sample_name + delim)
        is_write_row_info = True
        second_line = ""
        for idx, (sample_name, genotype_arr) in enumerate(alleles_window_sample_dict.items()):
            if is_write_row_info:
                window_writer.write("{}\t{}\t{}\t{}\t".format(chrom, start_window, end_window, "snow_sheep"))
                second_line += "{}\t{}\t{}\t{}\t".format(chrom, start_window, end_window, "reference")
                is_write_row_info = False
            delim = "\t"
            if len(alleles_window_sample_dict) - 1 == idx:
                delim = "\n"
            window_writer.write("{}{}".format(genotype_arr[0], delim))
            window_writer.flush()
            second_line += "{}{}".format(genotype_arr[1], delim)
        window_writer.write(second_line)
        window_writer.flush()


    def calculate_introgression(self):
        vcf = VCF()
        seek_position = 0
        chrom_sizes_dict = vcf.get_chromosome_sizes(self.vcf_file)

#        vcf_chrom_dict = vcf.read_vcf_file_per_chrom(self.vcf_file)
        last_idx = 0
        samples = vcf.get_sample_names(self.vcf_file)
        write_header = True
        with gzip.open(self.vcf_file, 'rt') as vcf_reader:
            try:
                for chromosome, size in chrom_sizes_dict.items():
                    window_writer = open(os.path.join(self.binned_output_folder, chromosome+"_window.csv"), 'w')
                    chunks = self.sliding_window_generator(size)
                    number_of_snps_per_window = []
                    vcf_arr = []
                    print("\nScreening: {}".format(chromosome))
                    for start_pos, end_pos in chunks:
                        vcf_arr, seek_position = vcf.read_lazy_chunks_per_chromsome(vcf_reader, samples, chromosome, start_pos, end_pos, seek_position, vcf_arr)
                       # last_idx, vcf_in_range = self.__get_vcfs_in_range(vcf_arr, start_pos, end_pos, last_idx)
                        sys.stdout.write("\rWindow:{}:{}-{}\tNb.SNPs: {}".format(chromosome, start_pos, end_pos,len(vcf_arr)))
                        sys.stdout.flush()
                        alleles_window_sample_dict = self.determine_alleles(vcf_arr, samples)
                        self.__write_window_to_file(alleles_window_sample_dict, window_writer, chromosome, start_pos, end_pos, write_header)
                        write_header = False
                        number_of_snps_per_window.append(len(vcf_arr))
                        #print("First/Last : {}:{}\t{}:{}".format(vcf_arr[0].chrom, vcf_arr[0].position, vcf_arr[-1].chrom, vcf_arr[-1].position))
                        if not seek_position:  # end of file reached
                            window_writer.close()
                            raise EOFError('End of file reached')  # a bit clumsy but outer loop needs to be broken
                    number_of_snps_per_window = np.asarray(number_of_snps_per_window)
                    np.savetxt("/data/local/medugorac/alex/nb_snps_{}.csv".format(chromosome), number_of_snps_per_window, delimiter="\t")
                    window_writer.close()
            except EOFError as e:
                print("\nDone")
        vcf_reader.close()


  #  def __calculate_window(self, idx, ):

