# Created by Dr. Alexander Graf
# Mail: graf@genzentrum.lmu.de

import sys
import os
import gzip
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
            window_writer.flush()
            second_line += "{}{}".format(genotype_arr[1], delim)
        window_writer.write(second_line)
        window_writer.flush()


    def calculate_introgression(self):
        vcf = VCF()
        seek_position = 0
        chrom_sizes_dict = vcf.get_chromosome_sizes(self.vcf_file)
        last_idx = 0
        samples = vcf.get_sample_names(self.vcf_file)
        write_header = True
        with gzip.open(self.vcf_file, 'rt') as vcf_reader:
            try:
                for chromosome, size in chrom_sizes_dict.items():
                    window_writer = open(os.path.join(self.binned_output_folder, chromosome+"_window.csv"), 'w')
                    chunks = self.sliding_window_generator(size)
                    vcf_arr = []
                    print("\nScreening: {}".format(chromosome))
                    for start_pos, end_pos in chunks:
                        vcf_arr, seek_position = vcf.read_lazy_chunks_per_chromsome(vcf_reader, samples, chromosome, start_pos, end_pos, seek_position, vcf_arr)
                       # last_idx, vcf_in_range = self.__get_vcfs_in_range(vcf_arr, start_pos, end_pos, last_idx)
                        sys.stdout.write("\r\tWindow:\t{}:{}-{}\tNb.SNPs: {}".format(chromosome, start_pos, end_pos,len(vcf_arr)))
                        sys.stdout.flush()
                        alleles_window_sample_dict = self.determine_alleles(vcf_arr, samples)
                        self.__write_window_to_file(alleles_window_sample_dict, window_writer, chromosome, start_pos, end_pos, write_header)
                        write_header = False
                        if not seek_position:  # end of file reached
                            window_writer.close()
                            raise EOFError('End of file reached')  # a bit clumsy but outer loop needs to be broken
                    window_writer.close()
            except EOFError as e:
                print("\n")
        vcf_reader.close()

