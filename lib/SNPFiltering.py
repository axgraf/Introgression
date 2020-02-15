# Created by Dr. Alexander Graf
# Mail: graf@genzentrum.lmu.de

import os
import sys
import gzip
from lib.VCF import SNP, VCF


class SNPFilter:

    def __init__(self, gvcf_file, samples_interest_names, prefix_output, output_folder, percentage_of_reference):
        self.gvcf_file = gvcf_file
        self.samples_interest_names = samples_interest_names
        self.prefix_output = prefix_output
        self.output_folder = output_folder
        self.percentage_of_reference = percentage_of_reference
        self.output_vcf = os.path.join(self.output_folder, self.prefix_output + "_specific_filtered.vcf.gz")
        self.log_file = os.path.join(self.output_folder, self.prefix_output + "_specific_filtered.log")
        self.__setup()

    def __setup(self):
        if not os.path.exists(self.output_folder):
            os.mkdir(self.output_folder)

    def __log(self, chromosome, specific_nb, total_nb, log_writer, log_to_file):
        if log_to_file: # assumes a switch in chromosomes -> reset counters
            log_writer.write("Chromosome:\t{}\t\tSpecific/Total:\t{}/{}\n".format(last_seen_chr, specific_nb, total_nb))
            log_writer.flush()
            sys.stdout.write("\rChromosome:\t{}\t\tSpecific/Total:\t{}/{}\n".format(chromosome, specific_nb, total_nb))
            return (0, 0)
        if total_nb % 10000 == 0:
            sys.stdout.write("\rChromosome:\t{}\t\tSpecific/Total:\t{}/{}".format(chromosome, specific_nb, total_nb))
            sys.stdout.flush()
        return (specific_nb, total_nb)

    def write_species_specific_vcf_gzip(self):
        vcf = VCF()
        samples = vcf.get_sample_names(self.gvcf_file)
        vcf_header = vcf.get_vcf_header(self.gvcf_file)
        total_nb = 0
        specific_snv = 0
        last_chr = None
        with gzip.open(self.gvcf_file, 'rt') as vcf_reader, gzip.open(self.output_vcf, 'wb') as output_writer, open(self.log_file, 'w') as log_writer:
            output_writer.write(vcf_header.encode('utf-8'))
            for line in vcf_reader:
                line = line.rstrip()
                if not line.startswith('#'):
                    snp = SNP(line, samples)
                    specific_snv, total_nb = self.__log(snp.chrom, specific_snv, total_nb, log_writer, False)
                    if last_chr != snp.chrom:
                        if last_chr != None:
                            specific_snv, total_nb = self.__log(snp.chrom, specific_snv, total_nb, log_writer, True)
                        last_chr = snp.chrom
                    if self.__is_not_INDEL(snp):
                        if not self.__is_variant_filtered(snp):  # general variant quality filter
                            if self.__is_variant_specific(snp, self.samples_interest_names):  # only species specific ones are kept
                                output_writer.write((snp.vcf_line + "\n").encode('utf-8'))
                                specific_snv += 1
                                last_chr = snp.chrom
                    total_nb += 1

            log_writer.write("Chromosome:\t{}\t\tSpecific/Total:\t{}/{}\n".format(last_chr, specific_snv, total_nb))
            log_writer.close()
            output_writer.close()
            vcf_reader.close()

    def __is_not_INDEL(self, snp: SNP):
        if len(snp.reference) == 1 and len(snp.alternate) == 1:
            return True
        else:
            return False

    @staticmethod
    def __is_variant_filtered(snp: SNP):  # return True if variant has bad qualities
        if hasattr(snp.info, 'QD'):
            if snp.info.QD < 2:  # quality by depth: variant confidence (from QUAL field) divided by the unfiltered depth of non-hom-ref samples. Intented to normalize the variant quality in order to avoid inflation caused when there is deep coverage.
                return True
        if snp.info.FS > 60:  # strand bias test: alternative seen more or less on forward or reverse strand than reference allele
            return True
        if snp.info.SOR > 3:  # strand odds ratio (strand bias test)
            return True
        if hasattr(snp.info, 'MQ'):
            if snp.info.MQ < 40:  # root mean square mapping quality over all reads at site
                return True
        if hasattr(snp.info, 'MQRankSum'):
            if snp.info.MQRankSum < -12.5:  # compare mapping qualities of reads supporting reference allele and alternate allele: >0: qualites of reads supporting the alternate allele are higher than those supporting reference allele. Values close to zero are best.
                return True
        if hasattr(snp.info, 'ReadPosRankSum'):
            if snp.info.ReadPosRankSum < -8.0:  # compares whether the positions of the reference and alternate alleles are different within the reads. Seeing an allele only near the ends of reads is indicative of error, because that is where sequencers tend to make the most errors.
                return True
        return False

    def __is_variant_specific(self, snp: SNP, names_of_species_of_interest_arr):
        is_homozygote_species_interest = False
        total_reference_samples = len(snp.samples_dict.keys()) - len(names_of_species_of_interest_arr)
        number_reference_samples = 0

        for sample_name, sample in snp.samples_dict.items():
            if sample_name in names_of_species_of_interest_arr:
                if sample.GT == '1/1' or sample.GT == '1|1':
                    is_homozygote_species_interest = True  # a single variant of the sample of interest must be homozygous for alternate allele
            elif sample.GT == '0/0' or sample.GT == '0|0':
                    number_reference_samples += 1
        if is_homozygote_species_interest:
            percentage_reference = number_reference_samples / total_reference_samples
            if percentage_reference > self.percentage_of_reference:  # default: at least 90% of the reference samples must be homozygous for reference allele
                return True
        return False
