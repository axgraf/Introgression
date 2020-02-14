# Created by Dr. Alexander Graf
# Mail: graf@genzentrum.lmu.de

import os
import gzip
from lib.VCF import SNP, VCF


class SNPFilter:

    def __init__(self, gvcf_file, snow_sheep_sample_names, prefix_output, output_folder, percentage_of_reference):
        self.gvcf_file = gvcf_file
        self.snow_sheep_sample_names = snow_sheep_sample_names
        self.prefix_output = prefix_output
        self.output_folder = output_folder
        self.percentage_of_reference = percentage_of_reference
        self.output_vcf = os.path.join(self.output_folder, self.prefix_output + "_species_specific.vcf.gz")
        self.log_file = os.path.join(self.output_folder, self.prefix_output + "_species_filter.log")
        self.__setup()

    def __setup(self):
        if not os.path.exists(self.output_folder):
            os.mkdir(self.output_folder)

    def write_snow_sheep_specific_vcf_gzip(self):
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
                    total_nb += 1
                    snp = SNP(line, samples)
                    if last_chr is None:
                        last_chr = snp.chrom
                    #  ##########  LOGGING ##########
                    if total_nb % 10000 == 0:
                        print("Chromosome:\t{}\tSpecific / Total:\t{}/{}".format(last_chr, specific_snv, total_nb), flush=True)
                    if last_chr != snp.chrom:
                        log_writer.write("Chromosome:\t{}\tSpecific / Total:\t{}/{}\n".format(last_chr, specific_snv, total_nb))
                        log_writer.flush()
                        print("\nChromosome:\t{}\tSpecific / Total:\t{}/{}".format(last_chr, specific_snv, total_nb), flush=True)
                        last_chr = snp.chrom
                        specific_snv = 0
                        total_nb = 0
                    #  ########## END LOGGING ##########
                    if self.__is_not_INDEL(snp):
                        if not self.__is_variant_filtered(snp):  # general variant quality filter
                            if self.__is_variant_specific(snp, self.snow_sheep_sample_names):  # only snow sheep specific ones are kept
                                output_writer.write((snp.vcf_line + "\n").encode('utf-8'))
                                specific_snv += 1
                                last_chr = snp.chrom
            log_writer.write("Chromosome:\t{}\tSpecific / Total:\t{}/{}\n".format(last_chr, specific_snv, total_nb))
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
                # print("Variant filtered: QD < 2: {}".format(snp.info.QD))
                return True
        if snp.info.FS > 60:  # strand bias test: alternative seen more or less on forward or reverse strand than reference allele
            # print("Variant filtered: FS > 60: {}".format(snp.info.FS))
            return True
        if snp.info.SOR > 3:  # strand odds ratio (strand bias test)
            # print("Variant filtered: SOR > 3: {}".format(snp.info.SOR))
            return True
        if hasattr(snp.info, 'MQ'):
            if snp.info.MQ < 40:  # root mean square mapping quality over all reads at site
                # print("Variant filtered: MQ < 40: {}".format(snp.info.MQ))
                return True
        if hasattr(snp.info, 'MQRankSum'):
            if snp.info.MQRankSum < -12.5:  # compare mapping qualities of reads supporting reference allele and alternate allele: >0: qualites of reads supporting the alternate allele are higher than those supporting reference allele. Values close to zero are best.
                # print("Variant filtered: MQRankSum < -12.5: {}".format(snp.info.MQRankSum))
                return True
        if hasattr(snp.info, 'ReadPosRankSum'):
            if snp.info.ReadPosRankSum < -8.0:  # compares whether the positions of the reference and alternate alleles are different within the reads. Seeing an allele only near the ends of reads is indicative of error, because that is where sequencers tend to make the most errors.
                # print("Variant filtered: ReadPosRankSum < -8.0: {}".format(snp.info.ReadPosRankSum))
                return True

    def __is_variant_specific(self, snp: SNP, snow_sheep_sample_names_arr):
        is_homozygote_species_interesst = False
        total_reference_samples = len(snp.samples_dict.keys()) - len(snow_sheep_sample_names_arr)
#        half_of_domestic_sheep_samples = total_domestic_sheep_samples / 2
        number_reference_samples = 0
        #        read_depths = []
        #        genotypes_domestic = []
        #        samples_names = []
        for sample_name, sample in snp.samples_dict.items():
            if sample_name in snow_sheep_sample_names_arr:
                if sample.GT == '1/1' or sample.GT == '1|1':
                    is_homozygote_species_interesst = True  # a single snow sheep must be homozygous for alternate allele
            else:
                #                if sample.GT == './.':
                #                    total_domestic_sheep_samples -= 1   # no genotype will not be counted
                if (sample.GT == '0/0' or sample.GT == '0|0'):  # and sample.DP >= 5:
                    #                    if hasattr(sample, 'GQ'):
                    #                        if sample.GQ is not ".":
                    #                            genotype_quality_confidence = 10 ** (-sample.GQ / 10)  # phred scaled
                    #                            if genotype_quality_confidence < 0.05:
                    number_reference_samples += 1
        #                genotypes_domestic.append(sample.GT)
        #                samples_names.append(sample_name)
        #                if hasattr(sample, 'DP'):
        #                    read_depths.append(sample.DP)
        #                else:
        #                    read_depths.append("NONE")
        if is_homozygote_species_interesst:  # and total_domestic_sheep_samples > half_of_domestic_sheep_samples:
            percentage_reference = number_reference_samples / total_reference_samples
            if percentage_reference > 0.9:  # at least 90% of the domestic sheeps must be homozygous for reference allele
                #                print("GOOD #######################")
                #                print("Genotypes:\t{}".format(' '.join(genotypes_domestic)))
                #                print("chr: {}\tpos: {}\tPercentage domestic:\t{}\ttotal found_domestic: {}".format(snp.chrom, snp.position, percentage_domestic_reference, total_domestic_sheep_samples))
                #                print("Read depth:\t{}".format(' '.join(map(str, read_depths))))
                #                print(" ")
                return True
        #           else:
        #               print("BAD #######################")
        #               print("\t\tGenotypes:\t{}".format(' '.join(genotypes_domestic)))
        #               print("\t\tchr: {}\tpos: {}\tPercentage domestic:\t{}\ttotal found_domestic: {}".format(snp.chrom, snp.position, percentage_domestic_reference, total_domestic_sheep_samples))
        #               print("\t\tRead depth:\t{}".format(' '.join(map(str, read_depths))))
        #               print(" ")
        return False
