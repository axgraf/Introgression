# Created by Dr. Alexander Graf
# Mail: graf@genzentrum.lmu.de

import sys
import argparse
from lib.VCF import VCF
from lib.SNPFiltering import SNPFilter
from lib.Introgression import Introgression


def __check_sample_names_in_vcf(vcf_file, sample_names):
    extracted_names = VCF().get_sample_names(vcf_file)
    for sample in sample_names:
        if not sample in extracted_names:
            raise LookupError("Sample name '{}' not found in VCF file '{}'".format(sample, vcf_file))

if __name__ == '__main__':

    FILTER_CMD = "filter"
    WINDOW_CMD = "window"
    parser = argparse.ArgumentParser(
        description="Introgression analysis pipeline. First variants should be filtered to obtain species specific variants. " \
                    "Second the number of genotypes of species specific and reference alleles are calculated in sliding windows. "
        , epilog="author: Dr. Alexander Graf (graf@genzentrum.lmu.de)")
    subparsers = parser.add_subparsers(help="sub-command help", dest="command")

    filter_parser = subparsers.add_parser(FILTER_CMD,
                                          help="Filters species specific homozygous variants from a list of provided sample names. " \
                                               "An alternate variant must be homozygous for at least one of the species specific samples. " \
                                               "The reference samples must be homozygous for the reference alle in a given percentage of the reference samples [default: 90 %%]. \n" \
                                               "Further variants are filtered out if:\n " \
                                               "QD: < 2; FS > 60; SOR > 3; MQ < 40; MQRankSum < -12.5; ReadPosRankSum < -8.0.   \n"
                                               "INDELs are discarded!\n")
    filter_parser.add_argument('-v',
                               '--vcf_file',
                               metavar="FILE",
                               help="VCF.gz file with all variants [must be zipped]",
                               required=True)

    filter_parser.add_argument('-l',
                               '--sample_names',
                               metavar="[Sample1,Sample2,...]",
                               help="List of samples to filter on. Names should be provided as listed in the vcf file. "\
                                    "Names should be separated with comma.",
                               required=True)

    filter_parser.add_argument('-p', '--percentage_reference',
                               help="Percentage of reference samples that must be homozygous to the reference allele [default: 0.9]",
                               metavar="[0-1]",
                               type=float,
                               choices=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],
                               default=0.9,
                               required=False)

    filter_parser.add_argument('-n',
                               '--name_prefix',
                               metavar="PREFIX",
                               help="Prefix to prepend each output file",
                               default="species",
                               required=False)

    filter_parser.add_argument('-o',
                               '--output_folder',
                               metavar="FOLDER",
                               help="Output folder",
                               required=True)

    window_parser = subparsers.add_parser(WINDOW_CMD,
                                          help="Calculates the number of species specific alleles and reference alleles in defined window and step size for all samples. ")
    window_parser.add_argument('-v',
                               '--vcf_file',
                               metavar="FILE",
                               help="VCF.gz file with filtered variants from filter step",
                               required=True)

    window_parser.add_argument('-w', '--window_size',
                               metavar="INT",
                               help="Windows size to calculate alleles in bp [default: 70000]",
                               type=int,
                               default=70000,
                               required=False)

    window_parser.add_argument('-s', '--step_size',
                               metavar="INT",
                               help="Step size to move to next window in bp [default: 10000]",
                               type=int,
                               default=10000,
                               required=False)

    window_parser.add_argument('-o',
                               '--output_folder',
                               metavar="FOLDER",
                               help="Output folder to store results per chromosome",
                               required=True)

    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.command == FILTER_CMD:

        sample_names = args.sample_names.split(',')
        try:  # check if provided sample names match with those in VCF file
            __check_sample_names_in_vcf(args.vcf_file, sample_names)
        except LookupError as e:
            print(e)
            exit(1)

        snp_filter = SNPFilter(args.vcf_file, sample_names, args.name_prefix,  args.output_folder, args.percentage_reference)
        snp_filter.write_species_specific_vcf_gzip()

    elif args.command == WINDOW_CMD:
        introgression = Introgression(args.vcf_file, args.window_size, args.step_size, args.output_folder)
        introgression.calculate_introgression()
    else:
        parser.print_help(sys.stderr)
        sys.exit(1)

