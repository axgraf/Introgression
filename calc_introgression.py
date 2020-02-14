# Created by Dr. Alexander Graf
# Mail: graf@genzentrum.lmu.de

import argparse
import gzip
from lib.SNPFiltering import SNPFilter
from lib.Introgression import Introgression

FILTER_CMD = "filter"
WINDOW_CMD = "window"

input_genotpyed_gvcf = '/data/local/medugorac/alex/all_sheeps_combined.genotyped.g.vcf.gz'
snow_sheep_samples = ['Ov5904', 'Ov6748', '100610', '100611', '100612', '100613']
output_tmp_vcf_folder = '/data/local/medugorac/alex/vcf_snp_filter_per_chromosome'
output_vcf_snow_sheep_specific = '/data/local/medugorac/alex/all_sheeps_combined.genotyped.snow_sheep_specific.vcf.gz'
#output_vcf_snow_sheep_specific = '/data/local/medugorac/alex/test.snp.gz'
output_log_file = '/data/local/medugorac/alex/introgression_snp_filter.log'

output_window_file = '/data/local/medugorac/alex/window_genotypes'

#snp_filter = SNPFilter(input_genotpyed_gvcf, snow_sheep_samples, output_tmp_vcf_folder, output_vcf_snow_sheep_specific, output_log_file)
#snp_filter.write_snow_sheep_specific_vcf_gzip()


snp_filter = SNPFilter(input_genotpyed_gvcf, snow_sheep_samples, output_tmp_vcf_folder, output_vcf_snow_sheep_specific, output_log_file)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Introgression analysis pipeline"
        , epilog="author: Dr. Alexander Graf (graf@genzentrum.lmu.de)")
    subparsers = parser.add_subparsers(help="sub-command help", dest="command")

    ################### RUN #########################
    filter_parser = subparsers.add_parser(FILTER_CMD, help="Filters species specific homozygous variants from a list of provided samples names. A variant must be homozygous for at least one of the "
                    "species specific samples. The reference samples must be homozygous for the reference alles with at least n % [default: 90%].\n")

    filter_parser.add_argument('-v',
                               '--vcf_file',
                               metavar="FILE",
                               help="VCF.gz file with all variants [must be zipped]",
                               required=True)
    filter_parser.add_argument('-l',
                               '--species_names',
                               metavar="FILE",
                               help="List of species to filter on. Names should be provided as listed in the vcf file. "
                                    "Names should be separated with comma. e.g. Sample1,Sample2,...",
                               required=True)
    filter_parser.add_argument('-p', '--percentage_reference', metavar="[0-1]",
                               help="Windows size to calculate alleles in bp [default: 0.9]",
                               type=float,
                               choices=range(0, 1),
                               default=0.9,
                               required=False)
    filter_parser.add_argument('-o',
                               '--output_folder',
                               metavar="FOLDER",
                               help="Output folder",
                               required=True)

    window_parser = subparsers.add_parser(WINDOW_CMD,
                                           help="Filter: Filters species specific homozygous variants from a list of provided samples names. A variant must be homozygous for at least one of the "
                                                "species specific samples. The reference samples must be homozygous for the reference alles with at least n % [default: 90%].\n")
    window_parser.add_argument('-v',
                               '--vcf_file',
                               metavar="FILE",
                               help="VCF.gz file with filtered variants from filter step",
                               required=True)
    window_parser.add_argument('-l',
                               '--species_names',
                               metavar="FILE",
                               help="List of species to filter on. Names should be provided as listed in the vcf file. "
                                    "Names should be separated with comma. e.g. Sample1,Sample2,...",
                               required=True)
    window_parser.add_argument('-w', '--window_size', metavar="INT",
                               help="Windows size to calculate alleles in bp [default: 70000]",
                               type=int,
                               default=70000,
                               required=False)
    window_parser.add_argument('-s', '--step_size', metavar="INT",
                               help="Step size to move to next window in bp [default: 10000]",
                               type=int,
                               default=10000,
                               required=False)
    window_parser.add_argument('-o',
                               '--output_folder',
                               metavar="FOLDER",
                               help="Output folder",
                               required=True)

    args = parser.parse_args()

    if args.command == FILTER_CMD:
        snp_filter = SNPFilter(input_genotpyed_gvcf, snow_sheep_samples, output_tmp_vcf_folder, output_vcf_snow_sheep_specific, output_log_file)
        snp_filter.write_snow_sheep_specific_vcf_gzip()
    elif args.command == WINDOW_CMD:
        #introgression = Introgression(output_vcf_snow_sheep_specific, snow_sheep_samples, 70000, 10000, output_window_file)
        introgression = Introgression(args.vcf_file, args.species_names, args.window_size, args.step_size, args.output_folder)
        introgression.calculate_introgression()

