# Created by Dr. Alexander Graf
# Mail: graf@genzentrum.lmu.de

import gzip
import re
import multiprocessing as mp


class VCF:

    @staticmethod
    def get_sample_names(gvcf_file):
        with gzip.open(gvcf_file, 'rt') as vcf_reader:
            for line in vcf_reader:
                if line.startswith("#CHROM"):
                    tabs = line.rstrip().split('\t')
                    vcf_reader.close()
                    return tabs[9:]
            vcf_reader.close()

    @staticmethod
    def get_vcf_header(vcf_file):
        vcf_header = ""
        with gzip.open(vcf_file, 'rt') as vcf_reader:
            for line in vcf_reader:
                line = line.rstrip()
                if line.startswith('#'):
                    vcf_header += line + "\n"
                else:
                    break
            vcf_reader.close()
        return vcf_header

    @staticmethod
    def get_chromosome_sizes(vcf_file):
        chrom_size_dict = {}
        with gzip.open(vcf_file, 'rt') as vcf_reader:
            for line in vcf_reader:
                if line.startswith("#"):
                    if line.startswith("##contig"):
                        chrom_name = re.findall(r'ID=([A-Za-z0-9_\\.]+),', line)
                        if chrom_name:
                            chrom_name = chrom_name[0]
                        length = re.findall(r'length=([0-9]+)>', line)
                        if length:
                            length = length[0]
                        if chrom_name and length:
                            chrom_size_dict[chrom_name] = int(length)
                else:
                    break
            vcf_reader.close()
        return chrom_size_dict

    def get_chromosome_queue(self, vcf_file):
        chromosome_queue = mp.Queue()
        for chrom, size in self.get_chromosome_sizes(vcf_file).items():
            chromosome_queue.put((chrom, size))
        chromosome_queue.put((None, None))  # signal end of queue
        return chromosome_queue



class SNP:

    def __init__(self, vcf_line, samples):
        self.vcf_line = vcf_line
        self.samples_dict = self.__init_samples(samples)
        self.__tabs = vcf_line.split('\t')
        self.chrom = self.__tabs[0]
        self.position = int(self.__tabs[1])
        self.reference = self.__tabs[3]
        self.alternate = self.__tabs[4]
        self.filter = self.__tabs[6]
        self.__convert_info()
        self.__convert_sample(self.samples_dict)

    @staticmethod
    def __init_samples(samples):
        samples_dict = {}
        for s in samples:
            samples_dict[s] = {}
        return samples_dict

    def __convert_info(self):
        info = Info()
        info_line = self.__tabs[7]
        tabs = info_line.split(';')
        for entry in tabs:
            key_value = entry.split('=')
            if len(key_value) == 2:
                value = self.__parse_value(key_value[1])
                setattr(info, key_value[0], value)
        self.info = info

    def __convert_sample(self, samples_dict):
        format_line = self.__tabs[8]
        format_tabs = format_line.split(':')
        samples_values = self.__tabs[9:]  ## start of sample genotypes
        for idx, sample_name in enumerate(samples_dict.keys()):
            sample_tab = samples_values[idx].split(':')
            sample = Sample(sample_name)
            for idx_format, format in enumerate(format_tabs):
                if idx_format < len(sample_tab):
                    value = self.__parse_value(sample_tab[idx_format])
                    setattr(sample, format, value)
            samples_dict[sample_name] = sample

    def __parse_value(self, value):
        if value.isdigit():
            value = int(value)
        else:
            try:
                value = float(value)
            except ValueError:
                if ',' in value:
                    value = value.split(',')
                pass
        return value


class Info:

    def __init__(self):
        pass


class Sample:

    def __init__(self, name):
        self.name = name
        pass
