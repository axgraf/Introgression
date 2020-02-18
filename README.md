## Install

The tool is written in Python 3. It requires following packages:

* `pysam `


## Running the Introgression pipeline

### Obtaining VCF file

Please follow the 'best practices workflow' for GATK using the Germline short variant discovery  

### Filter

 Filters species specific homozygous variants from a list of provided sample names. 
 An alternate variant must be homozygous for at least one of the species specific 
 samples. The reference samples must be homozygous for the reference alle in a 
 given percentage of the reference samples [default: 90 %]. 
 
 Further variants are filtered out if one of the following criteria apply: 
 
 ```
 QD: < 2 
 FS > 60; 
 SOR > 3; 
 MQ < 40; 
 MQRankSum < -12.5; 
 ReadPosRankSum < -8.0
 ``` 
 
 INDELs are discarded!



### Window

Calculates the number of species specific alleles 
and reference alleles in defined window and step size for all samples.
