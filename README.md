# vsrp

(Work in Progress)

The proposed pipeline builds on top of the data processing pipeline as 
described on [Michigan Imputation Server Docs](https://imputationserver.readthedocs.io/en/latest/prepare-your-data/)

Sample and annotation TSV files are converted into PLINK file formats, 
lifting over to hg19 if necessary, before passing into the aforementioned 
quality control pipeline. Subsequently, the PLINK files are split based on 
chromosome and converted to sorted vcf.gz file via plink and bcftools. 


Usage:

```
tsv2plink.py [-h] [--to-hg19] dir anno

Converts sample and annotation files from tsv to plink compatible file formats
(ped and map files), with the option to liftOver from hg18 to hg19

positional arguments:
  dir         Directory containing sample tsv files
  anno        Tsv annotation file

optional arguments:
  -h, --help  show this help message and exit
  --to-hg19   Convert chromosome and position from hg18 to hg19
```