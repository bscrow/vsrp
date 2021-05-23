"""
This is a bash pipeline that converts NetAffx annotation and sample files
(TSV format, from NA29 and/or NA32), into sorted, Gunzipped VCF files for
each chromosome - the required input file format for the Michigan
Imputation Server.

Steps:
------
1. Merge annotation files. If there is a conflict between the two files,
    NA32 version will be kept
2. Convert sample and annotation files from TSV format to PLINK formats
3. Perform quality control and convert to single chromosome VCF files using
    the pipeline recommended by Michigan Imputation Server.
4. Sort and gunzip VCF files using BCFTools

Requirements:
-------------
- Python 3.5 or newer
- plink and perl installed and executable added to path
- HRC-1000G-check-bim.pl (v4.3.0) & the tab delimited HRC reference,
    or wget and about 700Mb of disk space
- BCFTools installed and executable added to path
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys

from utils import merge_annotation, convert_to_plink, OUTDIR as TEMP_DIR

HRC_FILE = "https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip"
HRC_REF = "ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"


def parse_args():
    parser = argparse.ArgumentParser(
        description="This is a bash pipeline that converts NetAffx "
        "annotation and sample files (TSV format, from NA29 "
        "and/or NA32), into sorted, Gunzipped VCF files for"
        "each chromosome - the required input file format for "
        "the Michigan Imputation Server."
    )
    parser.add_argument(
        "dir32", help="Directory containing sample TSV files from NA32", type=str
    )
    parser.add_argument("anno32", help="TSV annotation file from NA32", type=str)
    parser.add_argument(
        "--dir29", help="Directory containing sample TSV files from NA29", type=str
    )
    parser.add_argument("--anno29", help="TSV annotation file for NA29", type=str)
    parser.add_argument(
        "--hrc",
        help="Directory containing HRC-1000G-check-bim-v4.3.0.zip & the "
        "tab delimited HRC reference. If not included, the required "
        "files will be downloaded into the current directory in hrc/.",
        type=str,
    )
    parser.add_argument("-o", help="Output directory of vcf.gz files", type=str)
    parser.add_argument("-r", help="Remove temporary files", action="store_true")
    args = parser.parse_args(args=None if sys.argv[1:] else ["--help"])

    # Checking HRC quality control files
    if args.hrc:
        if (
            not os.path.isdir(args.hrc)
            or "HRC-1000G-check-bim.pl" not in os.listdir(args.hrc)
            or "HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz" not in os.listdir(args.hrc)
        ):
            print("HRC files specified but not found. Exiting\n")
            return
    else:
        print("Downloading HRC files\n")
        if not os.path.exists("./hrc/"):
            os.mkdir("./hrc/")
        subprocess.run(["wget", HRC_FILE, "-P", "./hrc/"])
        subprocess.run(
            [
                "unzip",
                os.path.join("./hrc/", os.path.basename(HRC_FILE)),
                "-d",
                "./hrc/",
            ]
        )
        subprocess.run(["wget", HRC_REF, "-P", "./hrc/"])
        args.hrc = "./hrc/"

    # Step 1: Merge annotation files
    print("Merging annotation files\n")
    annotation = merge_annotation(args.anno32, args.anno29)

    # Step 2: Creating PLINK files
    print("Creating PLINK files\n")
    plink_file = convert_to_plink(args.dir32, args.dir29, annotation)

    # Step 3: Quality control via Michigan Imputation Server's pipeline
    print("Running Quality Control\n")
    subprocess.run(["plink", "--file", plink_file, "--make-bed", "--out", plink_file])
    subprocess.run(["plink", "--freq", "--bfile", plink_file, "--out", plink_file])
    subprocess.run(
        [
            "perl",
            os.path.join(args.hrc, "HRC-1000G-check-bim.pl"),
            "-b",
            plink_file + ".bim",
            "-f",
            plink_file + ".frq",
            "-r",
            os.path.join(args.hrc, "HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"),
            "-h",
        ]
    )
    subprocess.run(["sh", os.path.join(TEMP_DIR, "Run-plink.sh")])

    # Step 4: Sort and gunzip
    print("Preparing vcf.gz files\n")
    vcf_files = glob.glob(os.path.join(TEMP_DIR, "*.vcf"))
    if not os.path.exists(args.o):
        os.mkdir(args.o)
    for vcf_file in vcf_files:
        study_name, chrom = os.path.basename(vcf_file).split("-updated-")
        subprocess.run(
            [
                "bcftools",
                "sort",
                vcf_file,
                "-Oz",
                "-o",
                os.path.join(args.o, f"{study_name}_{chrom}.gz"),
            ]
        )

    # Remove intermediate files if desired
    if args.r:
        shutil.rmtree(TEMP_DIR, True)

    print(f"Done! Files are ready in {args.o}\n")


if __name__ == "__main__":
    parse_args()
