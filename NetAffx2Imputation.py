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

    # Separate parser for merge
    subparsers = parser.add_subparsers(title="Sub-commands", dest="command")
    parser_merge = subparsers.add_parser("merge")
    parser_merge.add_argument("anno29", help="CSV annotation file from NA29", type=str)
    parser_merge.add_argument("anno32", help="CSV annotation file from NA32", type=str)
    parser_merge.add_argument(
        "-o",
        help="Name of the output merged annotation file. Default is 'merged_annotation.tsv'",
        type=str,
        default="merged_annotation.tsv",
    )
    parser_run = subparsers.add_parser("run")
    parser_run.add_argument("anno", help="TSV annotation file", type=str)
    parser_run.add_argument(
        "samples", help="Directory containing sample TSV files", type=str
    )
    parser_run.add_argument(
        "group", help="File specifying the phenotype group of each sample", type=str
    )
    parser_run.add_argument(
        "--hrc",
        help="Directory containing HRC-1000G-check-bim-v4.3.0.zip & the "
        "tab delimited HRC reference. If not included, the required "
        "files will be downloaded into the current directory in hrc/.",
        type=str,
    )
    parser_run.add_argument(
        "-t",
        help="Allele difference threshold value. Default is 0.2.",
        type=float,
        default=0.2,
    )
    parser_run.add_argument("-o", help="Output directory of vcf.gz files", type=str)
    parser_run.add_argument("-r", help="Remove temporary files", action="store_true")

    args = parser.parse_args(args=None if sys.argv[1:] else ["--help"])

    if args.command == "merge":
        merge_annotation(args.anno32, args.anno29, args.o)
    elif args.command == "run":
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
        # print("Merging annotation files\n")
        # annotation = merge_annotation(args.anno32, args.anno29)

        # Step 2: Creating PLINK files
        print("Creating PLINK files\n")
        plink_file = convert_to_plink(args.samples, args.group, args.anno)

        # Step 3: Quality control via Michigan Imputation Server's pipeline
        print("Running Quality Control\n")
        subprocess.run(
            [
                "plink",
                "--file",
                plink_file,
                "--make-bed",
                "--out",
                plink_file,
                "--allow-no-sex",
            ]
        )
        subprocess.run(
            [
                "plink",
                "--freq",
                "--bfile",
                plink_file,
                "--out",
                plink_file,
                "--allow-no-sex",
            ]
        )
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
                "-t",
                str(args.t),
            ]
        )
        # Allow phenotypes despite samples not having sex
        with open(os.path.join(TEMP_DIR, "Run-plink.sh"), "r") as f:
            cmds = f.readlines()
        for cmd in cmds:
            if cmd.startswith("plink"):
                cmd = cmd.strip() + " --allow-no-sex\n"
        with open(os.path.join(TEMP_DIR, "Run-plink.sh"), "w") as f:
            f.writelines(cmds)
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
