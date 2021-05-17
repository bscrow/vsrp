"""
Command line program that takes in as argument a directory of tsv sample
files and an TSV annotation file and builds a PED and a MAP file for plink.
"""

import argparse
import glob
import os
from pathlib import Path
import subprocess
import sys


def update_genotype(ped_line, call_code, annotation):
    """
    Helper function that adds alleles for a SNP marker to a sample in
    the PED file, based on its call code.

    :param ped_line: list representing a line in the ped file
    :param call_code: string representing a call code for a marker
    :param annotation: list representing a sample from the annotation file
    """
    allele_cols = {
        "A": annotation[3],
        "B": annotation[4]
    }
    if call_code not in ("AA", "AB", "BA", "BB"):
        ped_line.append(f"0 0")
    else:
        ped_line.append(f"{allele_cols[call_code[0]]} {allele_cols[call_code[1]]}")


def convert_to_plink(samples_dir, anno, outdir):
    """
    Converts input tsv files into PED and MAP files for plink pipeline.
    These files will be named as <samples directory name>.ped/.map

    PED file: Each sample will be represented as a line in the file, with
        Family ID: samples directory name
        Individual ID: sample file name
        Paternal ID: 0 (missing)
        Maternal ID: 0 (missing)
        Phenotype: 0 (missing)
        Genotypes: Biallelic markers aligned to the order in the MAP file

    MAP file: Each SNP is represented by chromosome, rs ID and bp position.
    If rs ID is not available, Probe set ID is used instead.

    :param samples_dir: Directory of sample TSV files
    :param anno: file name of the annotation TSV file
    :param outdir: Name of directory to write output files in
    """
    if not (os.path.exists(samples_dir) and os.path.exists(anno)):
        raise Exception("Directory or Annotation files does not exist\n")
    ped_content, map_content = [], []
    samples = sorted(glob.glob(os.path.join(samples_dir, "*.tsv")))

    with open(anno, "r") as f:
        annotation = list(map(lambda x: x.strip().split("\t"), f.readlines()))
    annotation.sort(key=lambda x: (x[1], x[2]))  # sort by chromosome followed by position
    fam_id = os.path.basename(samples_dir)
    open(outdir + fam_id + ".ped", "w").close()
    for sample_file in samples:
        ped_line = [fam_id, os.path.split(sample_file)[1][:-4], "0", "0", "0", "0"]

        with open(sample_file, "r") as f:
            file = f.readlines()
        sample_cols = file[5].strip().split("\t")
        chrom_col = sample_cols.index("Chromosome")
        pos_col = sample_cols.index("Chromosomal Position")
        rsid_col = sample_cols.index("dbSNP RS ID")
        call_code_col = sample_cols.index("Call Codes")
        sample = list(map(lambda x: x.strip().split("\t"), file[6:]))
        sample.sort(key=lambda x: (x[chrom_col], x[pos_col]))  # sort by chromosome followed by position

        anno_pos, sample_pos = 0, 0
        while anno_pos < len(annotation) and sample_pos < len(sample):
            if (annotation[anno_pos][1], annotation[anno_pos][2]) == (sample[sample_pos][chrom_col], sample[sample_pos][pos_col]):
                update_genotype(ped_line, sample[sample_pos][call_code_col], annotation[anno_pos])
                # Replace probe set id with dbSNP RS ID in annotation list
                annotation[anno_pos][0] = sample[sample_pos][rsid_col]
                anno_pos += 1
                sample_pos += 1
            elif (annotation[anno_pos][1], annotation[anno_pos][2]) < (sample[sample_pos][chrom_col], sample[sample_pos][pos_col]):
                anno_pos += 1
                ped_line.append("0 0")
            else:
                sample_pos += 1

        with open(outdir + fam_id + ".ped", "a") as f:
            for i in range(len(ped_line)-1):
                f.write(ped_line[i] + "\t")
            f.write(ped_line[-1] + "\n")

    for line in annotation:
        map_content.append("\t".join([line[1], line[0], "0", line[2]]) + "\n")

    with open(outdir + fam_id + ".map", "w") as f:
        f.writelines(map_content)


def map2bed(mapfile, bed):
    with open(mapfile, "r") as f:
        lines = f.readlines()
    fo = open(bed, 'w')
    for ln in lines:
        chrom, rs, mdist, pos = ln.split("\t")
        chrom = 'chr' + chrom
        pos = int(pos)
        fo.write(f"{chrom}\t{pos-1}\t{pos}\t{rs}\n")
    fo.close()


def bed2map(output, unlifted, mapfile, newmap):
    with open(output, "r") as f:
        olines = list(map(lambda x: x.strip().split("\t"), f.readlines()))
    with open(unlifted, "r") as f:
        ulines = list(map(lambda x: x.strip().split("\t"), f.readlines()))
    mapping = dict()
    for l in olines:
        mapping[l[3]] = l[0][3:], l[2]
    for l in ulines:
        if l[0].startswith("#"):
            continue
        mapping[l[3]] = l[0][3:], "-"+l[2]
    with open(mapfile, "r") as f:
        mlines = f.readlines()
    fo = open(newmap, 'w')
    for ln in mlines:
        chrom, rs, mdist, pos = ln.split("\t")
        new_chrom, new_pos = mapping[rs]
        fo.write(f"{new_chrom}\t{rs}\t{mdist}\t{new_pos}\n")
    fo.close()


def remove_unused_markers(mapfile):
    """
    Removes annotation markers that is not found in any of the samples.
    Based on the implementation of convert_to_plink, if a marker is not
    used at all, it will have Probe set ID in place of RD ID.

    :param mapfile: path to the MAP file to filter.
    """
    with open(mapfile, "r") as f:
        lines = list(map(lambda x: x.split("\t"), f.readlines()))
    fo = open(mapfile, 'w')
    for ln in lines:
        if not ln[1].startswith("rs"):
            if not ln[3].startswith("-"):
                ln[3] = "-" + ln[3]
        fo.write("\t".join(ln))
    fo.close()


# def get_chrom_list(mapfile):
#     with open(mapfile, "r") as f:
#         lines = list(set(map(lambda x: x.split("\t")[0] + "\n", f.readlines())))
#     with open(mapfile[:-4]+"_chromosomes.txt", "w") as f:
#         f.writelines(lines)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Converts sample and annotation files from TSV to "
                    "plink compatible file formats (PED and MAP files), "
                    "with the option to liftOver from hg18 to hg19"
    )
    parser.add_argument("dir", help="Directory containing sample TSV files", type=str)
    parser.add_argument("anno", help="TSV annotation file", type=str)
    parser.add_argument(
        "--to-hg19",
        help="Convert chromosome and position from hg18 to hg19",
        action="store_true"
    )
    args = parser.parse_args(args=None if sys.argv[1:] else ["--help"])
    if args.dir.endswith("\\") or args.dir.endswith("/"):
        args.dir = args.dir[:-1]
    fam_id = os.path.basename(args.dir)
    outdir = fam_id + "plink/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    print("Creating PLINK files...")
    convert_to_plink(args.dir, args.anno, outdir)

    if args.to_hg19:
        print("Performing liftOver from hg18 to hg19...")
        map2bed(f"{outdir + fam_id}.map", f"{outdir + fam_id}.bed")
        path = Path(__file__).parent.absolute()
        subprocess.run([
            f"{path}/liftOver/liftOver",
            f"{outdir + fam_id}.bed",
            "./liftOver/hg18ToHg19.over.chain.gz",
            "output.bed",
            "unlifted.bed"
        ])
        bed2map("output.bed", "unlifted.bed", f"{outdir + fam_id}.map", f"{outdir + fam_id}.lifted.map")
        os.remove("output.bed")
        os.remove("unlifted.bed")
        os.remove(f"{outdir + fam_id}.bed")
        os.remove(f"{outdir + fam_id}.map")
        os.rename(f"{outdir + fam_id}.lifted.map", f"{outdir + fam_id}.map")

    print("Removing unused markers...")
    remove_unused_markers(f"{outdir + fam_id}.map")
    # get_chrom_list(f"{outdir + fam_id}.map")

    print("Done!")


if __name__ == "__main__":
    parse_args()
