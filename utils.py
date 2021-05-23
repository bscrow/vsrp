"""
Python file containing util functions to convert input TSV files into a
PED and a MAP file for PLINK.
"""

import glob
import os

OUTDIR = "temp/"


def merge_annotation(na32, na29):
    """
    Merges annotation data of NA29 and NA32 based on Probe Set ID. If a
    Probe Set ID exists in both NA29 and NA32, the more up-to-date NA32
    annotation is kept.

    :param na32: NA32 annotation file
    :param na29: NA29 annotation file, or None if only NA32 annotation file
        is available
    :return: list of merged annotation data
    """
    if not (os.path.exists(na32) and os.path.isfile(na32)):
        print("NA32 Annotation file does not exist\n")
        return
    with open(na32, "r") as f:
        na32_anno = list(map(lambda x: x.strip().split("\t"), f.readlines()))
    if not na29:
        return na32_anno
    if not (os.path.exists(na29) and os.path.isfile(na29)):
        print("NA29 Annotation file does not exist\n")
        return
    with open(na29, "r") as f:
        na29_anno = list(map(lambda x: x.strip().split("\t"), f.readlines()))
    na29_anno.sort(key=lambda x: x[0])
    na32_anno.sort(key=lambda x: x[0])
    combined_anno = []
    i, j = 0, 0
    while i < len(na29_anno) and j < len(na32_anno):
        if na29_anno[i][0] < na32_anno[j][0]:
            combined_anno.append(na29_anno[i])
            i += 1
        elif na29_anno[i][0] == na32_anno[j][0]:
            # NA32 annotations override NA29 annotations
            combined_anno.append(na32_anno[j])
            i += 1
            j += 1
        else:
            combined_anno.append(na32_anno[j])
            j += 1
    combined_anno.extend(na29_anno[i:])
    combined_anno.extend(na32_anno[j:])
    return combined_anno


def update_genotype(ped_line, call_code, annotation):
    """
    Helper function that adds alleles for a SNP marker to a sample in
    the PED file, based on its call code.

    :param ped_line: list representing a line in the ped file
    :param call_code: string representing a call code for a marker
    :param annotation: list representing a sample from the annotation file
    """
    allele_cols = {"A": annotation[3], "B": annotation[4]}
    if call_code not in ("AA", "AB", "BA", "BB"):
        ped_line.append("0 0")
    else:
        ped_line.append(f"{allele_cols[call_code[0]]} {allele_cols[call_code[1]]}")


def convert_to_plink(na32_dir, na29_dir, anno, study_name="mystudy", outdir=OUTDIR):
    """
    Converts input tsv files into PED and MAP files for plink pipeline.
    These files will be named as <samples directory name>.ped/.map

    PED file: Each sample will be represented as a line in the file, with
        Family ID: study_name
        Individual ID: sample file name
        Paternal ID: 0 (missing)
        Maternal ID: 0 (missing)
        Sex: 0 (unknown)
        Phenotype: 0 (missing)
        Genotypes: Biallelic markers aligned to the order in the MAP file

    MAP file: Each SNP is represented by chromosome, rs ID and bp position.

    :param na32_dir: Directory of sample TSV files
    :param na29_dir: Directory of sample TSV files
    :param anno: file name of the annotation TSV file
    :param study_name: Name of study to name the output files. Default
        name is "mystudy"
    :param outdir: Name of directory to write output files to. Default
        name is "temp"

    :return: Name of PLINK files without extensions
    """
    samples = []
    if na29_dir:
        if not (os.path.exists(na29_dir) and os.path.isdir(na29_dir)):
            print("NA29 sample directory does not exist\n")
            return
        samples.extend(sorted(glob.glob(os.path.join(na29_dir, "*.tsv"))))
    if not (os.path.exists(na32_dir) and os.path.isdir(na32_dir)):
        print("NA32 sample directory does not exist\n")
        return
    samples.extend(sorted(glob.glob(os.path.join(na32_dir, "*.tsv"))))

    anno.sort(key=lambda x: x[0])  # sort by Probe Set ID
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    file_name = os.path.join(outdir, study_name)

    ped_content, map_content = [], []
    open(file_name + ".ped", "w").close()
    # For each sample, find alleles via probeset ID in annotation file and
    # write to ped file
    for sample_file in samples:
        ped_line = [study_name, os.path.split(sample_file)[1][:-4], "0", "0", "0", "0"]
        with open(sample_file, "r") as f:
            file = f.readlines()
        for i in range(len(file)):
            if not file[i].startswith("#"):
                break
        col_titles = file[i].strip().split("\t")
        probeset_id_col = col_titles.index("Probe Set ID")
        # chrom_col = col_titles.index("Chromosome")
        # pos_col = col_titles.index("Chromosomal Position")
        rsid_col = col_titles.index("dbSNP RS ID")
        call_code_col = col_titles.index("Call Codes")
        sample = list(map(lambda x: x.strip().split("\t"), file[6:]))
        sample.sort(key=lambda x: x[probeset_id_col])  # sort by Probe Set ID

        anno_pos, sample_pos = 0, 0
        while anno_pos < len(anno) and sample_pos < len(sample):
            if anno[anno_pos][0] == sample[sample_pos][probeset_id_col]:
                update_genotype(
                    ped_line, sample[sample_pos][call_code_col], anno[anno_pos]
                )
                # Append dbSNP RS ID to annotation list if it has not been added
                if not anno[anno_pos][-1].startswith("rs"):
                    anno[anno_pos].append(sample[sample_pos][rsid_col])
                anno_pos += 1
                sample_pos += 1
            elif anno[anno_pos][0] < sample[sample_pos][probeset_id_col]:
                anno_pos += 1
                ped_line.append("0 0")
            else:
                sample_pos += 1
        while anno_pos < len(anno):
            anno_pos += 1
            ped_line.append("0 0")

        with open(file_name + ".ped", "a") as f:
            for i in range(len(ped_line) - 1):
                f.write(ped_line[i] + "\t")
            f.write(ped_line[-1] + "\n")

    # compile lines to be written to the map file. Unused variants are
    # marked for reemoval in the map file by adding "-" to its chrom pos.
    for line in anno:
        if not line[-1].startswith("rs"):
            if not line[2].startswith("-"):
                line[2] = "-" + line[2]
        map_content.append("\t".join([line[1], line[-1], "0", line[2]]) + "\n")

    with open(file_name + ".map", "w") as f:
        f.writelines(map_content)

    return file_name


# def remove_unused_markers(mapfile):
#     """
#     Removes annotation markers that is not found in any of the samples.
#     Based on the implementation of convert_to_plink, if a marker is not
#     used at all, it will have Probe set ID in place of RS ID.
#
#     :param mapfile: path to the MAP file to filter.
#     """
#     with open(mapfile, "r") as f:
#         lines = list(map(lambda x: x.split("\t"), f.readlines()))
#     fo = open(mapfile, 'w')
#     for ln in lines:
#         if not ln[1].startswith("rs"):
#             if not ln[3].startswith("-"):
#                 ln[3] = "-" + ln[3]
#         fo.write("\t".join(ln))
#     fo.close()
#
#
# def map2bed(mapfile, bed):
#     with open(mapfile, "r") as f:
#         lines = f.readlines()
#     fo = open(bed, 'w')
#     for ln in lines:
#         chrom, rs, mdist, pos = ln.split("\t")
#         chrom = 'chr' + chrom
#         pos = int(pos)
#         fo.write(f"{chrom}\t{pos-1}\t{pos}\t{rs}\n")
#     fo.close()
#
#
# def bed2map(output, unlifted, mapfile, newmap):
#     with open(output, "r") as f:
#         olines = list(map(lambda x: x.strip().split("\t"), f.readlines()))
#     with open(unlifted, "r") as f:
#         ulines = list(map(lambda x: x.strip().split("\t"), f.readlines()))
#     mapping = dict()
#     for l in olines:
#         mapping[l[3]] = l[0][3:], l[2]
#     for l in ulines:
#         if l[0].startswith("#"):
#             continue
#         mapping[l[3]] = l[0][3:], "-"+l[2]
#     with open(mapfile, "r") as f:
#         mlines = f.readlines()
#     fo = open(newmap, 'w')
#     for ln in mlines:
#         chrom, rs, mdist, pos = ln.split("\t")
#         new_chrom, new_pos = mapping[rs]
#         fo.write(f"{new_chrom}\t{rs}\t{mdist}\t{new_pos}\n")
#     fo.close()
#
#
# def parse_args():
#     parser = argparse.ArgumentParser(
#         description="Converts sample and annotation files from TSV to "
#                     "plink compatible file formats (PED and MAP files), "
#                     "with the option to liftOver from hg18 to hg19"
#     )
#     parser.add_argument("dir", help="Directory containing sample TSV files", type=str)
#     parser.add_argument("anno", help="TSV annotation file", type=str)
#     parser.add_argument("dir29", help="Directory containing sample TSV files from NA29", type=str)
#     parser.add_argument("anno29", help="TSV annotation file for NA29", type=str)
#     parser.add_argument(
#         "--to-hg19",
#         help="Convert chromosome and position from hg18 to hg19",
#         action="store_true"
#     )
#     args = parser.parse_args(args=None if sys.argv[1:] else ["--help"])
#     if args.dir.endswith("\\") or args.dir.endswith("/"):
#         args.dir = args.dir[:-1]
#     if not os.path.exists(OUTDIR):
#         os.makedirs(OUTDIR)
#     print("Creating PLINK files...")
#     convert_to_plink(args.dir, args.anno, OUTDIR)
#
#     if args.to_hg19:
#         print("Performing liftOver from hg18 to hg19...")
#         map2bed(f"{OUTDIR}.map", f"{OUTDIR}.bed")
#         path = Path(__file__).parent.absolute()
#         subprocess.run([
#             f"{path}/liftOver/liftOver",
#             f"{OUTDIR}.bed",
#             "./liftOver/hg18ToHg19.over.chain.gz",
#             "output.bed",
#             "unlifted.bed"
#         ])
#         bed2map("output.bed", "unlifted.bed", f"{OUTDIR}.map", f"{OUTDIR}.lifted.map")
#         os.remove("output.bed")
#         os.remove("unlifted.bed")
#         os.remove(f"{OUTDIR}.bed")
#         os.remove(f"{OUTDIR}.map")
#         os.rename(f"{OUTDIR}.lifted.map", f"{OUTDIR}.map")
#
#     print("Removing unused markers...")
#     remove_unused_markers(f"{OUTDIR}.map")
#
#     print("PLINK files created")
#
#
# if __name__ == "__main__":
#     parse_args()
