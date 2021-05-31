"""
Python file containing util functions to convert input TSV files into a
PED and a MAP file for PLINK.
"""

import glob
import os

OUTDIR = "temp/"


def merge_annotation(na32, na29, outfile):
    """
    Merges annotation data of NA29 and NA32 based on Probe Set ID. If a
    Probe Set ID exists in both NA29 and NA32, chromosome and position in
    NA32 annotation is kept. Alleles for NA29 and NA32 are both kept.

    Each row in the tab delimited output file is as follows:
    [Probe Set ID, rsID, Chromosome, Chromosomal Position,
    NA29 A Allele, NA29 B Allele, NA32 A Allele, NA32 B Allele]

    :param na32: NA32 annotation file
    :param na29: NA29 annotation file, or None if only NA32 annotation file
        is available
    :param outfile: Name of the output merged annotation file
    """
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "-": "-"}

    if not (os.path.exists(na32) and os.path.isfile(na32)):
        print("NA32 Annotation file does not exist\n")
        return
    if not (os.path.exists(na29) and os.path.isfile(na29)):
        print("NA29 Annotation file does not exist\n")
        return
    combined_anno = {}
    with open(na32, "r") as f:
        for line in f:
            if not line.startswith('"AX'):
                continue
            line = line.replace('"', "").split(",")
            if "---" in line[:3]:
                continue
            combined_anno[line[0]] = line[:4] + ["0", "0"] + line[8:10]
    with open(na29, "r") as f:
        for line in f:
            if not line.startswith('"AX'):
                continue
            line = line.replace('"', "").split(",")
            if "---" in line[:3]:
                continue
            if line[0] in combined_anno:
                # For now, remove probesets which are assigned different rsIDs
                # between NA29 and NA32
                if combined_anno[line[0]][1] != line[1]:
                    combined_anno.pop(line[0])
                    continue
                # If NA29 and NA32 have different alleles, check if this is
                # due to change in strand. If so, change NA29 strand
                if line[8:10] != combined_anno[line[0]][6:8]:
                    if combined_anno[line[0]][6:8] == [
                        complement[line[8]],
                        complement[line[9]],
                    ]:
                        line[8], line[9] = complement[line[8]], complement[line[9]]
                combined_anno[line[0]][4] = line[8]
                combined_anno[line[0]][5] = line[9]
            else:
                combined_anno[line[0]] = line[:4] + line[8:10] + ["0", "0"]
    with open(outfile, "w") as f:
        for entry in combined_anno:
            f.write("\t".join(combined_anno[entry]) + "\n")


# The default sorting algorithm is about 10X faster than this implementation
# def radix_sort(ls):
#     """
#     Helper function to sort a list based on the probeset ID at index 0,
#     which consist of a 8 digit number with a prefix "AX-".
#
#     :param ls: list to sort based on probeset ID at position pos
#     :return: sorted list
#     """
#     sorting_dict = {}
#     for i in range(10):
#         sorting_dict[str(i)] = []
#     for i in range(10, 2, -1):
#         for line in ls:
#             sorting_dict[line[0][i]].append(line)
#         ls = []
#         for j in range(10):
#             ls.extend(sorting_dict[str(j)])
#             sorting_dict[str(j)] = []
#     return ls


def update_genotype(ped_line, call_code, annotation, version):
    """
    Helper function that adds alleles for a SNP marker to a sample in
    the PED file, based on its call code.

    :param ped_line: list representing a line in the ped file
    :param call_code: string representing a call code for a marker
    :param annotation: list representing a sample from the annotation file
    :param version: string for sample annotation version (NA29/NA32)
    """
    allele_cols = {
        "NA29": {"A": annotation[4], "B": annotation[5]},
        "NA32": {"A": annotation[6], "B": annotation[7]},
    }
    if call_code not in ("AA", "AB", "BA", "BB"):
        ped_line.append("0 0")
    else:
        ped_line.append(
            f"{allele_cols[version][call_code[0]]} {allele_cols[version][call_code[1]]}"
        )


def convert_to_plink(
    samples_dir, group_file, anno_file, study_name="", outdir=OUTDIR
):
    """
    Converts input tsv files into PED and MAP files for plink pipeline.
    These files will be named as <samples directory name>.ped/.map

    PED file: Each sample will be represented as a line in the file, with
        Family ID: study_name, or sample file name if study_name is ""
        Individual ID: sample file name
        Paternal ID: 0 (missing)
        Maternal ID: 0 (missing)
        Sex: 0 (unknown)
        Phenotype: 1 (control) or 2 (case)
        Genotypes: Biallelic markers aligned to the order in the MAP file

    MAP file: Each SNP is represented by chromosome, rs ID and bp position.

    :param samples_dir: Directory of sample TSV files
    :param group_file: File specifying the phenotype group of each sample.
        If a sample file does not have an assigned group, the sample will
        be assigned 0 (missing)
    :param anno_file: file name of the annotation TSV file
    :param study_name: Name of study to name the output files. Default
        name is ""
    :param outdir: Name of directory to write output files to. Default
        name is "temp"

    :return: Name of PLINK files without extensions
    """
    if not (os.path.exists(samples_dir) and os.path.isdir(samples_dir)):
        print("Sample directory does not exist\n")
        return
    samples = glob.glob(os.path.join(samples_dir, "*.tsv"))
    if not (os.path.exists(group_file) and os.path.isfile(group_file)):
        print("Control sample directory does not exist\n")
        return
    groups = {}
    with open(group_file, "r") as f:
        for line in f:
            filename, group = line.strip().split("\t")
            groups[filename] = group
    if not (os.path.exists(anno_file) and os.path.isfile(anno_file)):
        print("Annotation file does not exist\n")
        return
    with open(anno_file, "r") as f:
        anno = list(map(lambda x: x.strip().split("\t"), f.readlines()))
    anno.sort(key=lambda x: x[0])

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    file_name = os.path.join(outdir, "raw")
    map_content = []
    open(file_name + ".ped", "w").close()
    # For each sample, find alleles via probeset ID in annotation file and
    # write to ped file
    for sample_file in samples:
        if os.path.basename(sample_file) in groups:
            group = groups[os.path.basename(sample_file)]
        else:
            group = "0"
        if not study_name:
            study_name = os.path.split(sample_file)[1][:-4]
        ped_line = [
            study_name,
            os.path.split(sample_file)[1][:-4],
            "0",
            "0",
            "0",
            group,
        ]
        version = "NA32"
        with open(sample_file, "r") as f:
            file = f.readlines()
        for i in range(len(file)):
            if file[i] == "#Axiom_GW_Hu_SNP.r2.na29.annot.db\n":
                version = "NA29"
            if not file[i].startswith("#"):
                break
        col_titles = file[i].strip().split("\t")
        probeset_id_col = col_titles.index("Probe Set ID")
        call_code_col = col_titles.index("Call Codes")
        sample = list(map(lambda x: x.strip().split("\t"), file[i + 1 :]))
        sample.sort(key=lambda x: x[probeset_id_col])  # sort by Probe Set ID

        anno_pos, sample_pos = 0, 0
        while anno_pos < len(anno) and sample_pos < len(sample):
            if anno[anno_pos][0] == sample[sample_pos][probeset_id_col]:
                update_genotype(
                    ped_line, sample[sample_pos][call_code_col], anno[anno_pos], version
                )
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

    # compile lines to be written to the map file.
    for line in anno:
        map_content.append("\t".join([line[2], line[1], "0", line[3]]) + "\n")

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
