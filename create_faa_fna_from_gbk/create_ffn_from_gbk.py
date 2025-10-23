from Bio import SeqIO # type: ignore
import argparse
import os
import sys
import itertools
from Bio.Seq import Seq # type: ignore
from pathlib import Path



def convert_gbk(genomes_folder, prefix, output_folder):
    output_subdir = os.path.join(output_folder, prefix)
    os.makedirs(output_subdir, exist_ok=True)

    root = Path(genomes_folder)
    gbk_paths = find_gbk_paths(root, prefix)
    if not gbk_paths:
        print(f"no .gbk files found that start with {prefix} and end with .gbk in {genomes_folder}")
        return
    # output .faa files to output subdir
    for gbk in gbk_paths: 
        print(f"\nProcessing: {gbk}")
        convert_to_ffn_output(gbk, output_subdir)
            
        print("\n" + "*-" * 40)


def find_gbk_paths(root, prefix):
    """
    First try flat search; if none, fall back to recursive.
    
    """
    flat = list(root.glob(f"{prefix}*.gbk"))  # no sorted()
    if flat:
        print(f"\n[info] Found {len(flat)} .gbk in {root} (flat search).")
        print("\n" + "*-" * 40)
        return flat

    rec = list(root.rglob(f"{prefix}*.gbk"))  # no sorted()
    if rec:
        print(f"\n[info] No flat matches; using recursive search. Found {len(rec)} under {root}.")
    else:
        print(f"\n[warn] No .gbk files starting with '{prefix}' found in {root} (flat or recursive).")
    print("\n" + "*-" * 40)
    return rec
                
                
def convert_to_ffn_output(gbk_pth, output_subdir):
    file = os.path.basename(gbk_pth) 
    ffn_file = os.path.join(output_subdir, file.replace(".gbk", ".fasta"))
    with open(gbk_pth, "r") as gbk_file: 
        for record in SeqIO.parse(gbk_file, "genbank"): 
            for feature in record.features: 
                if feature.type == "CDS": 
                    with open(ffn_file, "a") as ffn_out: 
                        locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                        translation = feature.qualifiers.get("translation",[""])[0]
                        protein_desc = feature.qualifiers.get("product", ["unknown"])[0]
                        if protein_desc in ["hypothetical protein", "unknown function"]: 
                            protein_desc = feature.qualifiers.get("interpro_product", ["hypothetical protein"])[0]
                        fasta_header = f"{locus_tag} {protein_desc} {record.id}"
                        cds_seq = feature.location.extract(record.seq)
                        protein_rec = SeqIO.SeqRecord(
                            Seq(str(cds_seq)), 
                            id=f"{locus_tag}",
                            description = fasta_header
                        )
                        SeqIO.write(protein_rec, ffn_out, "fasta")
                else: 
                    print(f"No translation for CDS {feature.get("locus_tag", "NULL")}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create .ffn files for all .gbk files found in passed directory.")

    parser.add_argument("--genomes_folder", required= True, help="Path to the input folder containing subdirectories with .gbk files.")
    parser.add_argument("--prefix", required= True, help="Prefix for the output file subdirectory")
    parser.add_argument("--output_folder", required= True, help="Path to the output folder where results will be saved.")
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    convert_gbk(args.genomes_folder, args.prefix, args.output_folder)