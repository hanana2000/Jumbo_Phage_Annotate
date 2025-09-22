import subprocess
import os
import argparse 
import sys
from Bio import SeqIO # type: ignore
from pathlib import Path


def retrieve_from_gbk(input_tsv, gbk1, gbk2, output_folder): 
    os.makedirs(output_folder, exist_ok=True) 
    output_tsv = os.path.join(output_folder, "SUMMARY_tsv_gbk_retreive_results.txt")
    gbk1_file_name = os.path.basename(gbk1).replace(".gbk", "").replace("_colour", "").replace("_NEW", "")
    gbk2_file_name = os.path.basename(gbk2).replace(".gbk", "").replace("_colour", "").replace("_NEW", "")

    with open(input_tsv, "r") as tsv_file, open(output_tsv, "w") as out_file:
        # create a feature library for gbk1 and gbk2
        feature_lib1 = create_feature_library(gbk1)
        feature_lib2 = create_feature_library(gbk2)
        # iterate through the tsv file hits
        for index, line in enumerate(tsv_file):
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            query_coords = (fields[6], fields[7])
            subject_coords = (fields[8], fields[9])
            hit_len, pident, evalue = int(fields[3]), float(fields[2]), fields[10]
            # Retrieve sequences from gbk files
            seq1 = retrieve_sequence_from_gbk(feature_lib1, query_coords)
            seq2 = retrieve_sequence_from_gbk(feature_lib2, subject_coords)
            if seq1 and seq2:
                out_file.write(f"\n******************************\nHit {index + 1}:\n\thit length: {hit_len}\n\tpercent identity: {pident}\n\tevalue: {evalue}\n******************************\n\n")
                out_file.write(f"for {gbk1_file_name}: \n\tcoords: {query_coords}")
                out_file.write(f"{seq1}\n")
                out_file.write(f"for {gbk2_file_name}: \n\tcoords: {subject_coords}")
                out_file.write(f"{seq2}\n")
            else: 
                out_file.write(f"\nCoordinates {query_coords} or {subject_coords} not found in respective gbk file coords.\n")
            out_file.write("-" * 50 + "\n")


def create_feature_library(gbk_file):
    feature_lib = {}
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features: 
            if feature.type == "CDS":
                ID = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                aa_length = len(feature.qualifiers.get("translation", [""])[0])
                product = feature.qualifiers.get("product", ["unknown"])[0]
                if "interpro_product" in feature.qualifiers: 
                    product += " | " + feature.qualifiers["interpro_product"][0]
                start = int(feature.location.start) + 1 # convert to 1-based
                end = int(feature.location.end)
                complement = False 
                if feature.location.strand == -1: # handle reverse strand
                    complement = True 
                feature_lib[product, (start,end)] = (start, end, complement, ID, aa_length)
    return feature_lib


def retrieve_sequence_from_gbk(gbk_lib, coords):
    out_message = ""
    start_save = []
    end_save = []
    hit_loc_start = []
    hit_loc_end = []
    for feature, (start, end, complement, ID, aa_length) in gbk_lib.items(): 
        if int(coords[0]) >= start and int(coords[0]) <= end:
            hit_loc_start.append(ID)
            start_save.append(feature[0])
            start_save.append((start, end))
            start_save.append("complement" if complement else "forward")
            start_save.append(ID)
            start_save.append(aa_length)
        #  out_message += f"\nhit started in\n\t {feature[0]} \n\t{start}-{end} \n\t{'complement' if complement else 'forward'} \n\t{ID}"
        if int(coords[1]) >= start and int(coords[1]) <= end:
            hit_loc_end.append(ID)
            end_save.append(feature[0])
            end_save.append((start, end))
            end_save.append("complement" if complement else "forward")
            end_save.append(ID)
            end_save.append(aa_length)
        #  out_message += f"\nhit ended in\n\t {feature[0]} \n\t{start}-{end} \n\t{'complement' if complement else 'forward'} \n\t{ID}\n"
    # right now this ONLY DISPLAYS THE FIRST HIT START AND END, even if multiple
    # can expand on later if needed
    if hit_loc_start == hit_loc_end and hit_loc_start != []:
        out_message += f"\nhit started and ended in \n\t{start_save[0]} \n\tlength: {start_save[4]} \n\t{start_save[1][0]}-{start_save[1][1]} \n\t{start_save[2]} \n\t{start_save[3]}\n"
    else: 
        out_message += f"\nhit started in \n\t{start_save[0]} \n\tlength: {start_save[4]} \n\t{start_save[1][0]}-{start_save[1][1]} \n\t{start_save[2]} \n\t{start_save[3]}"
        out_message += f"\nhit ended in \n\t{end_save[0]} \n\tlength: {end_save[4]} \n\t{end_save[1][0]}-{end_save[1][1]} \n\t{end_save[2]} \n\t{end_save[3]}\n"  
    return out_message if out_message else None



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="retreive all the hits from a tsv (tBlastx of two genomes) from corresponding gbk files")
    
    parser.add_argument("--input_tsv", required=True, help="Path to the input tsv containing the top hits between the two gbk files in outfmt 6.")
    parser.add_argument("--gbk1", required=True, help="Path to the .gbk file 1")
    parser.add_argument("--gbk2", required=True, help="Path to the .gbk file 2")
    parser.add_argument("--output_folder", required=True, help="Path to the output folder where results will be saved.")
    # if all arguments are not provided, print help message
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    retrieve_from_gbk(args.input_tsv, args.gbk1, args.gbk2, args.output_folder)