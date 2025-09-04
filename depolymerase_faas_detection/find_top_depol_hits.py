from Bio import SeqIO # type: ignore
import argparse
import os
import sys
from pathlib import Path
from collections import defaultdict
import re
import subprocess
import pandas as pd # type: ignore

"""
for each database in the target_genes folder, 
make a db 
iterate through all phages
find top hits 
move on to next db
repeat

"""

QUERY_COVER_THRESH = 0

def get_top_hits(target_genes, genomes_path, output_folder, prefix):
    os.makedirs(output_folder, exist_ok=True)
    databases_dir = os.path.join(output_folder, "databases")
    os.makedirs(databases_dir, exist_ok=True)
    output_results_folder = os.path.join(output_folder, prefix)
    os.makedirs(output_results_folder, exist_ok=True)

    print("\n" + "#@*"* 35 + '\n')
    
    if os.path.isdir(genomes_path):
        # Iterate through genomes and create database to blastp against 
        for db in os.listdir(genomes_path): 
            db_path = os.path.join(genomes_path, db)
            if not db.endswith(".fasta") and not db.endswith(".faa"): 
                print(f"[skip] {db}")
                print("\n" + "#@*"* 35 + '\n')
                continue
            print(f"Processing db: {db}")

            database_path = f"{databases_dir}/{db}_database.dmnd"
            command = ["diamond", "makedb", "--in", db_path, "-d", database_path]
            print(">>", " ".join(command), "\n")
            try:
                result = subprocess.run(command, check=True)
                print(f"{db} database created successfully!\n")
            except subprocess.CalledProcessError as e:
                print("There was an error creating database.\n")
                print(e)
                print("\n" + "#@*"* 35 + '\n')
                continue

            # Iterate through target genes .fasta files found
            target_gene_match(target_genes, output_results_folder, database_path, db)


def target_gene_match(target_genes, output_results_folder, database_path, db): 
    for target_gene_group in os.listdir(target_genes):
        target_genes_path = os.path.join(target_genes, target_gene_group)
        if not target_genes_path.endswith(".fasta"): 
            print(f"[skip] {target_genes_path} is not a .fasta file")
            print("\n" + "#@*"* 35 + '\n')
            continue 
        print(f"Processing target genes: {target_gene_group}")
        output_file = f"{output_results_folder}/{str(target_gene_group).strip(".fasta")}_{db}.tsv"
        command = ["diamond", "blastp", "-d", database_path, "-q", target_genes_path, "-o", output_file, "--query-cover", str(QUERY_COVER_THRESH), "--subject-cover", str(QUERY_COVER_THRESH), "--max-target-seqs", "5", "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
        try:
            result = subprocess.run(command, check=True)
            print(f"{target_gene_group} blastped successfully against {db}!\n")
            if os.path.getsize(output_file)== 0: 
                command = ["rm", output_file]
                try: 
                    subprocess.run(command, check=True)
                except subprocess.CalledProcessError as e: 
                    print(f"Problem deleting file {output_file}")
            else: 
                print(f"Results saved to {output_file}")                    
        except subprocess.CalledProcessError as e:
            print(f"There was an error blastping {target_gene_group} against {database_path}.\n")
            print(e)
            print("\n" + "#@*"* 35 + '\n')
            continue

        print("\n" + "#@*"* 35 + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="use DIAMOND to find top hit proteins for each target gene from reference fasta files")

    parser.add_argument("--target_genes", required=True, help="Path to the folder with subdirectories containing .fasta files with reference target genes as queries.")
    parser.add_argument("--genomes_path", required=True, help="Path to genomes, all in .faa or .fasta format with all proteins")
    parser.add_argument("--output_folder", required=True, help="Path to the output folder where results will be saved.")
    parser.add_argument("--prefix", required=True, help="Prefix for output subdirectory")
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    get_top_hits(args.target_genes, args.genomes_path, args.output_folder, args.prefix)

