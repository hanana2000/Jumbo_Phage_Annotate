import subprocess
import os
import argparse 
import sys
from Bio import SeqIO # type: ignore
from pathlib import Path

"""
0 white (RGB values: 255 255 255)
1 dark grey (RGB values: 100 100 100)
2 red (RGB values: 255 0 0)
3 green (RGB values: 0 255 0)
4 blue (RGB values: 0 0 255)
5 cyan (RGB values: 0 255 255)
6 magenta (RGB values: 255 0 255)
7 yellow (RGB values: 255 255 0)
8 pale green (RGB values: 152 251 152)
9 light sky blue (RGB values: 135 206 250)
10 orange (RGB values: 255 165 0)
11 brown (RGB values: 200 150 100)
12 pale pink (RGB values: 255 200 200)
13 light grey (RGB values: 170 170 170)
14 black (RGB values: 0 0 0)
15 mid red: (RGB values: 255 63 63)
16 light red (RGB values: 255 127 127)
17 pink (RGB values: 255 191 191)

"""


COLOURS = {
    "head and packaging": 6, # magenta
    "transcription regulation": 7, # yellow
    "other": 9, # light sky blue
    "lysis": 4, # blue
    "unknown function": 5, # cyan
    "connector": 2, # red
    "tail": 3, # green
    "DNA, RNA and nucleotide metabolism": 17, # pink
    "moron, auxiliary metabolic gene and host takeover": 10 # orange    
}


def add_colour_and_notes(input_folder, prefix, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    output_folder_results = os.path.join(output_folder, "colour_notes_added")
    os.makedirs(output_folder_results, exist_ok=True)
    
    print("\n" + "*-" * 40)
    # iterate through genome subdirs 
    root = Path(input_folder)
    gbk_paths = find_gbk_paths(root, prefix)
    if not gbk_paths:
        print(f"no .gbk files found that start with {prefix} and end with .gbk in {input_folder}")
        return
    # add color and notes to gbk file and output in output dir
    for gbk in gbk_paths: 
        parse_gbk(gbk, output_folder_results)
            
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


def parse_gbk(gbk_path, output_folder_results): 
    # read the records from the gbk file if found 
    records, list_funcs = [], []
    file = os.path.basename(gbk_path)
    with open(gbk_path, "r") as f: 
        # use bio to read each CDS for a product and interpro product qualifier
        for record in SeqIO.parse(f, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    if "function" in feature.qualifiers: 
                        ordered_qualifiers, curr_func = order_color_gbk_quals(feature)
                        feature.qualifiers = ordered_qualifiers
                        list_funcs.append(curr_func) # for checking all func types if missed one 
            records.append(record)

    new_file = file.replace("_NEW", "_colour") if "_NEW" in file else file.replace(".gbk", "_colour.gbk")
    SeqIO.write(records, f"{output_folder_results}/{new_file}", "genbank")       

    list_funcs_set = list(set(list_funcs))
    print("")
    print(f"Finished Processing: {file}\n")
    for item in list_funcs_set: print(item)     


def order_color_gbk_quals(feature): 
    """
    orders the .gbk qualifiers, adds a note that is the same 
    as function (or adds onto existing note), adds a colour based on 
    static color dict, and outputs as a dictionary 

    """
    ordered_qualifiers = {}
    curr_func = ""
    for qualifier, value in feature.qualifiers.items():
        if qualifier == "function": # then insert the new "note" qualifier
            ordered_qualifiers[qualifier] = value # use same qualifier and value for function 
            note_vals = feature.qualifiers.get("note", []) # if note already exists
          
            if isinstance(note_vals, str):
                note_vals = [note_vals]
            if isinstance(value, str):
                func_values = [value]
            else:
                func_values = value
            if note_vals: print(f"NOTE: {note_vals}") 
            merged = note_vals + [v for v in func_values if v not in note_vals]
            ordered_qualifiers["note"] = merged # insert the new "note" qualifier that is same as 
            # assign color value based on function
            if func_values[0] in COLOURS: ordered_qualifiers["colour"] = COLOURS[func_values[0]]
            # to list all function types
            curr_func = func_values[0]
        elif qualifier == "note": 
            continue # skip, handled when we hit function qualifier
        else:                                 
            ordered_qualifiers[qualifier] = value
    return ordered_qualifiers, curr_func


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run interproscan on all phage aa.fasta files of a given species from SPHAE output.")
    
    parser.add_argument("--input_folder", required=True, help="Path to the input folder containing subdirectories with .gbk files.")
    parser.add_argument("--prefix", required=True, help="Prefix to filter the aa.fasta files (e.g. PA-, KA-, Phage-).")
    parser.add_argument("--output_folder", required=True, help="Path to the output folder where results will be saved.")
    # if all arguments are not provided, print help message
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    add_colour_and_notes(args.input_folder, args.prefix, args.output_folder)