# a script to run interproscan from command line on all phage .aa files of a given species 
import subprocess
import os
import argparse 
import sys

def run_interproscan_batch(input_folder, prefix, output_folder, interpro_path):
    # Find all aa.fasta files in the input folder with the specified prefix
    aa_files = []
    for folder in os.listdir(input_folder):
        if folder.startswith(prefix) and folder.endswith("predict") and os.path.isdir(os.path.join(input_folder, folder)):
            Pholder = os.path.join(input_folder, folder)
            # print(f"Processing folder: {Pholder}")
            for filename in os.listdir(Pholder):
                if filename.endswith("aa.fasta") and filename.startswith(prefix):
                    aa_files.append(os.path.join(Pholder, filename))
    print(f"\nFound {len(aa_files)} .aa files with prefix '{prefix}' in '{input_folder}'.\n")
    print(f"Files found:\n")
    # Print the list of found files
    for file in aa_files:
        print(file)
    print("-"*50)
    
    # Run interproscan on each aa.fasta file
    for aa_file in aa_files:
        print(f"\nRunning interproscan on {aa_file}")
        # Call interproscan command here
        command = [interpro_path, "-i", aa_file, "-d", os.path.join(output_folder, os.path.basename(aa_file).replace("aa.fasta", "")), "-tempdir", os.path.join(output_folder, os.path.basename(aa_file).replace("aa.fasta", ""), "TEMP")]
        print(f"\nRunning command: {command}\n")
        command2 = ["rm", "-rf", os.path.join(output_folder, os.path.basename(aa_file).replace("aa.fasta", ""), "TEMP")]
        try:
            result = subprocess.run(command, check=True)
            print("InterProScan ran successfully!")
            subprocess.run(command2, check=True)
            print("Temporary files removed successfully!")
        except subprocess.CalledProcessError as e:
            print("There was an error running InterProScan.")
            print(e)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run interproscan on all phage aa.fasta files of a given species from SPHAE output.")
    
    parser.add_argument("--input_folder", help="Path to the input folder containing subdirectories with aa.fasta files (usually PROCESSING/genome-annotate in SPHAE output).")
    parser.add_argument("--prefix", help="Prefix to filter the aa.fasta files (e.g. PA-, KA-, Phage-).")
    parser.add_argument("--output_folder", help="Path to the output folder where results will be saved.")
    parser.add_argument("--interpro_path", help="path to the interproscan.sh script", default=None)
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 9:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    run_interproscan_batch(args.input_folder, args.prefix, args.output_folder, args.interpro_path)

