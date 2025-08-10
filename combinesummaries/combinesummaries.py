import os 
import argparse
import sys


def combine_summaries(input_folder, prefix):
    """
    Combine summary files
    if the output folder already exists, it will not be created again
    if it does not exist, it will be created

    """
    output_file = os.path.join(input_folder, "allsummaries.txt")
    with open(output_file, "w") as outfile:
    # Iterate through each folder in the input directory
        for folder in os.listdir(input_folder):
            # if the folder starts with PA and is a directory
            if folder.startswith(prefix) and os.path.isdir(os.path.join(input_folder, folder)):
                final_annotate_folder = os.path.join(input_folder, folder)
                for filename in os.listdir(final_annotate_folder):
                    if filename.endswith("_summary.txt"):
                        with open(os.path.join(final_annotate_folder, filename), "r") as infile:
                            print(f"Combining {filename} from {final_annotate_folder} into {input_folder}")
                            #test by writing file contents to the terminal
                            # print(infile.read())
                            print("\n")
                            outfile.write(f"Contents of {filename} from {final_annotate_folder}:\n")
                            # Write the contents of the file to the output file
                            outfile.write(infile.read() + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine summary files from SPHAE output in final_annotate folders. this will output a single txt file titled 'allsummaries.txt' in the input folder.")
    parser.add_argument("--input_folder", help="Path to the input folder containing subdirectories with summary files (usually titled 'final_annotate' by SPHAE).")
    parser.add_argument("--prefix", help="Prefix for the subdirectories to include (e.g. PA-, KA-, Phage-).")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    combine_summaries(args.input_folder, args.prefix)
