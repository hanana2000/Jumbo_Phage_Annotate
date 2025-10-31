import pandas as pd
from pathlib import Path
import argparse
import os
import sys


def convert_html_to_tsv(PhageDPO_results, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    tsv_results = f"{output_folder}/{os.path.basename(PhageDPO_results)}"
    print(f"TSV results will be saved to: {tsv_results}")
    os.makedirs(tsv_results, exist_ok=True)
    PhageDPO_results = Path(PhageDPO_results)
    for item  in PhageDPO_results.iterdir():
        if item.is_dir():
            print(f"now processing folder {item}")
            html_files = list(item.glob("*.html")) 
        
            for html_file in html_files:
                try:
                    # PhageDPO results are usually one table per HTML file
                    tables = pd.read_html(html_file)
                    if not tables:
                        print(f"No tables found in {html_file.name}")
                        continue

                    df = tables[0]  # usually the first table is the prediction summary
                    tsv_path = f"{tsv_results}/{item.stem + '.tsv'}"
                    df.to_csv(tsv_path, sep="\t", index=False)
                    print(f"Converted {html_file.name} → {os.path.basename(tsv_path)}")

                except Exception as e:
                    print(f"Error reading {html_file.name}: {e}")

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="conver html table produced by PhageDPO to tsv format.")

    parser.add_argument("--PhageDPO_results", required= True, help="Path to the input folder containing PhageDPO results in html format.")    
    parser.add_argument("--output_folder", required= True, help="Path to the output folder where results will be saved.")
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    convert_html_to_tsv(args.PhageDPO_results, args.output_folder)