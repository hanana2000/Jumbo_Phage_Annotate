import pandas as pd
from pathlib import Path
import argparse
import os
import sys


PERCENTAGE_THRESHOLD = 85.0


def get_top_DPO(PhageDPO_tsv, output_folder):
    os.makedirs(output_folder, exist_ok=True) 
    top_results = f"{output_folder}/{os.path.basename(PhageDPO_tsv)}"
    print(f"Top results will be saved to: {top_results}")
    os.makedirs(top_results, exist_ok=True)

    PhageDPO_tsv = Path(PhageDPO_tsv)
    for item in PhageDPO_tsv.iterdir(): 
        if item.is_file() and item.suffix == ".tsv":
            print(f"now processing file {item}")
            try:
                df = pd.read_csv(item, sep="\t")
                if 'model DPO Prediction (%)' not in df.columns:
                    print(f"'model DPO Prediction (%)' column not found in {item.name}, skipping.")
                    continue
                # get a dataframe of all rows with a score over 80%
                high_conf_df = df[df['model DPO Prediction (%)'] > PERCENTAGE_THRESHOLD]
                # top_df = df.loc[df['model DPO Prediction (%)'].idxmax()]  # get row with highest DPO_Score
                top_tsv_path = f"{top_results}/{item.stem + '_top_DPO.tsv'}"
                high_conf_df.to_csv(top_tsv_path, sep="\t", index=False)  # save as single-row DataFrame
                print(f"Extracted top DPO from {item.name} → {os.path.basename(top_tsv_path)}")

            except Exception as e:
                print(f"Error processing {item.name}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="find all DPO hits above a specified threshold from PhageDPO tsv results")

    parser.add_argument("--PhageDPO_tsv", required= True, help="Path to the input folder containing PhageDPO results in tsv format.")    
    parser.add_argument("--output_folder", required= True, help="Path to the output folder where results will be saved.")
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    get_top_DPO(args.PhageDPO_tsv, args.output_folder)


