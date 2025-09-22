# tsv_retreive_gbk.py

If you have run a tblastx between two phage genomes (.fasta nucleotide seqs), and have generated a tabular output file (outfmt 6), you can use this script to retrieve the annotation of the hits from the corresponding gbk files.

note: this script assumes that there is only ONE CDS feature that corresponds to the hit coordinates. If there are multiple overlapping CDS features, it will only retrieve the first one. 

## Usage: 

if you run the script without any arguments, it will display the help message:

```bash
python3 tsv_gbk_retreive.py

usage: tsv_gbk_retreive.py [-h] --input_tsv INPUT_TSV --gbk1 GBK1 --gbk2 GBK2 --output_folder OUTPUT_FOLDER

retreive all the hits from a tsv (tBlastx of two genomes) from corresponding gbk files

options:
  -h, --help            show this help message and exit
  --input_tsv INPUT_TSV
                        Path to the input tsv containing the top hits between the two gbk files in outfmt 6.
  --gbk1 GBK1           Path to the .gbk file 1
  --gbk2 GBK2           Path to the .gbk file 2
  --output_folder OUTPUT_FOLDER
                        Path to the output folder where results will be saved.

```

