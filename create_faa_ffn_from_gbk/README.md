# 🧬 Generate .faa or .ffn from .gbk

Creates a .faa and/or a .ffn file for every .gbk file in directory.

It can work on:

- A flat directory of .gbk files whose names start with a given prefix, or
- A directory tree (it will search recursively if no flat matches are found).

Requirements

- Python 3.8+
- Biopython

```bash
pip install biopython

```

## 🧪 Usage

if run without any arguments, this message will be displayed: 

```bash 
usage: create_faa_from_gbk.py [-h] --genomes_folder GENOMES_FOLDER --prefix PREFIX
                              --output_folder OUTPUT_FOLDER

Create .faa files for all .gbk files found in passed directory.

options:
  -h, --help            show this help message and exit
  --genomes_folder GENOMES_FOLDER
                        Path to the input folder containing subdirectories with .gbk files.
  --prefix PREFIX       Prefix for the output file subdirectory
  --output_folder OUTPUT_FOLDER
                        Path to the output folder where results will be saved.

```

example command: 

```bash 
python create_faa_from_gbk.py \
  --genomes_folder /path/to/gbk_root \
  --prefix PA- \
  --output_folder /path/to/output

```

If the gbk files are nested in subdirs, then the program will use recursive search. 

example recursive input:

```bash
└── 📁PA-_combined
    └── 📁PA-187_Pseudomonas_phage
        ├── PA-187_Pseudomonas_phage_NEW_summary.txt
        ├── PA-187_Pseudomonas_phage_NEW_top_interpro_hits.tsv
        ├── PA-187_Pseudomonas_phage_NEW.gbk
    └── 📁PA-277_Pseudomonas_phage
        ├── PA-277_Pseudomonas_phage_NEW_summary.txt
        ├── PA-277_Pseudomonas_phage_NEW_top_interpro_hits.tsv
        ├── PA-277_Pseudomonas_phage_NEW.gbk
    └── 📁PA-312_Pseudomonas_phage
        ├── PA-312_Pseudomonas_phage_NEW_summary.txt
        ├── PA-312_Pseudomonas_phage_NEW_top_interpro_hits.tsv
        ├── PA-312_Pseudomonas_phage_NEW.gbk

```

the program can also take flat input: 

```bash
└── 📁PA-_combined
    ├── PA-187_Pseudomonas_phage_NEW.gbk
    ├── TEST_file.txt
    ├── PA-277_Pseudomonas_phage_NEW.gbk
    └── PA-312_Pseudomonas_phage_NEW.gbk

```

the inputs are the same for the ffn version. 

## 📊 Outputs

New files:
```bash 
OUTPUT/Prefix/PA-319.faa, etc.

```

Console summary prints each file being processed in order.

for the .fna files, the output will be in .faa format, and the ffn files will be in .fasta format. 
you can rename the files to have the .fasta extension if desired.

# 🙋‍♀️ Author/ 📬 Contact

For questions of suggestions, contact: 

Hannah Kapoor
📧 hannahkapoor00@gmail.com 
