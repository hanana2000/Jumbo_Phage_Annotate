# 🧬🔎 DIAMOND BLASTP: Find Top Hits for Target Genes


This script batch-runs DIAMOND blastp: for each genome proteome (.faa/.fasta) it builds a DIAMOND database, then queries it with each target gene FASTA to collect the top hits (tabular output). It’s designed for flat folders, simple flags, and reproducible outputs. 

## 💡 What it does

- Builds a DIAMOND database for every proteome file in --genomes_path.
- For each target FASTA in --target_genes, runs diamond blastp against each database.
- Writes one TSV per (target file × genome file) into an output subfolder named with your --prefix, with top five hits for each target seq.
- Optionally filters by query/subject coverage via a constant at the top of the script.

## ⚙️ Requirements
 
- DIAMOND (v2+ recommended)
- Python 3.8+
- Python packages: biopython, pandas (installed but not required by core flow)]

```bash 
# Ubuntu (example)
sudo apt-get install diamond-aligner

# Python deps
pip install biopython pandas

```

## 🧪 Usage 

```bash
python find_top_hits.py \
  --target_genes /path/to/target_genes \
  --genomes_path /path/to/genomes_path \
  --output_folder /path/to/out \
  --prefix PA-

```

If you run find_top_depol_hits.py without required arguments you will get this message: 

```bash
usage: find_top_depol_hits.py [-h] --target_genes TARGET_GENES --genomes_path GENOMES_PATH
                              --output_folder OUTPUT_FOLDER --prefix PREFIX

use DIAMOND to find top hit proteins for each target gene from reference fasta files

options:
  -h, --help            show this help message and exit
  --target_genes TARGET_GENES
                        Path to the folder with subdirectories containing .faa files with reference target
                        genes as queries.
  --genomes_path GENOMES_PATH
                        Path to genomes, all in .faa format with all proteins
  --output_folder OUTPUT_FOLDER
                        Path to the output folder where results will be saved.
  --prefix PREFIX       Prefix for output subdirectory

```

If output folders specified do not exist, they will be created. 

### 📂 Expected inputs

- --target_genes: Flat directory containing .fasta files (each file can have one or many protein sequences). These are the queries.
- --genomes_path: Flat directory containing .faa or .fasta proteomes (all proteins per genome). These are the databases.

The target genes folder is expected to be in flat format with .fasta files contaning target sequences. Example: 

```bash
└── 📁depolymerase_faas_detection
    ├── targets_1.fasta
    ├── depol_detection.sh
    ├── targets_2.fasta
    ├── 📁zone_ids
    ├── targets_3.fasta
    ├── targets_4.fasta
    └── RANDOM.md
```

Only the fasta files will be detected and used as queries. 

The genomes_path is expected to be in flat format as well containing .faa or .fasta files with all proteins from each genome. Example: 

```bash 
└── 📁PA-
    ├── PA-187_Pseudomonas_phage_NEW.faa
    ├── PA-277_Pseudomonas_phage_NEW.faa
    ├── PA-312_Pseudomonas_phage_NEW.faa
    ├── PA-315_Pseudomonas_phage_NEW.faa
    ├── PA-319_Pseudomonas_phage_NEW.faa
    ├── PA-329_Pseudomonas_phage_NEW.faa
    ├── PA-337_Pseudomonas_phage_NEW.faa
    └── PA-711_Pseudomonas_phage_NEW.faa

```

These will be used as the DIAMOND databases to blastp against one at a time. 

The coverage threshold is set at the top of the script. By default it is set to 0, but set to 30 or 80 (or desired coverage) for a list of only close hits. 


## 📊 Output of find_top_depol_hits.py

The output will be stored in the output folder specified.  

you will get output structured as follows in TSVs named like: 

```bash 
└── 📁PA-
    ├── targets_1_PA-711_Pseudomonas_phage_NEW.faa.tsv
    ├── targets_3_PA-711_Pseudomonas_phage_NEW.faa.tsv
    ├── targets_1_PA-312_Pseudomonas_phage_NEW.faa.tsv
    └── targets_4_PA-319_Pseudomonas_phage_NEW.faa.tsv

```

Each file will begin with the query (target gene) file name, then the subject (genome) file name. It will contain all hits for that target-genome pair. 
Each TSV contains BLAST tabular outfmt 6 with extra fields. Example: 

```bash 
APNARBLE_CDS_0211	AGNWHBCZ_CDS_0009	100.0	284	0	0	1	284	1	284	2.0e-175	563.5	285	284
APNARBLE_CDS_0228	AGNWHBCZ_CDS_0027	99.1	456	4	0	1	456	1	456	1.3e-234	920.2	457	456
APNARBLE_CDS_0232	AGNWHBCZ_CDS_0030	98.7	532	7	0	1	532	1	532	1.9e-34	1046.2	532	532
APNARBLE_CDS_0232	AGNWHBCZ_CDS_0035	30.9	521	318	17	18	510	250	756	5.9e-283	194.1	532	787

```

Column meanings

- qseqid → The query sequence ID.
- sseqid → The subject (database) sequence ID.
- pident → Percent identity across the aligned region.
- length → Length of the alignment (number of positions compared).
- mismatch → Number of mismatched positions.
- gapopen → Number of gap openings in the alignment.
- qstart → Start coordinate of the alignment on the query.
- qend → End coordinate of the alignment on the query.
- sstart → Start coordinate of the alignment on the subject.
- send → End coordinate of the alignment on the subject.
- evalue → Expect value (statistical significance of the match). Extremely low = highly significant.
- bitscore → Alignment bit score (higher = better, more reliable alignment).
- qlen → Total length of the query sequence.
- slen → Total length of the subject sequence.


## 🧭 Behavior details

- Skips non-.faa/.fasta proteome files and non-.fasta target files with a [skip] message.
- Builds a fresh .dmnd for each proteome in genomes_path/.
- Writes empty results as no file (empties are removed).
- Uses --max-target-seqs 5 by default. this means only the top 5 hits will be output to the TSV (tweak inside the script if desired).

# 🙋‍♀️ Author/ 📬 Contact

For questions or suggestions, contact: 

Hannah Kapoor
📧 hannahkapoor00@gmail.com 
