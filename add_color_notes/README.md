# Colour + Notes Annotator for GenBank (.gbk)

Adds two things to every CDS feature in GenBank files:

- /note: ensures each CDS has a note that includes the value(s) from the existing function qualifier (preserving any prior notes).
- /colour: assigns an integer colour code based on a predefined functional category → colour map (Artemis gene color scheme) (compatible with EasyFig).

It can work on:

- A flat directory of .gbk files whose names start with a given prefix, or
- A directory tree (it will search recursively if no flat matches are found).

Requirements

- Python 3.8+
- Biopython

```bash
pip install biopython

```

## Colour map (integer codes)

These match EasyFig’s palette indices:

```bash
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

```

Color choices can be altered at the top of the add_color_notes_recurs.py script. Simply change the value to the desired color in the COLOURS dictionary. New categories can be added to the dictionary as well. 

ARTEMIS color values: 

```bash 

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

```

If a CDS’s function value is not found in the map, no /colour is added for that feature (it still gets /note).


## Usage

if run without any arguments, this message will be displayed: 

```bash 
usage: add_color_notes_recurs.py [-h] --input_folder INPUT_FOLDER --prefix PREFIX
                                 --output_folder OUTPUT_FOLDER

Run interproscan on all phage aa.fasta files of a given species from SPHAE output.

options:
  -h, --help            show this help message and exit
  --input_folder INPUT_FOLDER
                        Path to the input folder containing subdirectories with .gbk files.
  --prefix PREFIX       Prefix to filter the aa.fasta files (e.g. PA-, KA-, Phage-).
  --output_folder OUTPUT_FOLDER
                        Path to the output folder where results will be saved.

```

example command: 

```bash 
python add_colour_notes_recurs.py \
  --input_folder /path/to/gbk_root \
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
    ├── PA-277_Pseudomonas_phage_NEW.gbk
    └── PA-312_Pseudomonas_phage_NEW.gbk

```

## Outputs

New files:
```bash 
OUTPUT/colour_notes_added/PA-319_colour.gbk, etc.

```

Console summary prints each file processed and a list of unique function values encountered (helpful sanity check for missed categories).

