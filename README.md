# 🧬 Jumbo Phage Annotation and Lysogeny/ Biofilm Degradation Filtering

This repository is for a project annotating ~35 jumbo phages and filtering them for lysogenic markers + biofilm degrading enzymes. 
SPHAE and Interproscan are the main tools used: 

- https://github.com/linsalrob/sphae
- https://github.com/ebi-pf-team/interproscan 

to run the scripts provided, it is expected that SPHAE annotation was performed first. Then Interproscan can be run on all CDSs identified by SPHAE, and results can be combined.

## 🔡 subscripts: 

- The combinesummaries.py script will combine the .txt files output by SPHAE, with information such as lysogeny markers, virulence factors, and CRISPR genes. 
- The runinterprobatch.py script will take all the SPHAE predicted amino acid sequences and run them through InterProScan to detect any interdomain annotations that might have been missed. 
- The combineSPHAEinterpro.py script will take the results from SPHAE and InterProScan (after running the runinterprobatch.py script) and combine them into a single .gbk output file. It will also generate a summary file and a top_interpro_hits.tsv file for each genome.

example pipeline: 
```bash
# Step 1: Run SPHAE annotation: 
(sphae) user@MSI:~/folder$ sphae annotate --genome /path/to/PA_genomes/ --output example -k

# Step 2: Combine txt summaries (optional): 
(sphae) user@MSI:~/folder$ python combinesummaries.py --input_folder /path/to/PA_genomes/example/final-annotate/ --prefix PA-

# Step 3: Run InterProScan on SPHAE predicted proteins: 
(sphae) user@MSI:~/folder$ python3 runinterprobatch.py --input /path/to/PA_genomes/example/PROCESSING/genome-annotate/ --prefix PA- --output_folder /path/to/InterproScan_Results/ --interpro_path /path/to/interproscan.sh

# Step 4: Combine SPHAE and InterProScan results: 
(sphae) user@MSI:~/folder$ python3 combineSPHAEinterpro.py --interpro_folder /path/to/InterproScan_Results/PA- --prefix PA- --SPHAE_folder /path/to/PA_genomes/example/final-annotate/ --output_folder /path/to/InterproScan_Results/

```
- The find_top_depol_hits.py script batch-runs DIAMOND blastp. For each genome proteome (.faa/.fasta) it builds a DIAMOND database, then queries it with each target gene FASTA to collect the top hits (tabular output).
- The create_faa_from_gbk.py script creates a .faa file for every .gbk file in the passed directory. it can take a flat directory or nested directory. 
- The add_colour_notes.py script adds two things to every CDS feature in GenBank files:
    - /note: ensures each CDS has a note that includes the value(s) from the existing function qualifier (preserving any prior notes).
    - /colour: assigns an integer colour code based on a predefined functional category → colour map (Artemis gene color scheme) (compatible with EasyFig).
- The tsv_gbk_retrieve.py script will retrieve the annotation of the hits from the corresponding gbk files if you have run a tblastx between two phage genomes (.fasta nucleotide seqs), and have generated a tsv (outfmt 6).

## 🛠 PhageDPO Post-Processing Scripts:
PhageDPO is a tool for predicting depolymerase enzymes in phage genomes. This repository contains two scripts to help process the results from PhageDPO:
- The html_to_tsv.py script will convert the html results produced by PhageDPO to tsv format.
- The top_DPO_hits.py script will find the top DPO hits from the tsv files produced by html_to_tsv.py script.

See each subdirectory README.md for specific info on each script. 

# 🙋‍♀️ Author/ 📬 Contact

For questions of suggestions, contact: 

Hannah Kapoor
📧 hannahkapoor00@gmail.com 
