Jumbo Phage Annotation -- combine summaries script

This script expects that you have already run SPHAE annotation on all your phages, and that you have a consistent naming prefix for your files (e.g. all files begin with PA- or Phage-). 
see https://github.com/linsalrob/sphae for more information. 

requirements: 
- python 3.x


## 🚀 Usage

You can run the script from the command line as follows:

```bash
python combinesummaries.py --input_folder /path/to/final-annotate/ --prefix PA-

```
This script expects to be pointed to the SPHAE output file, usually titled "final_annotate". This should have individual subdirectories for each annotated genome within it, and in each subdirectory there should be a .gbk, .fasta, .png, .functions, and .txt file. The .txt file is the summary file that this script is looking for. 

```bash 
final-annotate/
├── PA-187_Pseudomonas_phage/
│ └── PA-187_Pseudomonas_phage_summary.txt
├── PA-329_Pseudomonas_phage/
│ └── PA-329_Pseudomonas_phage_summary.txt

```

if you run the command without any arguments, you will get a list of required/ accepted flags: 

```bash 

User@MSI:~/Folder$ python3 combinesummaries.py
usage: combinesummaries.py [-h] [--input_folder INPUT_FOLDER] [--prefix PREFIX]

Combine summary files from SPHAE output in final_annotate folders. this will output a single txt file titled 'allsummaries.txt' in the
input folder.

options:
  -h, --help            show this help message and exit
  --input_folder INPUT_FOLDER
                        Path to the input folder containing subdirectories with summary files (usually titled 'final_annotate' by
                        SPHAE).
  --prefix PREFIX       Prefix for the subdirectories to include (e.g. PA-, KA-, Phage-).
```

Output will be a single .txt file titled "allsummaries.txt" in the "final annotate/" folder, contaning concatenated results of all "_summary.txt" files.

For questions of suggestions, contact: 

Hannah Kapoor
📧 hannahkapoor00@gmail.com 


