# Convert html to tsv format and find top DPO hits from PhageDPO results 

PhageDPO results can be obtained from Phage Galaxy web UI, or by running PhageDPO locally. However, the results are given in html format. This repo contains scripts to convert the html results to tsv format, and then find the top DPO hits (above 85% prediction score) from the tsv files. 

requirements: 
- Python 3.x
- pandas

## Usage and Outputs: 

You must first run html_to_tsv.py to convert the html results to tsv format. Then run top_DPO_hits.py to find the top DPO hits from the tsv files.

example commands: 

```bash 
python3 html_to_tsv.py --PhageDPO_results /home/user/PhageDPO/results/Ph --output_folder ./tsv_results
python3 top_DPO_hits.py --PhageDPO_tsv ./tsv_results/Ph --output_folder /home/user/PhageDPO/top_DPO_hits

```

### Usage of html_to_tsv.py:


If the html_to_tsv.py script is run without any arguments, it will display this help message: 

```bash 
usage: html_to_tsv.py [-h] --PhageDPO_results PHAGEDPO_RESULTS --output_folder OUTPUT_FOLDER

conver html table produced by PhageDPO to tsv format.

options:
  -h, --help            show this help message and exit
  --PhageDPO_results PHAGEDPO_RESULTS
                        Path to the input folder containing PhageDPO results in html format.
  --output_folder OUTPUT_FOLDER
                        Path to the output folder where results will be saved.

```

If the specified output folders do not exist, the scripts will create them.

The input for the --PhageDPO results must be a folder containing subfolders with UNZIPPED html results from PhageDPO. Example: 

```bash 
└── 📁Ph
    └── 📁DPO_Prediction_Ph805
        ├── DPO_Prediction_html.html
    └── 📁DPO_Prediction_Ph818
        ├── DPO_Prediction_html.html
    └── 📁DPO_Prediction_Ph823
        ├── DPO_Prediction_html.html
    └── 📁DPO_Prediction_Ph848
        └── DPO_Prediction_html.html

```


### Output of html_to_tsv.py:

The html_to_tsv.py script will store the tsv files in a subfolder within the specified output folder with the SAME FOLDER NAME as the input folder. 
Feed this subfolder as the input for the top_DPO_hits.py script.

example file structure after running html_to_tsv.py (with tsv_results as output folder): 

```bash
└── 📁tsv_results
    └── 📁Ph
        ├── DPO_Prediction_Ph805.tsv
        ├── DPO_Prediction_Ph818.tsv
        ├── DPO_Prediction_Ph823.tsv
        └── DPO_Prediction_Ph848.tsv

```

you would use "Ph" as the input for the --PhageDPO_tsv argument in top_DPO_hits.py script.


### Usage of top_DPO_hits.py:

If the top_DPO_hits.py script is run without any arguments, it will display this help message: 

```bash 
usage: top_DPO_hits.py [-h] --PhageDPO_tsv PHAGEDPO_TSV --output_folder OUTPUT_FOLDER

find all DPO hits above a specified threshold from PhageDPO tsv results

options:
  -h, --help            show this help message and exit
  --PhageDPO_tsv PHAGEDPO_TSV
                        Path to the input folder containing PhageDPO results in tsv format.
  --output_folder OUTPUT_FOLDER
                        Path to the output folder where results will be saved.

```

If the specified output folders do not exist, the scripts will create them.

the input for the --PhageDPO_tsv must be the subfolder containing the tsv results from html_to_tsv.py script. See the html_to_tsv.py output section above for example file structure.


### Output of top_DPO_hits.py: 

The top_DPO_hits.py script will store the top DPO hits in a subfolder within the specified output folder with the SAME FOLDER NAME as the input folder in tsv format. Example output file: 

```bash
└── 📁top_DPO_hits
    └── 📁Ph
        ├── DPO_Prediction_Ph805_top_DPO.tsv
        ├── DPO_Prediction_Ph818_top_DPO.tsv
        ├── DPO_Prediction_Ph823_top_DPO.tsv
        └── DPO_Prediction_Ph848_top_DPO.tsv
```


Each output tsv file will contain only the DPO hits with prediction score above 85%. However, this threshold can be modified in the script if needed by changing the line at the top of the script: 

```python
PERCENTAGE_THRESHOLD = 85.0

```

This will filter the results to include only those with a prediction score greaten than the specified threshold. 


The tsv files will have two columns of CDS ID and model DPO Prediction (%). Example: 

```bash
CDS	ID	model DPO Prediction (%)
18	AGNWHBCZ_CDS_0018 homing endonuclease Ph805_Mycobacterium_phage	89.0
22	AGNWHBCZ_CDS_0022 baseplate wedge initiator Ph805_Mycobacterium_phage	87.0
25	AGNWHBCZ_CDS_0025 tail sheath Ph805_Mycobacterium_phage	97.0

```

If no DPO hits are found above the threshold, the output tsv file will contain only the header row.


# 🙋‍♀️ Author/ 📬 Contact

For questions of suggestions, contact: 

Hannah Kapoor
📧 hannahkapoor00@gmail.com 