# 🧬 Jumbo Phage Annotation and Lysogeny/ Biofilm Degradation Filtering

This repository is for a project annotating ~35 jumbo phages and filtering them for lysogenic markers + biofilm degrading enzymes. 
SPHAE and Interproscan are the main tools used: 

- https://github.com/linsalrob/sphae
- https://github.com/ebi-pf-team/interproscan 

to run the scripts provided, it is expected that SPHAE annotation was performed first. Then Interproscan can be run on all CDSs identified by SPHAE, and results can be combined.

## 🔡 subscripts: 

- The combinesummaries.py script will combine the .txt files output by SPHAE, with information such as lysogeny markers, virulence factors, and CRISPR genes. 
- The interproscan runinterprobatch.py script will take all the predicted amino acid sequences by SPHAE and run them through interproscan to detect any interdomain annotations that might have been missed. 


# 🙋‍♀️ Author/ 📬 Contact

For questions of suggestions, contact: 

Hannah Kapoor
📧 hannahkapoor00@gmail.com 
