from Bio import SeqIO
import argparse
import os
import sys
import itertools

MIN_EVALUE = 1e-10
NEXT_BEST_EVALUE = 1e-15


def combine_results(interpro_folder, prefix, SPHAE_folder, output_folder):
    """
    main function: combines all following functions and is called in main

    """
    # Iterate through all subdirectories in the SPHAE folder that begin with the specified prefix
    # Create libraries one phage at a time
    print("\n"+"-"* 60)
    print(" ")  
    output_subfolder = os.path.join(output_folder, f"{prefix}_combined")
    os.makedirs(output_subfolder, exist_ok=True)  # Create output subfolder if it does not exist

    # Iterate through all subdirectories in the SPHAE folder
    for folder in os.listdir(SPHAE_folder):

        # check if the .gbk file is found, if not, move on to the next subdirectory
        gbk_path = os.path.join(SPHAE_folder, folder, f"{folder}.gbk")
        if os.path.exists(gbk_path):
            print(f"Processing GenBank file: {gbk_path}\n")
        else: 
            print(f"\nGenBank file not found: {gbk_path}\n")
            print("#@*"* 45)
            print("#@*"* 45 + '\n')
            continue

        # find the InterProScan results for the current phage
        interpro_path = find_interpro_path(interpro_folder, folder)
        if interpro_path:
            print("-" * 60)
            print(f"Processing InterProScan results: {interpro_path}")
        else:
            print(f"No InterProScan results found for {folder}. Skipping.")
            print("-"* 50)
            continue # Skip to the next folder if no interproscan results found

        # read interpro results for phage we are currently processing
        all_lines_in_phage = read_interpro_tsv(interpro_path)
        print(f"\nFound {len(all_lines_in_phage)} lines in InterProScan results for {folder}.\n")

        # select most specific hits for each locus with the lowest e-value
        top_interpro_hits = select_top_hits(all_lines_in_phage)
        interpro_prod_lib = {line[0].split(":")[1]:line[5] for line in top_interpro_hits}
        interpro_evalue_lib = {line[0].split(":")[1]:line[8] for line in top_interpro_hits}    

        # create a subdirectory for each genome in the output folder if it does not exist 
        # create the output subfolder if it does not exist
        print("-" * 60)
        phage_results_subfolder = os.path.join(output_subfolder, folder)
        os.makedirs(phage_results_subfolder, exist_ok=True)
        # changed to new.gbk for testing and comparison
        new_gbk_path = os.path.join(phage_results_subfolder, f"{folder}_NEW.gbk")

        print(f"\nProcessing folder: {folder}")

        # write to the new_gbk_path file 
        # if the CDS is in the interpro_product_lib then edit the product and source qualifiers
        # glean all SPHAE information by iterating through .gbk ONCE
        all_loci, SPHAE_prod_lib, SPHAE_func_lib, just_SPHAE, interpro_hypo_SPHAE_lib, interpro_annot_SPHAE_lib, yes_interpro, yes_SPHAE, no_SPHAE, fully_hypothetical = write_newgbk(new_gbk_path, interpro_prod_lib, gbk_path)

        # write to a summary text file and note how many CDS were annotated by either SPHAE/Interpro 
        # and how many CDS remain entirely hypothetical 
        write_summary_file(phage_results_subfolder, folder, all_loci, yes_interpro, yes_SPHAE, interpro_hypo_SPHAE_lib, just_SPHAE, interpro_annot_SPHAE_lib, fully_hypothetical)

        # write a .tsv file that has all the top Interpro scan hits and e-values
        write_tophits_tsv(interpro_prod_lib, interpro_evalue_lib, phage_results_subfolder, folder)

        """
        interpro_hypo_SPHAE_lib (CDS that are SPHAE hypothetical, but interpro annotated)
        interpro_annot_SPHAE_lib (CDS that are annotated by both InterProScan and SPHAE)
        just_SPHAE (annotated by just SPHAE)
        fully_hypothetical (not annotated by either SPHAE or InterProScan)

        ALL CDS should fall within one of these categories. 

        """
        
        print("\n")
        print("-"* 120)
        print("#@*"* 40)
        print("#@*"* 40)
        print("-"* 120)
        print("\n")


def select_top_hits(all_lines_in_phage): 
    """
    select the top interpro hit by evalue with the most specific description 

    """
    # remove all ambiguous hits from the unfiltered InterProScan result
    filtered_interpro_hits = remove_ambig_hits(all_lines_in_phage)
    # sort by locus in case out of order (use [-1] so it’s robust to extra colons)
    filtered_interpro_hits.sort(key=lambda r: r[0].split(':')[-1])
    # create top hits list with only hits of the highest e-value 
    current_locus = filtered_interpro_hits[0][0].split(':')[1]
    current_locus_hits = []
    top_interpro_hits = []
    for item in filtered_interpro_hits: 
        locus = item[0].split(':')[1]
        if locus == current_locus: # if part of the same group, append and continue 
            current_locus_hits.append(item)
        else: # if new group, flush the previous group 
            # filter to find the top interpro hit for (old) current  and append to master list of top hits
            top_interpro_hits = next_best_hit(current_locus_hits, top_interpro_hits)
            # reset the current_locus_hits for the next locus
            current_locus_hits = [item]
            # change current locus to new locus (DO NOT REUSE until next iteration)"
            current_locus = locus
    # flush the final group
    if current_locus_hits:
        top_interpro_hits = next_best_hit(current_locus_hits, top_interpro_hits)
    return top_interpro_hits


def next_best_hit(current_locus_hits, top_interpro_hits): 
    """
    finds the top interpro hit of the lowest e value but highest specificity

    """
    best_hit = min(current_locus_hits, key=lambda x: float(x[8].lower()))
    # best_hit = min(current_locus_hits, key=lambda x: float(x[8].lower().split("e")[0]) * 10 ** float(x[8].lower().split("e")[1]))
    current_locus_hits.remove(best_hit)  # Remove the best hit from the filtered list
    if current_locus_hits:
        next_best_hit = min(current_locus_hits, key=lambda x: float(x[8].lower()))
        # check if the next best hit contains the same product as the best hit, but is more specific
        # if so, then change the best hit to the next best hit
        # example: best hit is "GGDEF", while next best hit is "GGDEF diguanylate cyclase"
        if best_hit[5] in next_best_hit[5] and len(best_hit[5]) < len(next_best_hit[5]): 
            if float(next_best_hit[8].lower()) <= NEXT_BEST_EVALUE: 
                best_hit = next_best_hit
    top_interpro_hits.append(best_hit)
    return top_interpro_hits


def remove_ambig_hits(all_lines_in_phage):
    """
    Remove ambiguous hits from the InterProScan results. 

    """
    filtered_interpro_hits = []
    for line in all_lines_in_phage: 
        # do not keep lines that do not have a product or e-value
        if line[8] == "-" or line[5] == '-':
            continue
        e_value = float(line[8].lower())
        if e_value <= MIN_EVALUE: # do not keep hits with e-value greater than global min value
            filtered_interpro_hits.append(line)
    return filtered_interpro_hits


def write_newgbk(new_gbk_path, interpro_prod_lib, gbk_path):
    """
    Writes the records to a new GenBank file and returns 
    libraries of SPHAE/ interpro information 
    while iterating the genbank file only once
    
    """
    all_loci, SPHAE_prod_lib, SPHAE_func_lib, just_SPHAE, interpro_hypo_SPHAE_lib, interpro_annot_SPHAE_lib, yes_interpro = [], {}, {}, [], {}, {}, []
    yes_SPHAE, no_SPHAE, fully_hypothetical = [], [], []
    records = list(SeqIO.parse(gbk_path, "genbank"))
    with open(new_gbk_path, "w") as gbk_file:
        for record in records: 
            for feature in record.features: # iterate through all features in gbk 
                if feature.type == "CDS":
                    locus, SPHAE_product = feature.qualifiers.get("locus_tag", [""])[0], feature.qualifiers["product"][0].lower()
                    all_loci, yes_SPHAE, no_SPHAE, SPHAE_prod_lib, SPHAE_func_lib = populate_SPHAE_libs(all_loci, SPHAE_product, yes_SPHAE, no_SPHAE, SPHAE_prod_lib, SPHAE_func_lib, feature, locus)
                    if locus in interpro_prod_lib: # if there is a InterProScan hit then alter CDS
                        populate_SPHAE_interpro_libs(interpro_prod_lib, yes_interpro, interpro_hypo_SPHAE_lib, interpro_annot_SPHAE_lib, SPHAE_product, locus)
                        # make a new library of ordered qualifiers wiht interpro results added
                        ordered_qualifiers = order_gbk_qualifiers(feature, interpro_prod_lib, locus)
                        feature.qualifiers = ordered_qualifiers
                        # add "InterProScan" as a source for the CDS
                        if feature.qualifiers["source"][0]:
                            feature.qualifiers["source"] = f"{feature.qualifiers['source'][0]}, InterProScan"
                        else: feature.qualifiers["source"] = f"InterProScan"
                    elif SPHAE_product not in ["hypothetical protein", "unknown function"]: 
                        just_SPHAE.append(locus) # if not interpro hit, and not SPHAE hypothetical, append to just SPHAE loci
                    else: # no interpro or SPHAE product
                        fully_hypothetical.append(locus)
            SeqIO.write(record, gbk_file, "genbank") # write the current record to the new .gbk file
    return all_loci, SPHAE_prod_lib, SPHAE_func_lib, just_SPHAE, interpro_hypo_SPHAE_lib, interpro_annot_SPHAE_lib, yes_interpro, yes_SPHAE, no_SPHAE, fully_hypothetical


def populate_SPHAE_libs(all_loci, SPHAE_product, yes_SPHAE, no_SPHAE, SPHAE_prod_lib, SPHAE_func_lib, feature, locus):
    """
    populate the SPHAE libraries 

    """
    all_loci.append(locus)
    if SPHAE_product in ["hypothetical protein", "unknown function"]: no_SPHAE.append(locus)
    else: yes_SPHAE.append(locus)
    SPHAE_prod_lib[locus], SPHAE_func_lib[locus] = feature.qualifiers.get("product", ["hypothetical protein"])[0], feature.qualifiers.get("function", ["unknown function"])[0]
    return all_loci, yes_SPHAE, no_SPHAE, SPHAE_prod_lib, SPHAE_func_lib
    

def populate_SPHAE_interpro_libs(interpro_prod_lib, yes_interpro, interpro_hypo_SPHAE_lib, interpro_annot_SPHAE_lib, SPHAE_product, locus): 
    """
    populate the interpro / SPHAE interpro libraries 

    """
    yes_interpro.append(locus)
    # populate the SPHAE hypothetical and SPHAE annotated interpro hits 
    if SPHAE_product in ["hypothetical protein", "unknown function"]: 
        interpro_hypo_SPHAE_lib[locus] = interpro_prod_lib[locus]
    else: interpro_annot_SPHAE_lib[locus] = interpro_prod_lib[locus]
    return yes_interpro, interpro_hypo_SPHAE_lib, interpro_annot_SPHAE_lib


def order_gbk_qualifiers(feature, interpro_prod_lib, locus): 
    """
    orders the .gbk qualifiers and outputs as a dictionary 

    """
    ordered_qualifiers = {}
    for qualifier, value in feature.qualifiers.items():
        if qualifier == "product": # insert the new "interpro_product" qualifier
            # if it is a hypothetical protein, append interpro product to existing product 
            if value[0].lower() in ["hypothetical protein", "unknown function"]:
                ordered_qualifiers[qualifier] = [value[0] + " | " + interpro_prod_lib[locus]]
            # if it is not a hypothetical protein, just use the original product name
            else: 
                ordered_qualifiers[qualifier] = [value[0]]
            # insert the new "interpro_product" qualifier for ALL CDS with interpro results
            ordered_qualifiers["interpro_product"] = [interpro_prod_lib[locus]]
        else:                                 
            ordered_qualifiers[qualifier] = value
    return ordered_qualifiers


def write_tophits_tsv(interpro_prod_lib, interpro_evalue_lib, phage_results_subfolder, folder): 
    """
    Write a .tsv file containing all top InterProScan hits.

    """
    interpro_hits = []
    for hit in interpro_prod_lib.items():
        interpro_hits.append((hit[0], hit[1], interpro_evalue_lib[hit[0]], "InterproScan"))
    with open(f"{phage_results_subfolder}/{folder}_NEW_top_interpro_hits.tsv", "w") as tsv_file:
        tsv_file.write("Locus\tProduct\tE-value\tSource\n")
        for hit in interpro_hits:
            tsv_file.write(f"{hit[0]}\t{hit[1]}\t{hit[2]}\t{hit[3]}\n")


def write_summary_file(phage_results_subfolder, folder, all_loci, yes_interpro, yes_SPHAE, interpro_hypo_SPHAE_lib, just_SPHAE, interpro_annot_SPHAE_lib, fully_hypothetical):
    """
    Write relevant stats to a summary .txt file

    """
    with open(f"{phage_results_subfolder}/{folder}_NEW_summary.txt", "w") as summary_file:
        summary_file.write(f"Summary for {folder}:\n")
        summary_file.write(f"Total CDS: {len(all_loci)}\n")
        summary_file.write(f"Annotated by InterProScan: {len(yes_interpro)}\n")
        summary_file.write(f"Annotated by SPHAE: {len(yes_SPHAE)}\n")
        summary_file.write(f"Annotated by just InterProScan (NEW): {len(interpro_hypo_SPHAE_lib)}\n")
        summary_file.write(f"Annotated by just SPHAE: {len(just_SPHAE)}\n")
        summary_file.write(f"Annotated by both InterProScan and SPHAE: {len(interpro_annot_SPHAE_lib)}\n")
        summary_file.write(f"Fully hypothetical (no SPHAE or Interpro hits): {len(fully_hypothetical)}\n")
        summary_file.write("\nNewly annotated CDS:\n") 
        for loci in interpro_hypo_SPHAE_lib:
            summary_file.write(f"\t{loci} -> {interpro_hypo_SPHAE_lib[loci]}\n")


def create_interpro_lib(interpro_hits):
    """
    Creates a library of InterProScan  from a list of interpro lines.

    """
    interpro_lib = {}
    for line in interpro_hits:
        locus_tag = line[0].split(":")[1]
        interpro_lib[locus_tag] = line[5]  # Store the InterProScan product for each locus tag
    return interpro_lib


def find_interpro_path(interpro_folder, folder):
    """
    Finds the path to the InterProScan results for the given phage folder.
    
    """
    interpro_path = None
    for filename in os.listdir(interpro_folder):
        if filename.startswith(folder) and os.path.isdir(os.path.join(interpro_folder, filename)): 
            for file in os.listdir(os.path.join(interpro_folder, filename)):
                if file.endswith(".tsv"):
                    interpro_path = os.path.join(interpro_folder, filename, file)
    return interpro_path


def read_interpro_tsv(interpro_path):
    """
    Reads the InterProScan TSV file and returns a list of lines.

    """
    all_lines_in_phage = []
    with open(interpro_path, "r") as interpro_file:
        lines = interpro_file.readlines()
        # Initialize lists for hypothetical and annotated proteins
        # there is no header line in the interpro file
        for line in lines:
            parts = line.strip().split("\t")
            all_lines_in_phage.append(parts)
    return all_lines_in_phage




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine InterProScan results with SPHAE results in a GBK file, with supplementary information contained in summary.txt and top_interpro_hits.tsv.")

    parser.add_argument("--interpro_folder", required= True, help="Path to the input folder containing subdirectories with interproscan .tsv files.")
    parser.add_argument("--prefix", required= True, help="Prefix to filter the aa.fasta files (e.g. PA-, KA-, Phage-).")
    parser.add_argument("--SPHAE_folder", required= True, help="Path to the SPHAE results folder (usually labeled 'final-annotate/' by SPHAE).")
    parser.add_argument("--output_folder", required= True, help="Path to the output folder where results will be saved.")
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 1:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    combine_results(args.interpro_folder, args.prefix, args.SPHAE_folder, args.output_folder)
