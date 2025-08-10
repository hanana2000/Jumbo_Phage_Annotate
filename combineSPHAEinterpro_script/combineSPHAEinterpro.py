from Bio import SeqIO
import argparse
import os
import sys
import itertools

def combine_results(interpro_folder, prefix, SPHAE_folder, output_folder):
    # Iterate through all subdirectories in the SPHAE folder that begin with the specified prefix
    # Create libraries one phage at a time
    print("\n"+"-"* 60)
    print(" ")  
    output_subfolder = os.path.join(output_folder, f"{prefix}_results_combined")
    os.makedirs(output_subfolder, exist_ok=True)  # Create output subfolder if it does not exist

    # Iterate through all subdirectories in the SPHAE folder
    for folder in os.listdir(SPHAE_folder):

        # check if the .gbk file is found, if not, move on to the next subdirectory
        gbk_path = os.path.join(SPHAE_folder, folder, f"{folder}.gbk")
        if os.path.exists(gbk_path):
            print(f"Processing GenBank file: {gbk_path}")
        else: 
            print(f"\nGenBank file not found: {gbk_path}")
            print("\n")
            print("#@*"* 45)
            print("#@*"* 45 + '\n')
            continue

        # create SPHAE locus_tag pair dictionaries for products and functions
        # each CDS will have one SPHAE entry in each dictionary
        locus_function = create_gbk_lib("function", folder, SPHAE_folder, prefix, gbk_path)
        locus_product = create_gbk_lib("product", folder, SPHAE_folder, prefix, gbk_path)
        # create lists of loci that are and aren't annotated by SPHAE
        all_loci = [loci[0] for loci in locus_product.items()]
        no_SPHAE = [loci[0] for loci in locus_product.items() if loci[1] in ["hypothetical protein", "unknown function"]]
        yes_SPHAE = [loci[0] for loci in locus_product.items() if loci[1] not in ["hypothetical protein", "unknown function"]]
        fully_hypothetical = []
        print(f"\nLoaded {len(all_loci)} products from GenBank file.\n")

        # print product and function information
        print(f"SPHAE products and functions for {folder} (top):")
        for locus in all_loci[:5]:
            print(f"Product: {locus} -> {locus_product[locus]} \t Function: {locus_function[locus]}")
        print("\n")

        # Print all the locus tags 
        print("Locus tags: (top)")
        print_all_byfour(all_loci[:10])

        # # print all loci not annotated by SPHAE
        # print_filter_yesno(all_loci, yes_SPHAE, no_SPHAE, "SPHAE")

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
        print(f"\nFound {len(all_lines_in_phage)} lines in InterProScan results for {folder}.")

        # Filter the InterProScan results based on the locus tags from the SPHAE library
        # select most specific hits for each locus with the lowest e-value
        interpro_sort = interpro_sort_func(all_lines_in_phage, locus_product, all_loci)
        top_interpro_hits, no_interpro, yes_interpro = interpro_sort[0], interpro_sort[1], interpro_sort[2]
        interpro_prod_lib = create_interpro_lib(top_interpro_hits)  # Store the InterProScan product for each locus tag
        # Create a dictionary to store locus:e-value pairs for each locus tag  
        interpro_evalue_lib = {line[0].split(":")[1]:line[8] for line in top_interpro_hits}

        # # print loci annotated by interpro and loci not annotated by interpro
        # print_filter_yesno(all_loci, yes_interpro, no_interpro, "InterProScan")

        # Find all InterproScan hits that were hypothetical with unknown function in the SPHAE library
        # and all InterproScan hits that are also annotated by SPHAE
        hypo_annot = split_interpro_results(locus_product, locus_function, top_interpro_hits)
        interpro_hypo_SPHAE, interpro_annot_SPHAE = hypo_annot[0], hypo_annot[1]
        interpro_hypo_SPHAE_lib, interpro_annot_SPHAE_lib = create_interpro_lib(interpro_hypo_SPHAE), create_interpro_lib(interpro_annot_SPHAE)
        
        # get lists off Loci that are not annotated by either InterProScan or SPHAE
        # or were just annotated by SPHAE, can also print out results in four column table
        fully_hypothetical, just_SPHAE = fully_hypo_or_justSPHAE(all_loci, interpro_prod_lib, interpro_hypo_SPHAE_lib, interpro_annot_SPHAE_lib, no_SPHAE, yes_SPHAE, locus_product, locus_function)  # List to store fully hypothetical loci
        # print_fully_hypo_or_justSPHAE(all_loci, interpro_prod_lib, interpro_hypo_SPHAE_lib, interpro_annot_SPHAE_lib, no_SPHAE, yes_SPHAE, locus_product, locus_function)


        """
        AT THIS POINT we now have: 

        interpro_hypo_SPHAE_lib (CDS that are SPHAE hypothetical, but interpro annotated)
        interpro_annot_SPHAE_lib (CDS that are annotated by both InterProScan and SPHAE)
        just_SPHAE (annotated by just SPHAE)
        fully_hypothetical (not annotated by either SPHAE or InterProScan)

        ALL CDS should fall within one of these categories. 

        """

        # create a subdirectory for each genome in the output folder if it does not exist 
        # create the output subfolder if it does not exist
        print("-" * 60)
        phage_results_subfolder = os.path.join(output_subfolder, folder)
        os.makedirs(phage_results_subfolder, exist_ok=True)
        new_gbk_path = os.path.join(phage_results_subfolder, f"{folder}.gbk")
        
        # write to a summary text file and note how many CDS were annotated by either SPHAE/Interpro 
        # and how many CDS remain entirely hypothetical 
        write_summary_file(phage_results_subfolder, folder, all_loci, yes_interpro, yes_SPHAE, interpro_hypo_SPHAE_lib, just_SPHAE, interpro_annot_SPHAE_lib, fully_hypothetical)

        # write a .tsv file that has all the top Interpro scan hits and e-values
        write_tophits_tsv(interpro_prod_lib, interpro_evalue_lib, phage_results_subfolder, folder)

        # write to the new_gbk_path file 
        # if the CDS is in the interpro_product_lib then edit the product and source qualifiers 
        write_newgbk(new_gbk_path, interpro_prod_lib, gbk_path)


        print("\n")
        print("-"* 120)
        print("#@*"* 40)
        print("#@*"* 40)
        print("-"* 120)
        print("\n")


def write_newgbk(new_gbk_path, interpro_prod_lib, gbk_path):
    """
    Writes the records to a new GenBank file.
    
    """
    records = list(SeqIO.parse(gbk_path, "genbank"))
    with open(new_gbk_path, "w") as gbk_file:
        for record in records: 
            for feature in record.features: # iterate through all features in gbk 
                if feature.type == "CDS":
                    locus = feature.qualifiers.get("locus_tag", [""])[0]
                    if locus in interpro_prod_lib: # if there is a InterProScan hit then alter CDS
                        # make a new library of ordered qualifiers 
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
                        feature.qualifiers = ordered_qualifiers
                        # add "InterProScan" as a source for the CDS
                        if feature.qualifiers["source"][0]:
                            feature.qualifiers["source"] = f"{feature.qualifiers['source'][0]}, InterProScan"
                        else:
                            feature.qualifiers["source"] = f"InterProScan"
            SeqIO.write(record, gbk_file, "genbank") # write the current record to the new .gbk file



def write_tophits_tsv(interpro_prod_lib, interpro_evalue_lib, phage_results_subfolder, folder): 
    """
    Write a .tsv file containing all top InterProScan hits.

    """
    interpro_hits = []
    for hit in interpro_prod_lib.items():
        interpro_hits.append((hit[0], hit[1], interpro_evalue_lib[hit[0]], "InterproScan"))
    with open(f"{phage_results_subfolder}/{folder}_top_interpro_hits.tsv", "w") as tsv_file:
        tsv_file.write("Locus\tProduct\tE-value\tSource\n")
        for hit in interpro_hits:
            tsv_file.write(f"{hit[0]}\t{hit[1]}\t{hit[2]}\t{hit[3]}\n")


def write_summary_file(phage_results_subfolder, folder, all_loci, yes_interpro, yes_SPHAE, interpro_hypo_SPHAE_lib, just_SPHAE, interpro_annot_SPHAE_lib, fully_hypothetical):
    """
    Write relevant stats to a summary .txt file

    """
    with open(f"{phage_results_subfolder}/{folder}_summary.txt", "w") as summary_file:
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


def fully_hypo_or_justSPHAE(all_loci, interpro_prod_lib, interpro_hypo_SPHAE_lib, interpro_annot_SPHAE_lib, no_SPHAE, yes_SPHAE, locus_product, locus_function):  # List to store fully hypothetical loci
    """
    Create two lists: fully_hypothetical and just_SPHAE
    these contain locus keys of CDS that are not annotated by either SPHAE or Interpro
    and CDS that were just annotated by SPHAE

    """
    fully_hypothetical = []  # List to store fully hypothetical loci
    just_SPHAE = []  # List to store loci that are only annotated by SPHAE
    for locus in all_loci:
        if locus not in interpro_prod_lib:
            if locus in no_SPHAE: #No SPHAE or Interpro hit
                fully_hypothetical.append(locus)
            elif locus in yes_SPHAE: #Just SPHAE
                just_SPHAE.append(locus)
    return fully_hypothetical, just_SPHAE


def print_fully_hypo_or_justSPHAE(all_loci, interpro_prod_lib, interpro_hypo_SPHAE_lib, interpro_annot_SPHAE_lib, no_SPHAE, yes_SPHAE, locus_product, locus_function):  # List to store fully hypothetical loci
    """
    Print out a four column table listing just SPHAE hits, just Interpro hits. 
    Interpro and SPHAE hits, or CDS with no hits

    """
    print("\nJust SPHAE \t Interpro no SPHAE \t Interpro and SPHAE \t No SPHAE or Interpro hit")
    print("-"* 100)
    fully_hypothetical = []  # List to store fully hypothetical loci
    just_SPHAE = []  # List to store loci that are only annotated by SPHAE
    for locus in all_loci:
        if locus in interpro_prod_lib:
            if locus in interpro_hypo_SPHAE_lib: #Interpro no SPHAE
                print (f"\t \t {locus}")
                print(f"\t \t Interpro result: {interpro_hypo_SPHAE_lib[locus]}")
                print(f"\t \t SPHAE product: {locus_product[locus]}, function: {locus_function[locus]}")
            elif locus in interpro_annot_SPHAE_lib: #Interpro and SPHAE
                print (f"\t \t \t \t \t {locus}")
                print(f"\t \t \t \t \t Interpro result: {interpro_annot_SPHAE_lib[locus]}")
                print(f"\t \t \t \t \t SPHAE product: {locus_product[locus]}, function: {locus_function[locus]}")
        elif locus not in interpro_prod_lib:
            if locus in no_SPHAE: #No SPHAE or Interpro hit
                print(f"\t \t \t \t \t \t \t \t {locus}")
                fully_hypothetical.append(locus)
                print(f"\t \t \t \t \t \t \t \t SPHAE product: {locus_product[locus]}, function: {locus_function[locus]}")
            elif locus in yes_SPHAE: #Just SPHAE
                print(locus)
                just_SPHAE.append(locus)
                print(f"SPHAE product: {locus_product[locus]}, function: {locus_function[locus]}")
    return fully_hypothetical, just_SPHAE


def create_interpro_lib(interpro_hits):
    """
    Creates a library of InterProScan hits.

    """
    interpro_lib = {}
    for line in interpro_hits:
        locus_tag = line[0].split(":")[1]
        interpro_lib[locus_tag] = line[5]  # Store the InterProScan product for each locus tag
    return interpro_lib


def print_filter_yesno(all_loci, yes_annot, no_annot, keyword):
    """
    Prints all locus tags that are not annotated by specified program
    and loci that are in 2 columns.

    """
    print(f"Loci annotated by {keyword}:\t Loci not annotated by {keyword}:\n" + '-'* 60)
    for locus in all_loci:
        if locus in no_annot:
            print(f" \t \t \t \t {locus}")
        elif locus in yes_annot: 
            print(f"{locus}")
    print("\n")


def print_all_byfour(all_loci): 
    """
    Prints all locus tags in groups of four.

    """
    four_count = -1
    line = ""
    for locus in all_loci: 
        four_count += 1
        if four_count % 4 == 0:
            print(line)
            line = locus + '\t'
        else:
            line += locus + '\t'
        if four_count == len(all_loci) - 1:
            print(line)
    print("\n")


def print_interproSPHAE_lib(library, locus_product, locus_function, keyword): 
    """
    prints out interpro annotated or interpro hypothetical libraries.

    """
    identifier = "hypothetical" if keyword == "unannotated" else "annotated"
    print(f"InterProScan hits that were {identifier} in the SPHAE library:")
    print("-"* 50)
    for line in library:
        print(f"{line[0].split(':')[1]} with interpro product: {line[5]}, e-value: {line[8]}")
        print(f"SPHAE product: {locus_product.get(line[0].split(':')[1])}, function: {locus_function.get(line[0].split(':')[1])}")
        print("+" * 10)
    print("\n")


def split_interpro_results(locus_product, locus_function, top_interpro_hits):
    """
    Find all InterproScan hits that were labeled hypothetical protein 
    in the SPHAE library vs all InterproScan hits that were annotated by SPHAE

    """
    interpro_hypo_SPHAE = []
    interpro_annot_SPHAE = []
    for line in top_interpro_hits:
        locus_tag = line[0].split(":")[1]
        if locus_product[locus_tag].lower() in ["hypothetical protein", "unknown function"]:
            interpro_hypo_SPHAE.append(line) # Hypothetical SPHAE hits labeled by interpro
        else: 
            interpro_annot_SPHAE.append(line) # Annotated SPHAE hits labeled by interpro
    return interpro_hypo_SPHAE, interpro_annot_SPHAE


def interpro_sort_func(all_lines_in_phage, locus_product, all_loci):
    """
    Sorts InterProScan hits into top hits and no hits based on e-value.

    """
    top_interpro_hits = [] # List to store the best InterProScan hits
    no_interpro = []  # List to store locus tags with no InterProScan hits
    yes_interpro = [] # List to store locus tags with InterProScan hits
    # for locus, product in itertools.islice(locus_product.items(), 3):
    for locus in all_loci:
        # remove all ambiguous hits from the unfiltered InterProScan results
        # this will remove hits with e-value greater than 1e-10 or no product
        filtered_interpro = remove_ambiguous_hits(all_lines_in_phage, locus)
        # Find the hit with the lowest e-value for each locus tag
        if filtered_interpro:
            best_hit = best_interpro_hit(filtered_interpro, locus)
            top_interpro_hits.append(best_hit)
            yes_interpro.append(locus)
        else: 
            no_interpro.append(locus)
    return top_interpro_hits, no_interpro, yes_interpro


def remove_ambiguous_hits(all_lines_in_phage,locus):
    """
    Removes ambiguous hits from the unfiltered InterProScan results.
    (i.e. hits with no e value, hits with no product)

    """
    # Convert e-value to float for comparison
    # e_value = float(line[8]) will not work if e-value is a string like "1.96E-18"
    filtered_interpro = [line for line in all_lines_in_phage if line[0].split(":")[1] == locus]
    to_remove = []
    for line in filtered_interpro:
        # Check if the e-value is ambiguous or no scientific notation (e.g. 8.543)
        # also check if the product is "-" (empty)
        if line[8] == "-" or 'E' not in line[8] or line[5] == '-':
            to_remove.append(line)
            continue                
    # Remove the lines that were marked for removal
    # This is done to work around modifying the list while iterating over it, 
    # avoiding unexpected behavior
    for line in to_remove:
        filtered_interpro.remove(line)
    to_remove = []  # Reset the list
    # Filter out hits with e-value greater than 1e-10
    for line in filtered_interpro:
        # Convert e-value to float for comparison
        e_value = (float(line[8].split("E")[0]) * (10 ** float(line[8].split("E")[1])))
        if e_value > 1e-10: # Filter out hits with e-value greater than 1e-10
            to_remove.append(line)
    for line in to_remove:
        filtered_interpro.remove(line)
    return filtered_interpro


def best_interpro_hit(filtered_interpro, locus): 
    """
    Finds the most specific hit with the lowest e-value for each locus tag.

    """
    best_hit = min(filtered_interpro, key=lambda x: float(x[8].split("E")[0]) * (10 ** float(x[8].split("E")[1])))
    filtered_interpro.remove(best_hit)  # Remove the best hit from the filtered list
    if filtered_interpro:
        next_best_hit = min(filtered_interpro, key=lambda x: float(x[8].split("E")[0]) * (10 ** float(x[8].split("E")[1])))
        # check if the next best hit contains the same product as the best hit, but is more specific
        # if so, then change the best hit to the next best hit
        # example: best hit is "GGDEF", while next best hit is "GGDEF diguanylate cyclase"
        if best_hit[5] in next_best_hit[5] and len(best_hit[5]) < len(next_best_hit[5]):
            best_hit = next_best_hit
    return best_hit


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


def create_gbk_lib(keyword, folder, SPHAE_folder, prefix, gbk_path):
    """
    Creates a library from the GenBank file for the specified phage folder.
    The library maps locus tags to products or functions (or any other qualifier) based on the keyword.
    
    """
    fallback = "hypothetical protein" if keyword == "product" else "unknown function"
    library = {}
    if folder.startswith(prefix) and os.path.isdir(os.path.join(SPHAE_folder, folder)): 
        # Read the GenBank file and create a lookup dictionary for products
        # and functions based on locus_tag
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    item = feature.qualifiers.get(f"{keyword}", f"{fallback}")[0]
                    locus_tag_value = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                    library[locus_tag_value] = item.lower()
    return library



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split interproscan results into two folders -- one for SPHAE labeled hypothetical proteins and one for SPHAE annotated proteins.")
    
    parser.add_argument("--interpro_folder", help="Path to the input folder containing subdirectories with interproscan .tsv files.")
    parser.add_argument("--prefix", help="Prefix to filter the aa.fasta files (e.g. PA-, KA-, Phage-).")
    parser.add_argument("--SPHAE_folder", help="Path to the SPHAE results folder (usually labeled 'final-annotate/' by SPHAE).")
    parser.add_argument("--output_folder", help="Path to the output folder where results will be saved.")
    # if all 4 arguments are not provided, print help message
    if len(sys.argv) < 7:
        parser.print_help(sys.stderr)
        print("\n")
        sys.exit(1)
    args = parser.parse_args()

    combine_results(args.interpro_folder, args.prefix, args.SPHAE_folder, args.output_folder)
