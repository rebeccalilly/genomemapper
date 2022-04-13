from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from progress.bar import Bar
import os
import re

if not (os.path.exists("results.xml")):
    sequence_file = input("Enter the name of your sequence file in FASTA format:")
    sequence_data = open(sequence_file).read()
    bar = Bar("Processing sequence data... ")
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)
    bar.finish()

    with open("results.xml","w") as save_file:
        blast_results = result_handle.read()
        save_file.write(blast_results)

results = open("results.xml","r")
records = list(NCBIXML.parse(results)) #parse results (xml file) into a list
one_query = records[0] #from a list to a record class that holds all BLAST output

accession_numbers = [] #create empty list to store the accession_numbers of each organism

for i in range(len(one_query.alignments)): #loop through every query to obtain the accession_number
    one_hit = one_query.alignments[i]
    accession_numbers.append(one_hit.accession)


#parsing GenBank data
from Bio import Entrez
Entrez.email = input("Enter your email address (so that NCBI can contact you if there's a problem):")
#bar = Bar("Parsing location for each accession numbers from GenBank... ", max = len(accession_numbers))

location_list = []
for i in accession_numbers:
    handle = Entrez.efetch(db = 'nucleotide', rettype = 'gb', retmode = 'text', id = i)
    record = handle.read()
    info = record.split("\n")
    location_line = [i for i in info if i.startswith('  JOURNAL   Submitted') or i.startswith('FEATURES') or i.startswith('COMMENT')]
    location_line_indices = []
    for j in location_line:
        index_number = info.index(j)
        location_line_indices.append(index_number)

    index_1 = location_line_indices[0]
    index_2 = location_line_indices[-1]

    location_info = info[index_1:index_2]
    first_line = location_info[0]
    first_line_split = first_line.split(") ")
    cleaned_first_line = first_line_split[1]
    location = cleaned_first_line

    for h in location_info[1:]:
        cleaned_line = h.lstrip(" ")
        location += cleaned_line

    if "COMMENT" in location:
        updated_location = location.split("COMMENT")
        location = updated_location[0]

    if "URL" in location:
        updated_location = location.split("URL")
        location = updated_location[0]

    # parse author info
    author_line = [a for a in info if a.startswith('REFERENCE   2') or a.startswith('  TITLE     Direct Submission')]
    author_line_indices = []
    for q in author_line:
        index_number_ = info.index(q)
        author_line_indices.append(index_number_)

    index_reference = author_line_indices[0]
    index_author = index_reference + 1
    index_title = author_line_indices[-1]
    author_info = info[index_author:index_title]

    try:
        first_author_line = author_info[0]
        first_author_split = first_author_line.split("THORS   ")
        cleaned_first_author_line = first_author_split[1]
        author = cleaned_first_author_line

        for k in author_info[1:]:
            cleaned_line = k.lstrip(" ")
            author += cleaned_line


    except IndexError as error:
        author_line_indices = []
        author_line = [g for g in info if g.startswith('  AUTHORS') or g.startswith('  TITLE     Direct Submission')]
        for l in author_line:
            index_number_ = info.index(l)
            author_line_indices.append(index_number_)

        index_author = author_line_indices[0]
        index_title = author_line_indices[-1]
        author_info = info[index_author:index_title]

        first_author_line = author_info[0]
        first_author_split = first_author_line.split("THORS   ")
        cleaned_first_author_line = first_author_split[1]
        author = cleaned_first_author_line

        for k in author_info[1:]:
            cleaned_line = k.lstrip(" ")
            author += cleaned_line



    print(f"For accession number {i}, the journal was published by {author} from {location}")
    location_list.append(location)

    #bar.next()
#bar.finish()


address_line = []
regex_ = "([\sa-zA-Z0-9-]+,[\sa-zA-Z]+[\s0-9a-zA-Z-]+,[a-zA-Z\s]+)$"
for i in location_list:
    result = re.findall(regex_, i)
    final_result = " ".join(result)
    address_line.append(final_result.strip())

print(address_line)
