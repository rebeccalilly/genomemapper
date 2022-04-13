from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from progress.bar import Bar
import os
import re
import geopy.geocoders
from geopy.geocoders import Nominatim
from Bio import Entrez

def runBlast(sequence_file: str) -> list:
    '''
    Function to run BlastN using the given sequence file as input and return a
    list of the accession numbers for 50 top hits

    Parameters:
        sequence_file: str -- name of DNA sequence file in FASTA format

    Returns:
        accession numbers for 50 top hits, as list
    '''
    if not (os.path.exists("results.xml")):
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

    return accession_numbers

def getLocations(accession_numbers: list) -> list:
    '''
    Function to use Entrez to parse location and first author data for the each
    of 50 top hits from BlastN output and return that data in a list

    Parameters:
        accession_numbers: list -- accession numbers for 50 top hits

    Results:
        location of institution and first author associated with 50 top hits
        from BlastN output
    '''
    #parsing GenBank data
    Entrez.email = input("Enter your email address (so that NCBI can contact you if there's a problem):")
    bar = Bar("Parsing location and author information for each accession numbers from GenBank... ", max = len(accession_numbers))

    location_list = []
    author_list = []
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

        location_list.append(location)

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

            author_list.append(author)


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

            author_list.append(author)
        bar.next()
    location_author_list = [location_list, author_list]
    bar.finish()
    return location_author_list

def getAddress(locations: list):
    '''
    '''
    address_list = []
    regex_ = "([\sa-zA-Z0-9-]+,[\sa-zA-Z]+[\s0-9a-zA-Z-]+,[a-zA-Z\s]+)$"
    for i in locations:
        result = re.findall(regex_, i)
        final_result = " ".join(result)
        cleaned_result = final_result.strip()
        split_address = cleaned_result.split(",")
        first, second = split_address[0], split_address[2]
        address = first, second
        valid_address = ', '.join(address)
        address_list.append(valid_address)
    return address_list

def getLatLong(address: str) -> list:
    '''
    Function to use Nominatim from the geopy library to return the latitude and
    longitude of a given address, as list

    Args:
        address: str -- specific address to find latitude and longitude for

    Returns:
        tuple containing latitude and longitude of address
    '''
    geolocator = Nominatim(user_agent="dcs211_randers2/1")
    location = geolocator.geocode(address)
    lat = location.latitude
    long = location.longitude
    return [lat,long]

def getLatLongLists(address_list: list) -> list:
    '''
    '''
    lat_long_list = []
    bar = Bar("Retrieving latitude and longitude for each address using Nominatim...", max = len(address_list))
    for i in address_list:
        try:
            lat_long = getLatLong(i)
        except AttributeError as error:
            print(f"Unable to fetch the following address using Nominatim: {i}")
            lat_long = "Invalid Address"
        lat_long_list.append(lat_long)
        bar.next()
    bar.finish()
    return lat_long_list

def main():
    accession_numbers = runBlast("blast.fasta")
    location_author_list = getLocations(accession_numbers)
    address_list = getAddress(location_author_list[0])
    print(address_list)
    latlong_list = getLatLongLists(address_list)
    print(latlong_list)

if __name__ == "__main__":
    main()
