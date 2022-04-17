from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from progress.bar import Bar
import pandas as pd
import geopandas
import matplotlib.pyplot as plt
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
    # if the xml file which contains the blast outputs does not exist, use qblast to invoke the NCBI BLAST server over the internet and save the outputs on the local machine as an xml file:
    if not (os.path.exists("results.xml")):
        sequence_data = open(sequence_file).read()
        bar = Bar("Processing sequence data... ")
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data) #blastn = program, "nt" = database to search against, in this case, "nt" stands for the nucleotide databse
        bar.finish()

        with open("results.xml","w") as save_file: # write the blastn output to a xml file
            blast_results = result_handle.read()
            save_file.write(blast_results)

    results = open("results.xml","r")
    records = list(NCBIXML.parse(results)) # parse results (xml file) into a list
    one_query = records[0] #from a list to a record class that holds all BLAST output

    accession_numbers = [] #create empty list to store the accession_numbers of each organism

    for i in range(len(one_query.alignments)): #loop through every query to obtain the accession_number
        one_hit = one_query.alignments[i]
        accession_numbers.append(one_hit.accession)

    return accession_numbers

def getLocationsAuthors(accession_numbers: list) -> list:
    '''
    Function to use Entrez to parse location and first author data for the each
    of 50 top hits from BlastN output and return that data in a list

    Parameters:
        accession_numbers: list -- accession numbers for 50 top hits

    Results:
        location of institution and first author associated with 50 top hits
        from BlastN output
    '''

    Entrez.email = input("Enter your email address (so that NCBI can contact you if there's a problem): ")
    bar = Bar("Processing accession numbers in GenBank:", max = len(accession_numbers))

    # create two empty lists that store the location and author informations
    location_list = []
    author_list = []

    # parse GenBank data for every matching sequence (accession number)
    for i in accession_numbers:
        handle = Entrez.efetch(db = 'nucleotide', rettype = 'gb', retmode = 'text', id = i)
        record = handle.read()
        info = record.split("\n") # split the text file by new line
        location_line = [i for i in info if i.startswith('  JOURNAL   Submitted') or i.startswith('FEATURES') or i.startswith('COMMENT')] # extract lines that starts with '  JOURNAL   Submitted' which contains the location information, and lines that starts with 'FEATURES' or 'COMMENT' which are the lines right after the location information
        location_line_indices = []

        for j in location_line:
            index_number = info.index(j) #extract the index of the locations line
            location_line_indices.append(index_number)

        index_1 = location_line_indices[0] # first location line
        index_2 = location_line_indices[-1] # last location line

        location_info = info[index_1:index_2]
        first_line = location_info[0]
        first_line_split = first_line.split(") ") # all location lines begin with 'JOURNAL Submitted (date)', so split the line accordingly and only store the relevant information
        cleaned_first_line = first_line_split[1]
        location = cleaned_first_line

        for h in location_info[1:]: # for the rest of the location lines
            cleaned_line = h.lstrip(" ") # delete the extra spaces in front of the first letter
            location += cleaned_line # append the info to the first line

        # filter out all possible irrelevant information:
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
           index_author = index_reference + 1 # author information starts on the next line
           index_title = author_line_indices[-1] # title line is the line right after the author info
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


        except IndexError as error: # the line reference 2 may not exist if there is only one reference
            author_line_indices = []
            author_line = [g for g in info if g.startswith('  AUTHORS') or g.startswith('  TITLE     Direct Submission')] # search for line that starts with 'AUTHORS' instead of 'Reference 2'
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

def getAddress(locations: list) -> list:
    '''
    Function to use regex to pull the valid address from each of the more specific
    locations parsed by GenBank

    Args:
        locations: list of all locations parsed from GenBank

    Returns:
        list of valid addresses

    * Note: this regex will not work for all locations returned from GenBank,
            but invalid addresses will be caught by getLatLongLists
    '''
    address_list = []
    regex_ = "([\sa-zA-Z0-9-]+,[\sa-zA-Z]+[\s0-9a-zA-Z-]+,[a-zA-Z\s]+)$" # regular expression to get the specific address that exludes the institution names and in some cases the street name as well
    for i in locations:
        result = re.findall(regex_, i)
        final_result = " ".join(result)
        cleaned_result = final_result.strip()
        split_address = cleaned_result.split(",")
        first, second = split_address[0], split_address[2] # only need the street/city and the country. The zip code info which belongs to split_address[1] is eliminated as Nominatim is unable to process zip code.
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
    Function to get the latitude and longitude for each address and return them
    as a list containing the latitude list and the longitude list

    Args:
        address_list: list -- list of all addresses to be processed

    Returns:
        list containing a list of latitudes for every address and a list of
        longitudes for every address
    '''
    lat_list = []
    long_list = []
    bar = Bar("Processing each address in Nominatim:", max = len(address_list))
    for i in address_list:
        try:
            lat_long = getLatLong(i)
        except AttributeError as error:
            print(" ")
            print(f"Unable to fetch the following address using Nominatim: {i}")
            new_address = input("Enter valid address in 'city, country' format: ") # in case a foreign street name was the problem that Nominatim cannot process the address info
            try:
                lat_long = getLatLong(new_address)
            except AttributeError as error:
                print(" ")
                print("Invalid Address entered")
                lat_long = ["Invalid Address", "Invalid Address"]
        lat_list.append(lat_long[0])
        long_list.append(lat_long[1])
        bar.next()
    bar.finish()
    return [lat_list, long_list]

def makeDataDict(lat_long_list: list, location_author_list: list):
    blast_dict = {'authors': location_author_list[1], 'address': location_author_list[0], 'Latitude': lat_long_list[0], 'Longitude': lat_long_list[1]}
    return blast_dict
    '''
#Making graph interactive
def makeInteractive(blast_data_dict: dict):
#Converting dictionary to dataframe and then to GeoDataFrame
    df = pd.DataFrame(blast_data_dict)
    gdf = geopandas.GeoDataFrame(
        df, geometry = geopandas.points_from_xy(df.Longitude, df.Latitude))
#mapping using .explore()
    # world = geopandas.read_file(geopandas.geopandas.datasets.get_path(gdf))
    # ax = world.plot(colot = 'white', edgecolor = black)
    author = blast_data_dict[0]
    address = blast_data_dict[1]
    #other tiles
    gdf.explore("geometry", cmap = 'Set2', Legend = False,
                tooltip = False, popup = ['author','address'])

    # gdf.plot(ax=ax, color = 'red')
    # plt.show()
    '''

def main():
    accession_numbers = runBlast("blast.fasta")
    location_author_list = getLocationsAuthors(accession_numbers)
    address_list = getAddress(location_author_list[0])
    # print(address_list)
    latlong_list = getLatLongLists(address_list)
    # print(latlong_list)
    '''
    blast_data_dict = makeDataDict(latlong_list, location_author_list)
    makeInteractive(blast_data_dict)
    '''

if __name__ == "__main__":
    main()
