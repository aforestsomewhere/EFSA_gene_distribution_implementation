#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Download reference genomes and type strains from
NCBI's FTP server to assist in the creation of
the database for taxonomic classification of the 
MoPS pipeline.

USAGE EXAMPLES:

- Basic usage

To download a species from the from NCBI's Bacteria Refseq database, run the following:

    python download.py -g bacteria -n "Klebsiella pneumoniae" -o path/to/output_directory

You can also download all Refseq genomes for the genus _Bacillus_ by running:

    python download.py -g bacteria -n "Bacillus" -o path/to/output_directory

- Multiple species

You can also provide multiple species names or genera for download by running:

    python download.py -g bacteria -n "Klebsiella pneumoniae,Salmonella enterica,Bacillus" -o path/to/output_directory

The script also accepts a file containing 
the names of the species or genera, 
one organism per line, for example:

    python download.py -g bacteria -n my_species.txt -o path/to/output_directory

- Type material

In order to obtain type strains you have to run the following:

    python download.py -g bacteria -n "Klebsiella pneumoniae" -o path/to/output_directory --type-material type

The script offers the following options:

- type - deposited in at least two different culture collections;
- synonym - the sequences were derived from synonym type material;
- neotype - replacement culture for a type that has been lost;
- reference - the sequences were derived from reference material where type material 
              never was available and is not likely to ever be available.

Sources:

- Type material in the NCBI Taxonomy Database - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383940/
- Assembly summary report README - https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_assembly_summary.txt

- Parallel downloads

The scripts has uses by default 2 threads to perform the download of genomes.
If you have a good connection you can increase it by running (be careful to not get your IP banned!):

    python download.py -g bacteria -n "Klebsiella pneumoniae" -o path/to/output_directory --threads 4

- Retries

By default, the amount of retries for a failed download is 7,
however it may be increased/decreased by running:

    python download.py -g bacteria -n "Klebsiella pneumoniae" -o path/to/output_directory --retries 10

- File extensions and NCBI database

The script provides multiple options for file extensions
that are available in the FTP server, which are:

- genomic.fna.gz
- assembly_report.txt
- cds_from_genomic.fna.gz
- feature_count.txt.gz
- feature_table.txt.gz
- genomic.gbff.gz
- genomic.gff.gz
- genomic.gtf.gz
- protein.faa.gz
- protein.gpff.gz
- rna_from_genomic.fna.gz
- translated_cds.faa.gz

You can change the file extension by running:

    python download.py -g bacteria -n "Klebsiella pneumoniae" -o path/to/output_directory --file-extension protein.faa.gz

You are also able to download genomes from Genbank if 
your initial query does not have matches on Refseq, just run:

    python download.py -g bacteria -n "Klebsiella pneumoniae" -o path/to/output_directory --ncbi_db genbank

## Output

The downloaded genomes will be on the NCBI database directory
and also a tab-delimited metadata file for each downloaded genome.

"""

import os
import csv
import sys
import time
import socket
import hashlib
import argparse
import requests
import itertools
import urllib.request
from io import StringIO
from urllib.parse import urlparse
from appdirs import user_cache_dir
from tqdm import tqdm
from tqdm.contrib.concurrent import thread_map

# set socket timeout for urllib calls
socket.setdefaulttimeout(30)

# Get the user's cache dir in a system-independent manner
CACHE_DIR = user_cache_dir(appname="type_strain_assembly_report", appauthor="efsa")

# URL for NCBI's ftp server
FTP_URI = "https://ftp.ncbi.nih.gov/genomes"


def md5sum(filename):
    """Calculate the md5sum of a file and return the hexdigest.

    Parameters
    ----------
    filename : str
        Path to a file.

    Returns
    -------
    hash_md5_hexdigest: str
        The hexadecimal equivalent of the md5 hash.
    """
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as handle:
        for chunk in iter(lambda: handle.read(4096), b""):
            hash_md5.update(chunk)

    hash_md5_hexdigest = hash_md5.hexdigest()
    return hash_md5_hexdigest


def download_summary_file(FTP_URI, ncbi_db, group, use_cache=False):
    """Downloads the assembly report file from
    the input database.
    Alternatively, reads a previously downloaded
    assembly report.

    Parameters
    ----------
    FTP_URI : str
        URL for NCBI's ftp server.

    ncbi_db : str
        NCBI database to download from.

    group : str
        Taxonomic group.

    Returns
    -------
    assembly_report_content: StringIO object
        Contains the data from the downloaded
        assembly report in a StringIO object.
    """

    assembly_report_file = "{0}_{1}_assembly_summary.txt".format(ncbi_db, group)
    assembly_report_url = "{0}/{1}/{2}/assembly_summary.txt".format(
        FTP_URI, ncbi_db, group
    )

    # Save the assembly report file in cache
    # This will avoid downloading the file again
    # in subsequent runs
    cache_path = os.path.join(CACHE_DIR, assembly_report_file)
    if use_cache and os.path.exists(cache_path):
        print(f"Using cache: {cache_path}")
        with open(cache_path, "r") as cp:
            return StringIO(cp.read())

    req = requests.get(assembly_report_url)

    try:
        os.makedirs(CACHE_DIR)
    except OSError as err:
        if err.errno != 17:
            raise

    with open(cache_path, "w") as cpw:
        cpw.write(req.text)

    assembly_report_content = StringIO(req.text)

    return assembly_report_content


def clean_summary_file(assembly_report_content):
    """Cleans the assembly report.

    Parameters
    ----------
    assembly_report_content : str
        URL for NCBI's ftp server.

    Returns
    -------
    clean_content: list
        Contains the clean data
        from the downloaded assembly report.
    """

    # Skips first line because it only contains general information
    next(assembly_report_content)

    headers = next(assembly_report_content)[2:].split("\t")

    clean_content = []

    for row in assembly_report_content.readlines():
        clean_row = row.rstrip("\n").split("\t")
        row_data = {key: value for key, value in zip(headers, clean_row)}
        clean_content.append(row_data)

    return clean_content


def in_organism_name_list(species, organism_list):
    """Search for a species name in a list of organisms.

    Parameters
    ----------
    species : str
        Name of a species.

    organism_list : list
        List of organisms.

    Returns
    -------
    True if species is in list of organisms, False otherwise
    """

    for organism in organism_list:
        if species.startswith(organism):
            return True
        elif species.startswith(organism.capitalize()):
            return True
    return False


def filter_content(content, organism_name, type_material):
    """Filters content from the assembly_report file.

    Parameters
    ----------
    content : str
        An url to download a file.

    organism_name : str
        The name of the organissm to filter from the
        content file or a path to a file containing a
        list organism names to be downloaded.

    type_material : list
        List containing the relation to type material
        for the assembly.

    Returns
    -------
    new_entries : list
        A list containing the entries that match the
        provided parameters.
    """

    # Type material terms to filter the content
    type_material_filter_terms = {
        "type": "assembly from type material",
        "synonym": "assembly from synonym type material",
        "neotype": "assembly designated as neotype",
        "reference": "assembly designated as reftype",
    }

    # Check if organism_name is a file or not
    organism_name_list = []

    if os.path.isfile(organism_name):
        with open(organism_name, "r") as handle:
            input_list = handle.read().splitlines()
        organism_name_list.extend(input_list)
    elif isinstance(organism_name, str):
        organism_name_list.extend(organism_name.split(","))

    new_entries = []
    for entry in content:
        # Check if entries have the provided type material
        if type_material != "any":
            type_material_filter_term = [
                type_material_filter_terms[term] for term in [type_material]
            ]
            if (
                not entry["relation_to_type_material"]
                or entry["relation_to_type_material"] not in type_material_filter_term
            ):
                continue

        # Check if the entries have the organism name
        if organism_name_list and not in_organism_name_list(
            entry["organism_name"], organism_name_list
        ):
            continue

        # Check if the entries have a ftp url
        if entry["ftp_path"] == "na":
            continue

        # We don't want genomes that were
        # excluded from refseq
        if entry["excluded_from_refseq"] != "na":
            continue

        new_entries.append(entry)

    return new_entries


def get_md5_hash(content_hit_ftp_path):
    """Download and parse the md5checksums
    file for each hit.

    Parameters
    ----------
    content_hit_ftp_path : str
        URL for the md5checksums file
        in NCBI's FTP server.'

    Returns
    -------
    checksums_list: list
        A list of dictionaries containing
        the md5 hashes and the respective
        filenames.
    """

    # List to store download errors
    error_list = []

    # Build the md5 file full url
    md5_url = "{}/md5checksums.txt".format(content_hit_ftp_path)

    # Download the md5checksums file
    md5_req = requests.get(md5_url)

    # Get the content of the response from the request
    md5_req_text = md5_req.text

    # Parse the file and return a dictionary
    # with hashes and filenames
    checksums_list = []

    filename_error =  content_hit_ftp_path.split("/")[-1]

    for line in md5_req_text.split("\n"):
        # Skip empty lines
        if line == "":
            continue

        try:
            checksum, filename = line.split()
        except ValueError:
            if filename_error not in error_list:
                print("error downloading file: ", content_hit_ftp_path.split("/")[-1])
                error_list.append(filename_error)
            continue
        if filename.startswith("./"):
            filename = filename.replace("./", "")
        checksums_list.append({"checksum": checksum, "file": filename})

    return checksums_list


def has_file_changed(checksums_list, output_directory, file_extension):
    """Checks if the checksum of a given file has changed.

    Parameters
    ----------
    checksums_list : list
        List of checksums downloaded
        from the FTP server.

    output_directory : str
        Path to the output directory.

    file_extension : str
        Type of file to find.

    Returns
    -------
    True if file has changed, False otherwise.

    """

    for checksum in checksums_list:
        if not checksum["file"].endswith(file_extension):
            # not the file we want
            continue

        if file_extension == "genomic.fna.gz":
            if ("cds_from" in checksum["file"]) or ("rna_from" in checksum["file"]):
                # also not the file we want
                continue

        filename = checksum["file"]
        expected_checksum = checksum["checksum"]

    out_filename = os.path.join(output_directory, filename)

    # if file doesn't exist, it has changed
    if not os.path.isfile(out_filename):
        return True

    actual_checksum = md5sum(out_filename)

    return expected_checksum != actual_checksum


def download_file(url, file_name, retry):
    """Accepts a URL to download a file.

    Parameters
    ----------
    url : str
        An url to download a file.

    file_name : str
        The name of the file to be downloaded.

    retry : int
        Maximum number of retries if download fails.

    Returns
    -------
    response : str
        A string indicating that the download failed or
        an object with the response information for a
        successful download.
    """

    tries = 0
    while tries < retry:
        try:
            response = urllib.request.urlretrieve(url, file_name)
            break
        except Exception:
            response = "Failed: {0}".format(file_name)
            tries += 1
            print("Retrying {0} ...{1}".format(file_name.split("/")[-1], tries))
            time.sleep(1)

    return response


def write_metadata(metadata, downloaded_files, output_directory, file_extension):
    """Writes a tab-delimited file with
    metadata for each downloaded genome.

    Parameters
    ----------
    metadata: list
        Metadata for each downloaded genome.

    downloaded_files: list
        Names of the files that were
        actually downloaded.

    output_directory : str
        Path to the output directory.

    Returns
    -------
    None
        Writes a tab-delimited file in the
        output directory.
    """

    metadata_header = metadata[0].keys()
    metadata_filename = "{0}/metadata.tsv".format(output_directory)

    # Print the names of the files that were not downloaded
    for m in metadata:
        file_to_be_downloaded = "{0}_{1}".format(
            urlparse(m["ftp_path"]).path.rsplit("/")[-1], file_extension
        )
        if file_to_be_downloaded not in downloaded_files:
            print(
                "The following file was not downloaded: {0}, {1}, {2}".format(
                    m["organism_name"], m["assembly_accession"], m["asm_name"]
                )
            )

    if os.path.exists(metadata_filename):
        # Open the existing metadata file and append new data
        with open(metadata_filename, "a+") as metadata_append_obj:
            metadata_append_writer = csv.DictWriter(
                metadata_append_obj, fieldnames=metadata_header, dialect="excel-tab"
            )
            metadata_append_writer.writerows(metadata)
    else:
        # Create a new metadata file
        with open(metadata_filename, "w") as metadata_file:
            writer = csv.DictWriter(
                metadata_file, fieldnames=metadata_header, dialect="excel-tab"
            )
            writer.writeheader()
            writer.writerows(metadata)


def main(
    group,
    organism_name,
    output_directory,
    file_extension,
    ncbi_db,
    type_material,
    threads,
    retries,
):

    # Create output directory
    full_output_directory = os.path.join(output_directory, ncbi_db)
    if not os.path.isdir(full_output_directory):
        os.makedirs(full_output_directory)

    # Download assembly report file from the input NCBI database
    assembly_report_content = download_summary_file(FTP_URI, ncbi_db, group)

    # Parse the assembly report file
    clean_content = clean_summary_file(assembly_report_content)

    # Get data that matches the organism_name provided
    content_hits = filter_content(clean_content, organism_name, type_material)

    if len(content_hits) < 1:
        sys.exit(
            "There were no matches with the filter you provided. Please check your options."
        )

    content_hits_ftp_path = [hit_ftp_path["ftp_path"] for hit_ftp_path in content_hits]

    tqdm.write(
        "Found {0} hit(s). Downloading the MD5 checksum file(s).".format(
            len(content_hits)
        )
    )
    checksums_list = thread_map(
        get_md5_hash, content_hits_ftp_path, max_workers=threads
    )

    # Flatten the list for better manipulation
    checksums_list_flattened = list(itertools.chain(*checksums_list))

    # Create filenames and download urls
    ftp_urls = []
    assembly_ids = []
    content_hits_to_be_downloaded = []

    for hit in content_hits:

        # Check if entries have been already downloaded
        if has_file_changed(
            checksums_list_flattened, full_output_directory, file_extension
        ):

            # Retrieves the last part of the url path
            # This allows the correct URL to be used
            url_download = urlparse(hit["ftp_path"]).path.rsplit("/")[-1]

            ftp_url = "{0}/{1}_{2}".format(
                hit["ftp_path"], url_download, file_extension,
            )
            assembly_id_output_path = "{0}/{1}_{2}".format(
                full_output_directory, url_download, file_extension,
            )

            ftp_urls.append(ftp_url)
            assembly_ids.append(assembly_id_output_path)
            content_hits_to_be_downloaded.append(hit)

    if len(ftp_urls) < 1:
        sys.exit("The genomes that matched your filter have already been downloaded.")

    # Download files
    tqdm.write(
        "Starting download of {0} new genomes with {1} threads.".format(
            len(ftp_urls), threads
        )
    )

    res = thread_map(
        download_file,
        ftp_urls,
        assembly_ids,
        itertools.repeat(retries),
        max_workers=threads,
    )

    # Get the names of the files that were actually downloaded
    downloaded_files = [os.path.split(dl_file[0])[1] for dl_file in res]

    # Write metadata file
    write_metadata(
        content_hits_to_be_downloaded,
        downloaded_files,
        full_output_directory,
        file_extension,
    )


def parse_arguments():

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-g",
        "--group",
        type=str,
        required=True,
        choices=["bacteria", "fungi", "viral"],
        dest="group",
        help="Taxonomic group to search for.",
    )

    parser.add_argument(
        "-n",
        "--organism-name",
        type=str,
        required=True,
        dest="organism_name",
        help="Organism name. For example: 'Bacillus mycoides'.",
    )

    parser.add_argument(
        "-o",
        "--output-directory",
        type=str,
        required=True,
        dest="output_directory",
        help="Path to the directory to which downloaded files will be stored.",
    )

    parser.add_argument(
        "--fe",
        "--file-extension",
        type=str,
        required=False,
        choices=[
            "genomic.fna.gz",
            "assembly_report.txt",
            "cds_from_genomic.fna.gz",
            "feature_count.txt.gz",
            "feature_table.txt.gz",
            "genomic.gbff.gz",
            "genomic.gff.gz",
            "genomic.gtf.gz",
            "protein.faa.gz",
            "protein.gpff.gz",
            "rna_from_genomic.fna.gz",
            "translated_cds.faa.gz",
        ],
        default="genomic.fna.gz",
        dest="file_extension",
        help="Choose file type to download through extension.",
    )

    parser.add_argument(
        "--ncbi_db",
        type=str,
        required=False,
        choices=["refseq", "genbank"],
        default="refseq",
        dest="ncbi_db",
        help="The script can search for the files to "
        "download in RefSeq or Genbank or both "
        "(will only search in Genbank if download "
        "from RefSeq fails).",
    )

    parser.add_argument(
        "--type-material",
        type=str,
        required=False,
        choices=["any", "type", "synonym", "neotype", "reference"],
        default="any",
        dest="type_material",
        help="The relation to type material for the assembly.",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=2,
        dest="threads",
        help="Number of threads for download.",
    )

    parser.add_argument(
        "-r",
        "--retries",
        type=int,
        required=False,
        dest="retries",
        default=7,
        help="Maximum number of retries when a download fails.",
    )

    parser.add_argument(
        "-m",
        "--metaonly",
        type=int,
        required=False,
        dest="retries",
        default=7,
        help="Create only the metadata file without downloading the genomes",
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
