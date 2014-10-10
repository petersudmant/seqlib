from gtf_to_genes import *
import logging
import argparse
import pdb

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_gtf")
    parser.add_argument("--fn_out_index")
    args = parser.parse_args() 

    
    pdb.set_trace()
    logger = logging.getLogger("test")
    index_file       = "/net/cpp-mirror/databases/ftp.ensembl.org/gtf.index"
    search_path_root = "/net/cpp-mirror/databases/ftp.ensembl.org"
    search_path_root = "/net/cpp-mirror/databases/ftp.ensembl.org/pub/release-56/gtf/homo_sapiens/"
    regex_input          = r"(.+\/)(([^.]+)\..+\.(.+)\.gtf(?:\.gz)?)$"

    # put cache file in same directory as GTF file
    cache_file_pattern   = r"\1\2.cache"

    #
    # uncomment this line to put cache file in same directory index file
    #
    #cache_file_pattern   = r"{INDEX_FILE_PATH}/\2.cache"

    #
    # Unique identifier per GTF file
    # e.g. "Anolis_carolinensis:56"
    #
    identifier_pattern   = r"\3:\4"

    index_gtf_files(index_file,
                    search_path_root,
                    regex_input,
                    cache_file_pattern,
                    identifier_pattern,
                    True,
                logger)
