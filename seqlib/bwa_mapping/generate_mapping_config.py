import os
import argparse
import json

"""
Generates a config file of the following format:
{
    "bwa_index": "/global/home/users/sweram19/scratch/references/UCSC/pan_troglodytes/panTro6/index/bwa/panTro6.fa",
    "input_sample_path": "/global/home/users/sweram19/scratch/chimp/SRP018689/",
    "input_samples": {
        "SRS396861_Damian_Pan-troglodytes-ellioti_M": [
            "SRR748133_1.fastq",
            "SRR748133_2.fastq",
            "SRR748134_1.fastq",
            "SRR748134_2.fastq",
        ],
        "SRS396866_Taweh_Pan-troglodytes-ellioti_M": [
            "SRR748156_1.fastq",
            "SRR748156_2.fastq",
            "SRR748157_1.fastq",
            "SRR748157_2.fastq",
        ]
    }
}

To be used with the following Snakefile:
/global/home/users/psudmant/code/snakemake-workflows/bio/ngs/workflows/bwa_map/run_maps_Snakefile

Note: Make sure input samples directory is only filled with directories, one for each sample. If any other files
or directories are in the input samples directory, this file will error.
"""

__author__ = "Swetha Ramesh"


def make_json(bwa, bwa_path, input_samples_path, config_name):
    j_out = {bwa: bwa_path}
    j_out['input_sample_path'] = input_samples_path
    samples = {}
    cwd = os.getcwd()
    dirs = os.listdir(input_samples_path)
    for d in dirs:
        files = []
        for f in os.listdir(os.path.join(input_samples_path, d)):
            files.append(f)
        samples[d] = files
    j_out['input_samples'] = samples
    FOUT = open("{cwd}/{config_name}.json".format(cwd=cwd, config_name=config_name),'w')
    FOUT.write(json.dumps(j_out, indent=4, separators=(",", ": ")))
    FOUT.close()


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config_name", default="config")
    parser.add_argument("--index", required=True, help="Full path to bwa index file")
    parser.add_argument("--samples", required=True, help="Full path to sample directories")
    o = parser.parse_args()


    make_json("bwa_index", o.index, o.samples, o.config_name)
