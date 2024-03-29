import os
import pandas as pd
from snakemake.utils import validate
import datetime
import time

report: "../report/workflow.rst"

configfile: 'config/config.yaml'
validate(config, schema="../schemas/config.schema.yaml")

# Setup vars
config_threads = int(config["threads"])
foldername = os.path.basename(os.getcwd())

# Load in metadata
metadata_file = "config/metadata.tsv"
metadata_df = pd.read_csv(metadata_file, sep = "\t").set_index('sample_name', drop=False)
# validate(metadata_df, schema="schemas/metadata_schema.yaml")

# Setup samplesheet
samples = metadata_df.index.tolist()

# Get contrasts
contrasts = config["deseq2"]["contrasts"].split(",")

# Get the raw fastq files to start with
def get_fastq(wildcards):
    return metadata_df.loc[(wildcards.sample), ["r1", "r2"]]

def get_deg_input(wildcards):
	return ["results/counts/{}/quant.sf".format(sample) for sample in samples]

def get_gmt_file(wildcards):
	return 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/{}.all.v{}.symbols.gmt'.format(wildcards.collection, config['gsea']['msigdb_version'])
