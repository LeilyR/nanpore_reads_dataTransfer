#snakefile
from npr.snakehelper import retRule
import pandas as pd

# create a pandas dataframe of samples to get the sample : project relationship
metadata = dict(config["data"])
del metadata["projects"]
del metadata["samples"]
metadata = pd.DataFrame(metadata).T

# wildcard mapping from wildcard name -> metadata column
wc_mapping = {
    "project" : "Sample_Project",
    "sample_id" : "Sample_ID",
    "sample_name" : "Sample_Name"
}       


def expand_project_path(path, metadata=metadata, wc_mapping=wc_mapping):
    """
    Given a path containing {project}, {sample_id}, and {sample_name}, 
    find the corresponding columns in the metadata object
    """
    mapping = {k : metadata[v] for k, v in wc_mapping.items()}
    return expand(path, zip, **mapping)

include: retRule("0_prepare.smk", config)
include: retRule("1_basecalling.smk", config)
include: retRule("2_rename.smk", config)
include: retRule("3_qc.smk", config)
include: retRule("4_data_transfer.smk", config)
include: retRule("5_mapping.smk", config)
include: retRule("6_multiqc.smk", config)

rule all:
    input:
        'flags/0_prepare.done',
        'flags/1_basecalling.done',
        'flags/2_renamed.done',
        'flags/3_qc.done',
        'flags/4_transfer.done',
        'flags/5_mapping.done'
        'flags/6_multiqc.done'
