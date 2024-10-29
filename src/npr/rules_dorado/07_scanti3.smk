'''
add QC reports for RNAseq post-alignment 
DeepSeq suggest to use SQANTI3 to get an idea about
Avg. reads/exon
Avg. completeness/gene
No. of transcripts found
many other stats
Could be good to have this as part of the reports
https://github.com/ConesaLab/SQANTI3

To run scanti3 we also need to run flair to generate a gtf file
https://github.com/BrooksLabUCSC/flair
'''

# define source and target pattern
source_fq = sample_dat + ".fastq.gz"
logpat = sample_log + "_scanti3.log"
target = sample_ana + ".isoforms.gtf"
bchpat = sample_bch + "_flair.tsv"

rule modbed_final:
    input: expand_project_path(source_fq)
    output: touch("flags/07_flair.done")
    
rule bam2modbed:
    input:
        flag = "flags/07_flair.done",
        fq = source_fq
    output:
        bed = target
    log:
        logpat
    params:
        rna = config["info_dict"]["barcode_kit"]['protocol']
        genome =  config['genome'].get(org, None)
    benchmark:
        bchpat
    threads: 10
    conda:
        "ont-ppp-flair"
    shell:'''
        echo "run flair" 2>> {log}
        if [[ "{params.rna}" == "rna" ]]; then
            echo flair 1234 -t {threads} -r {input.fq} -g {params.genome} -o {output.bed} 2>> {log}
            flair 1234 -t {threads} -r {input.fq} -g {params.genome} -o {output.bed} 2>> {log}
        else
            echo "flair skipped" 2>> {log}
        fi
    '''