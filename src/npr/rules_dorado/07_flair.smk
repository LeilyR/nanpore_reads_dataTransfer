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
isoforms = sample_ana + ".isoforms.bed"
aligned = sample_ana + ".flair.aligned.bed"
corected = sample_ana + "_all_corrected.bed"
logpat = sample_log + ".flair.log"
bchpat = sample_bch + ".flair.tsv"

rule modbed_final:
    input: expand_project_path(source_fq)
    output: touch("flags/07_flair.done")
    
rule bam2modbed:
    input:
        flag = "flags/07_flair.done",
        fq = source_fq
    output:
        aligned = aligned
        corrected = corrected
        isof = isoforms
    log:
        logpat
    params:
        rna = config["info_dict"]["barcode_kit"]['protocol']
        genome =  config['genome'].get(org, None)
    benchmark:
        bchpat
    threads: 16
    conda:
        "ont-ppp-flair"
    shell:'''
        echo "run flair" 2>> {log}
        if [[ "{params.rna}" == "RNA" || "{params.rna}" == "cDNA" ]]; then
            echo flair 123 -t {threads} -r {input.fq} -g {params.genome} -o {output.bed} 2>> {log}
            flair align -t {threads} -r {input.fq} -g {params.genome} -o {output.aligned} 2>> {log}
            #echo flair correct -q {output.aligned} -g  {params.genome}  2>> {log}
            #flair correct -t {threads} -q {output.aligned} -g  {params.genome}  2>> {log}
            echo flair collapse -t {threads} -g {params.genome} -q {output.bed} -r {input.fq} -o {output.isoforms}  2>> {log}
            #flair collapse -t {threads} -g {params.genome} -q {output.corrected} -r {input.fq} -o {output.isoforms}  2>> {log}
            flair collapse -t {threads} -g {params.genome} -q {output.aligned} -r {input.fq} -o {output.isoforms}  2>> {log}
        else
            echo "no flair" 2>> {log}
        fi
    '''