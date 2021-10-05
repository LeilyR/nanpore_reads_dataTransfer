# nanoporeReads_dataTransfer
A pipeline to transfer the Nanopore reads to the end users

# Usage
```bash
module load anaconda3;
python pipeline.py -i flowcell -r organism
```

# Example input path:
```bash
../path/to/flowcell/
├── drift_correction_FAL82553_acfbfd23.csv
├── duty_time_FAL82553_acfbfd23.csv
├── fast5
├── final_summary_FAL82553_acfbfd23.txt
├── mux_scan_data_FAL82553_acfbfd23.csv
├── report_FAL82553_20200319_1033_c25a4aea.md
├── Samplesheet.csv
├── sequencing_summary_FAL82553_acfbfd23.txt
└── throughput_FAL82553_acfbfd23.csv
```
# Example output path for bioinfo. core:
```bash
../output_path/to/flowcell/
├── drift_correction_FAL82553_acfbfd23.csv
├── duty_time_FAL82553_acfbfd23.csv
├── fast5
├── final_summary_FAL82553_acfbfd23.txt
├── mux_scan_data_FAL82553_acfbfd23.csv
├── report_FAL82553_20200319_1033_c25a4aea.md
├── Samplesheet.csv
├── sequencing_summary_FAL82553_acfbfd23.txt
├── throughput_FAL82553_acfbfd23.csv
├── FASTQC_Project_number_User_PI
│   └── Sample_its_name
|       ├── NanoComp_20200408_1503.log
│       ├── NanoComp_lengths.png
│       ├── NanoComp_log_length.png
│       ├── NanoComp_number_of_reads.png
│       ├── NanoComp_OverlayHistogram.html
│       ├── NanoComp_OverlayHistogram_Normalized.html
│       ├── NanoComp_OverlayLogHistogram.html
│       ├── NanoComp_OverlayLogHistogram_Normalized.html
│       ├── NanoComp_quals.png
│       ├── NanoComp-report.html
│       ├── NanoComp_total_throughput.png
│       └── NanoStats.txt
└── Project_number_User_PI
    └── Sample_its_name
        └── Sample_name.fastq.gz
```
# Example output path for a user:
```bash
../user_path/to/flowcell/
├── FASTQC_Project_number_User_PI
│   └── Sample_its_name
|       ├── NanoComp_20200408_1503.log
│       ├── NanoComp_lengths.png
│       ├── NanoComp_log_length.png
│       ├── NanoComp_number_of_reads.png
│       ├── NanoComp_OverlayHistogram.html
│       ├── NanoComp_OverlayHistogram_Normalized.html
│       ├── NanoComp_OverlayLogHistogram.html
│       ├── NanoComp_OverlayLogHistogram_Normalized.html
│       ├── NanoComp_quals.png
│       ├── NanoComp-report.html
│       ├── NanoComp_total_throughput.png
│       └── NanoStats.txt
├── Project_number_User_PI
|   └── Sample_its_name
|       └── Sample_name.fastq.gz
└── Analysis_number_User_PI
    └── mapping_on_organism
        ├── organism_genome.fa
        ├── organism_genome.mmi
        ├── sample.bam
        └── sample.bam.bai
```
