# CNV Pipeline

The CNV method implemented here combines the information from on-target and off-target fragments to assign a final copy number index or Z-score to a set of pre-defined genomic regions.


# Input: 

The script takes:

(1) Directory of barcode-deduped bams (on-target information) for target samples
(2) Directory of samtools-deduped bams (off-target information) for target samples 
(3) Path to the selector (only the targeted part- the padded part is created automatically)
(4) Path to human genome reference FASTA file 
(5) Path to human genome reference bed file (2-column file with chr order and size)
(6) Directory of barcode-deduped bams for control cohort 1
(7) Directory of samtools-deduped bams for control cohort 1
(8) Directory of barcode-deduped bams for control cohort 2
(9) Directory of samtools-deduped bams for control cohort 2
(10) Path to the annotation file (format: 4-column bed file with header as: "#chr"    "start"   "end"     "name") 
(11) Directory in which the final Z-score files and plots would be stored (optional- if not provided everything will be written in the **target sample barcode-deduped directory**, i.e. number 1 above)
(12) Off-target bin size (optional- default: 100kb)
(13) Path to targeted background for control cohort 1 (required **if 6 is not provided**)
(14) Path to off-target background for control cohort 1 (required **if 7 is not provided**)
(15) Path to targeted background for control cohort 2 (required **if 8 is not provided**)
(16) Path to off-target background for control cohort 2 (required **if 9 is not provided**)
(17) Maximum fragment length (optional- default: 1000bp)
(18) Minimum mapping quality (optional- default: 20)

# Output:

The algorithm generates various intermediate files: raw coverage --> GC-bias corrected coverage --> log-transformed normalized coverage --> sample-specific STD normalized --> annotated by the annotation file  
