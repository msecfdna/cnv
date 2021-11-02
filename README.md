# CNV Pipeline

The CNV method implemented here combines the information from on-target and off-target fragments to assign a final copy number index or Z-score to a set of pre-defined genomic regions.


# Input: 

The script takes:

(1) Directory of barcode-deduped bams (on-target information) for target samples<br /> 
(2) Directory of samtools-deduped bams (off-target information) for target samples <br /> 
(3) Path to the selector (only the targeted part- the padded part is created automatically)<br /> 
(4) Path to human genome reference FASTA file <br /> 
(5) Path to human genome reference bed file (2-column file with chr order and size)<br /> 
(6) Directory of barcode-deduped bams for control cohort 1<br /> 
(7) Directory of samtools-deduped bams for control cohort 1<br /> 
(8) Directory of barcode-deduped bams for control cohort 2<br /> 
(9) Directory of samtools-deduped bams for control cohort 2<br />
(10) Path to the annotation file (format: 4-column bed file with header as: "#chr"    "start"   "end"     "name")<br />
(11) Number of processes in multiprocessing using Pool (this is to speed up the coverage calculation part which could be very slow- default: 16)
(11) Directory in which the final Z-score files and plots would be stored (optional- if not provided everything will be written in the **target sample barcode-deduped directory**, i.e. number 1 above)<br />
(12) Off-target bin size (optional- default: 100kb)<br /> 
(13) Path to targeted background for control cohort 1 (required **if 6 is not provided**)<br /> 
(14) Path to off-target background for control cohort 1 (required **if 7 is not provided**)<br /> 
(15) Path to targeted background for control cohort 2 (required **if 8 is not provided**)<br /> 
(16) Path to off-target background for control cohort 2 (required **if 9 is not provided**)<br /> 
(17) Maximum fragment length (optional- default: 1000bp)<br /> 
(18) Minimum mapping quality (optional- default: 20)<br /> 

# Output:

The algorithm generates various intermediate files: raw coverage --> GC-bias corrected coverage --> log-transformed normalized coverage --> sample-specific STD normalized --> annotated by the annotation file  

# Example:

The current version does argument parsing and can be run as follows:

python3 cnv-pipeline.py -b \\<path to barcode-dedupe bams\\> -B <path to samtools deduped bams> -s <path to selector file> -g <human genome reference FASTA> -G <human genome reference bed> -m <path to barcode deduped bams for control cohort I> -M <path to samtools deduped bams for control cohort I> -n <path to barcode deduped bams for control cohort II> -N <path to samtools deduped bams for control cohort II> -i <sample information file to match barcode and samtools deduped files> -a <path to annotation file> -x <path to the reference centromeres to filter> -r <number of processes to multiprocess the coverage calculations> -o <output directory to write the final Z-scores and .png files> 
