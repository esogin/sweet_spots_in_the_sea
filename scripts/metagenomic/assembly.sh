#Assembly Script
#EM Sogin
#Seagrass rhizosphere microbiome
#Updated March 2020

## All libraries were concatinated prior to cleaning
lib=ALL_LIB;
ref=/path/to/illumina/adaptors/adapters.fa;

## Remove low quality reads and adaptors
bbduk.sh ref=${ref} ktrim=l mink=11 hdist=1 in1=all_reads_R1.fastq.gz in2=all_reads_R2.fastq.gz out=all_reads_ktriml.fq.gz;
bbduk.sh ref=${ref} ktrim=r trimq=2 qtrim=rl minlength=50 mink=11 hdist=1 in=all_reads_ktriml.fq.gz out=all_reads_q2_ktrimmed.fq.gz;

## Normalize read set to have a target depth of 100X and remove sequences less then 2X, considered them errors
bbnorm.sh -Xmx400g  in=all_reads_q2_ktrimmed.fq.gz out=highfreq_kmers.fq.gz target=100 min=2;

#Run metaspades for error correction and assembly
spades.py -o spades --12 highfreq_kmers.fq.gz -t 64 -m 1.8T --phred-offset 33 --only-error-correction;

#Assemble using megahit
megahit -1 spades/corrected/highfreq_kmers_1.00.0_0.cor.fastq.gz \
		-2 spades/corrected/highfreq_kmers_2.00.0_0.cor.fastq.gz \
		-r spades/corrected/highfreq_kmers__unpaired.00.0_0.cor.fastq.gz \
		-t 48 -o megahit --out-pre 3847 --k-min 21 --k-max 151 --k-step 10
#END