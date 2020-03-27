#Coassembly binning
#EM Sogin
#Seagrass rhizosphere microbiome
#Updated March 2020

#Bin coassembly using concoct, metabat, maxbin2, refine MAGs in dastool

##Map individual library reads back to coassembly
mkdir bam/
for lib in $libraries; 
	do
	#clean libraries
		bbduk.sh ref=${ref} ktrim=l mink=11 hdist=1 in1=${lib}_R1.fastq.gz in2=${lib}_R2.fastq.gz out=${lib}_ktriml.fq.gz;
		bbduk.sh ref=${ref} ktrim=r trimq=2 qtrim=rl minlength=50 mink=11 hdist=1 in=${lib}_reads_ktriml.fq.gz out=${lib}_reads_q2_ktrimmed.fq.gz;
		spades.py -o ${lib}_spades --12 ${lib}_reads_q2_ktrimmed.fq.gz -t 24 --phred-offset 33 --only-error-correction;
	#map to coassembly
		in1=${lib}_spades/corrected/${lib}_1.00.0_0.cor.fastq.gz;
		in2=${lib}_spades/corrected/${lib}_2.00.0_0.cor.fastq.gz
		bbmap.sh in=$in1 in2=$in2 ref=contigs.fa covstats="$lib"_COVSTATS \
				 out=bam/"$lib".bam scafstats="$lib"_scafstats.txt \
				 statsfile="$lib"_stderr;
		samtools sort -o bam/${lib}_sorted.bam bam/${lib}.bam
		samtools index bam/${lib}_sorted.bam;
	done

#Metabat binning
jgi_summarize_bam_contig_depths --minContigLength 1000 --outputDepth depth.txt \
								--pairedContigs paired.txt bam/*sorted.bam 
metabat2 -i contigs.fa -a depth.txt -o result/metabat -v 

cd result/
path/to/DAS_Tool-master/src/Fasta_to_Scaffolds2Bin.sh -e fa > ../metabat.scaffolds2bins.tsv

#Concoct binning
conda activate concoct_env
	cut_up_fasta.py contigs.fa  -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
	concoct_coverage_table.py contigs_10K.bed bam/*sorted.bam > coverage_table.tsv
	concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -t 48 -b result/
	merge_cutup_clustering.py result/clustering_gt1000.csv > result/clustering_merged.csv
	mkdir result/fasta_bins
	extract_fasta_bins.py contigs.fa result/clustering_merged.csv --output_path result/fasta_bins
conda deactivate concoct_env

cd result/
sed "s/,/\t/g" clustering_merged.csv > concoct.s2b.tsv
tail -n +2 concoct.s2b.tsv > file.tmp && mv file.tmp concoct.s2b.tsv

#maxbin2 binning
run_MaxBin.pl -thread 48  -contig contigs.fa -out maxbin \
			-abund 3847_A_COVSTATS_forMaxBin.txt -abund2 3847_B_COVSTATS_forMaxBin.txt \
			-abund3 3847_C_COVSTATS_forMaxBin.txt -abund4 3847_D_COVSTATS_forMaxBin.txt \
			-abund5 3847_E_COVSTATS_forMaxBin.txt -abund6 3847_F_COVSTATS_forMaxBin.txt \
			-abund7 3847_G_COVSTATS_forMaxBin.txt -abund8 3847_H_COVSTATS_forMaxBin.txt \
			-abund9 3847_I_COVSTATS_forMaxBin.txt
rsync -a /scratch/sogin/tmp.$JOB_ID /opt/extern/bremen/symbiosis/sogin/Data/SedimentMG/processed_reads/libraries/library_3847/coassembly/binning_v2/maxbin/
cd result/
Path/to/DAS_Tool-master/src/Fasta_to_Scaffolds2Bin.sh -e fasta > ../maxbin2.scaffolds2bin.tsv

#Prep data for DAS_Tool
## select only headers that were binned in assembly 
cut -f 1 concoct.s2b.tsv > fah.concont.txt
cut -f 1 maxbin2.scaffolds2bin.tsv > fah.maxbin.txt
cut -f 1 metabat.s2b.tsv > fah.metabat.txt
cat fah.concont.txt fah.maxbin.txt fah.metabat.txt |sort | uniq > fasta_headers.txt
seqtk subseq contigs.fa fasta_headers.txt > contigs.use.fa 

## Run DAS_Tool
DAS_Tool -i data/concoct.s2b.tsv,data/maxbin2.scaffolds2bin.tsv,data/metabat.s2b.tsv \
		 -l concoct,maxbin,metabat -c data/contigs.use.fa  -o result/DASToolRun \
		 --threads 48 --search_engine diamond --write_bins 1 

#Check quality of bins
#Run CheckM
checkm tree bins/ checkm_tree/ -t 24 -x .fa
checkm tree_qa ./checkm_tree -o 2 -f tree_qa 
checkm lineage_set ./checkm_tree markers
checkm analyze markers ./bins checkm_analyze -t 24 -x .fa
checkm qa markers ./checkm_analyze -o 2 -t 24 -f checkm_qa_results.txt

#Check for tRNAs with aragorn
# all in one line:
for n in *.fa; do echo $n; aragorn -t -fon $n | grep '>' | cut -d' ' -f2 \
	|sed 's/([a-z]*)//g'| sed 's/tRNA-//g' | sed 's/?(//g' | sed 's/)//g' \
	| sed 's/|/\n/g' | sort -u | wc -l; done > trnas_total

#Get 16S from select MAGs
#Quality trim all reads to Q27
bbduk.sh in1=all_reads_R1.fq.gz in2=all_reads_R2.fq.gz out=clean_r.fq.gz qtrim=r trimq=27
bbduk.sh in=clean_r.fq.gz out=clean_rl.fq.gz qtrim=l trimq=27

for bin in $select_mags; 
	do
		bbmap.sh in=clean_rl.fq.gz ref=$bin nodisk t=24 minid=1 \
				covstats=covstats_1.txt outm=${bin}_mapped.fastq.gz;
		spades.py -k 21,33,55,77,99 -m 400 -t 24 -o ${bin}_spades --12 ${bin}_mapped.fastq.gz;
		cd ${bin}_spades/;
		barrnap  -o 16s_rrna.fasta contigs.fasta;
	done

#END