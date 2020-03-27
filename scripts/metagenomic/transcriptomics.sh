#Transcritome mapping
#EM Sogin
#Seagrass rhizosphere microbiome
#Updated March 2020

ref=/path/to/illumina/adapters.fa;
sortme=/path/to/sortmerna-2.1b

#Clean transcriptomic lirabries, quality trim, error correct, filter rRNA
for lib in $libraries;
	do
		bbduk.sh ref=${ref} ktrim=l mink=11 hdist=1 in="$lib"_R1.fastq.gz \
				 in2="$lib"_R1.fastq.gz out="${lib}"_ktriml.fq.gz;
		bbduk.sh ref=${ref} ktrim=r trimq=2 qtrim=rl minlength=50 mink=11 hdist=1 \
				 in="${lib}"_ktriml.fq.gz out="${lib}"_q2_ktrimmed.fq.gz;
		spades.py -o "${lib}"_spades --12 "${lib}"_q2_ktrimmed.fq.gz \
				  -t 48 --phred-offset 33 --only-error-correction;
		reformat.sh in1="${lib}"_spades/corrected/"${lib}"_q2_ktrimmed_1.00.0_0.cor.fastq.gz \
				    in2="${lib}"_spades/corrected/"${lib}"_q2_ktrimmed_2.00.0_0.cor.fastq.gz \
				    out="${lib}"_corrected.fastq.gz
		gunzip "${lib}"_corrected.fastq.gz
		${sortme}/sortmerna --ref \
			${sortme}/rRNA_databases/silva-bac-16s-id90.fasta,${sortme}/index/silva-bac-16s-db:\
			${sortme}/rRNA_databases/silva-bac-23s-id98.fasta,${sortme}/index/silva-bac-23s-db:\
			${sortme}/rRNA_databases/silva-arc-16s-id95.fasta,${sortme}/index/silva-arc-16s-db:\
			${sortme}/rRNA_databases/silva-arc-23s-id98.fasta,${sortme}/index/silva-arc-23s-db:\
			${sortme}/rRNA_databases/silva-euk-28s-id98.fasta,${sortme}/index/silva-euk-28s:\
			${sortme}/rRNA_databases/rfam-5s-database-id98.fasta,${sortme}/index/rfam-5s-db:\
			${sortme}/rRNA_databases/rfam-5.8s-database-id98.fasta,${sortme}/index/rfam-5.8s-db \
			--reads "${lib}"_corrected.fastq --num_alignments 1 \
			--fastx --aligned "${lib}"_reads_rRNA --other "${lib}"_reads_non_rRNA --log -a 24 -m 64000 --paired_in -v
		gzip *fastq
	done

#Prep high QC MAGs for mapping
#rename gene file fasta headers  to include name of mag in header
bins=$(echo *.contigs_genes.fa | sed 's/\.contigs_genes.fa//g')
for MAG in $bins; 
    do
        sed -i "s/k151/$MAG/g" "$MAG".contigs_genes.fa 
    done
cat *fa > all_genes_DNA_MAP.fa

#Count transcripts against all high QC bins in analysis
ref=all_genes_DNA_MAP.fa
for lib in $libraries;
	do;
		#Split iterleaved reads
			reformat.sh in="${lib}"_reads_non_rRNA.fastq.gz \
						out1="${lib}"_non_rRNA_R1.fq.gz out2="${lib}"_non_rRNA_R2.fq.gz
			r1="${lib}"_non_rRNA_R1.fq.gz
			r2="${lib}"_non_rRNA_R2.fq.gz
			kallisto index -i ref.idx $ref
			kallisto quant -i ref.idx -o "$lib"_output -b 100 $r1 $r2
	done
#END