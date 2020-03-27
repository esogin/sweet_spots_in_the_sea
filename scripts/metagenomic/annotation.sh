#Bin annotation
#EM Sogin
#Seagrass rhizosphere microbiome
#Updated March 2020

#Call genes

for i in *fa; do 
	prodigal -i $i -o coords/${i%%.fa}_coords.gbk  \
			 -a proteins/${i%%.fa}_orfs.faa \
			 -d genes/${i%%.fa}_genes.fa;
done

#Annotate proteins
for Proteins in $(ls *.faa | sed 's/\.faa//g'); 
	do
	#Diamond search against cutsom uniprot download
		diamond blastp \
		    --tmpdir /tmp/ \
		    -p 20 \
		    -q "$Proteins".faa \
		    --db /scratch/dmeier/Uniprot.dmnd \
		    --daa "$Proteins".daa \
		    --evalue 1E-15 \
		    --salltitles \
		    -k 1;
	#Generating a tabular Output from the diamond alignment file, specifying the columns for the output at the end.
	#The column "salltitles" includes protein function and organism name
		diamond view \
		    --daa "$Proteins".daa \
		    --out "$Proteins"_diamond_Uniprot.tsv \
		    -f 6 qseqid sseqid salltitles pident evalue

	#Emapper search 
		emapper.py \
		    --cpu 1 \
		    -m diamond \
		    --data_dir /apps/eggnogmapper/2.0.0/data/ \
		    --dmnd_db /apps/eggnogmapper/2.0.0/data/eggnog_proteins.dmnd \
		    --qtype seq \
		    --seed_ortholog_evalue 0.00001 \
		    -i "$Proteins".faa \
		    --output_dir ./Emapper_Diamond \
		    -o "$Proteins"
	#Hmmscan serach agains dbCan V7 database
		hmmscan \
		    --cpu 1 \
		    -E 1E-15 \
		    --domtblout dbCAN_annotations/"$Proteins".out \
		    --tblout dbCAN_annotations/"$Proteins".out \
		    dbCAN-HMMdb-V7.hmm \
		   	"$Proteins".faa
    #InterPro scan 
	    /Path/to/interpro/interproscan-5.36-75.0/interproscan.sh \
		    -i "$Proteins".faa \
		    -f tsv \
		    -o "$Proteins".csv \
		    -appl Pfam,Phobius \
		    -T /tmp/
   done

 #END  