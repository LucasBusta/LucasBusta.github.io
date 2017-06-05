#!/bin/bash

#OIFS="$IFS"
#IFS=$'\n'
version=3.0
bold=$(tput bold)
normal=$(tput sgr0)

#Usage
	Usage="

  \\    /\\
   )  ( ') TRANSCriptome Analysis Toolbox  
  (  /  )   (TRANSCAT) Version $version
   \\(__)|

	Usage: 
		bash transcat.sh -T <filename> [-Q <filename>] [instructions]
		example: bash transcat.sh -T N_benth_transcriptome -Q FAD3s -cbgst
	
	Options:
		Required:
	    -T <filename> Specify ${bold}TRANSCRIPTOME${normal} (must be nucleotide .fasta, for example, a list of contigs)

	    Optional:
	    -Q <filename> Specify ${bold}QUERY${normal} (must be protein .fasta, for example, target protein sequence(s))

	    Instructions:
	    -c 			  ${bold}CLEAN${normal} the TRANSCRIPTOME (i.e. simplify names, linearize etc.)
	    -g 			  ${bold}GET${normal} DATA from the SRA
	    -s 			  Run ${bold}SALMON${normal} on READS v. REFERENCE
	    -b 			  ${bold}BLAST${normal} the QUERY against the TRANSCRIPTOME
	    -t 			  Genterate ${bold}HEAT TREE${normal} using blast results (prerequisites: ldbt)

	    Help:

	    --help 		  Displays this helpful help screen
	 
	 Required steps:
		1. Create directories: 'reference', and 'scripts'
		2. Place this script in 'scripts'
		3. Place TRANSCRIPTOME (as nucleotide .fasta) in 'reference'
		
	 Steps for additional functions:
	 	1. Place QUERY protein sequence (as .fasta) in 'reference'
	 	2. Place list of READS (as 'reads') in 'reference' interleaved with sample NAMES
	 		example:	SRR685298
	 				 	Leaf
	 				 	SRR696884
	 				 	Stem
	 				 	SRR696915
	 				 	Flower"

#Set options
	T=0
	Q=0
	c=0
	g=0
	s=0
	b=0
	t=0

	optarg=0

	#colons come after options that require arguments
	shrt_opts=T:Q:cgsbt
	long_opts=help

	opts=`getopt -n transcat.sh -o $shrt_opts -ul $long_opts -- $@`

	eval set -- "$opts"

	while [ $1 != -- ]; do
		case "$1" in
		-T ) T=$2;;
		-Q ) Q=$2;;
		-c ) c=1;;
		-g ) g=1;;
		-s ) s=1;;
		-b ) b=1;;
		-t ) t=1;;
		--help)
			echo ""
			echo "$Usage"
			echo ""
			exit 0
			;;
		esac
		shift
	done
	shift	# move past '--' from getopt call

#Welcome screen
	
	echo "
  \\    /\\
   )  ( ') TRANSCriptome Analysis Toolbox  
  (  /  )   (TRANSCAT) Version $version
   \\(__)|"
	echo ""

	#Identify the TRANSRIPTOME
	echo TRANSCRIPTOME = "../reference/${T}.fasta"
	if [ -e "../reference/${T}.fasta" ]; then
	    echo "Found ../reference/${T}.fasta"
		else
		echo "Cannot find ../reference/${T}.fasta"
	fi

	#Identify the QUERY
	if [ "$Q" != "0" ]; then
		echo QUERY = "../reference/${Q}.fasta"
		if [ -e "../reference/${Q}.fasta" ]; then
		    echo "Found ../reference/${Q}.fasta"
			else
			echo "Cannot find ../reference/${Q}.fasta"
		fi
	fi
	echo ""

#[-c]CLEAN the TRANSCRIPTOME
	if [ "$c" == "1" ]; then

		#simplify TRANSCRIPTOME  headers
			echo "Cleaning ../reference/${T}.fasta:"
			# 	echo "Simplifying ../reference/${T}.fasta headers..."
			# 	#cutting off up 'til "LOC":  (^ means beginning of line, .* is wildcard)
			# 		# sed -i 's/^.*LOC/>LOC/' rna.fasta
			# 	#cutting after ", mRNA":  ($ means end of line)
			# 		# sed -i 's/ .*$//' rna.fasta
			# 	#replace all spaces:
			# 		sed -i 's/ /_/g' ../reference/${T}.fasta
			# 	#replace all |:
			# 		sed -i 's/|/_/g' ../reference/${T}.fasta
	fi

#[-b]BLAST the QUERY against the TRANSCRIPTOME
	if [ "$b" == "1" ]; then
		if [ "$Q" != "0" ]; then

		#linearizing the TRANSCRIPTOME
			echo "Linearizing ../reference/${T}.fasta..."
			sed -e 's/\(^>.*$\)/#\1#/' "../reference/${T}.fasta" | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > "../reference/${T}_l.fasta"
			rm ../reference/${T}.fasta
			mv ../reference/${T}_l.fasta ../reference/${T}.fasta

		#Make the transcriptome BLASTable
			echo "Making BLASTable database from ../reference/${T}.fasta..."
			module load blast
			makeblastdb -in "../reference/${T}.fasta" -parse_seqids -dbtype nucl
			
		#Linearize the QUERY
			echo ""
			echo "Linearizing ../reference/${Q}.fasta..."
			sed -e 's/\(^>.*$\)/#\1#/' "../reference/${Q}.fasta" | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > "../reference/${Q}_l.fasta"
			rm ../reference/${Q}.fasta
			mv ../reference/${Q}_l.fasta ../reference/${Q}.fasta

		#Blast query against database
			echo "BLASTing ../reference/${Q}.fasta against ../reference/${T}.fasta..."
			tblastn -db "../reference/${T}.fasta" -query "../reference/${Q}.fasta" -out "../reference/${Q}_${T}_hits.out" -outfmt "6 sacc"
			echo "Removing duplicate hits..."
			awk '!a[$0]++' "../reference/${Q}_${T}_hits.out" > "../reference/${Q}_${T}_hits_s.out"
			rm ../reference/${Q}_${T}_hits.out
			mv ../reference/${Q}_${T}_hits_s.out ../reference/${Q}_${T}_hits.out
			nhits=$(cat "../reference/${Q}_${T}_hits.out" | wc -l)
			echo "$nhits Hits"
			echo "Preparing output file..."

		#Create "results.fasta" (output fasta of BLAST hits) contents
			for i in $(seq 1 $nhits); do
				hit_name=$(sed -n "${i}p" "../reference/${Q}_${T}_hits.out");
				hit_line=$(grep -n $hit_name "../reference/${T}.fasta"| cut -d : -f 1);
				#echo "Hit Line: $hit_line"
				seq_line=$((1+$hit_line));
				#echo "Seq Lin: $seq_line"
				hit_header=$(sed -n "${hit_line}p" "../reference/${T}.fasta");
				#echo "Hit: $hit_name"
				seq=$(sed -n "${seq_line}p" "../reference/${T}.fasta");
				#echo "Seq: $seq"
				hits_fasta[$((2 * ${i} - 1))]=$hit_header;
				hits_fasta[$((2 * ${i}))]=$seq;
			done

		#print "results.fasta"
			printf '%s\n' "${hits_fasta[@]}" > "../reference/${Q}_${T}_hits.fasta"
			echo "Filtering out hits less than 800bp..."
			awk '!/^>/ { next } { getline seq } length(seq) >= 800 { print $0 "\n" seq }' "../reference/${Q}_${T}_hits.fasta" > "../reference/${Q}_${T}_hits_f.fasta"
			rm ../reference/${Q}_${T}_hits.fasta
			mv ../reference/${Q}_${T}_hits_f.fasta ../reference/${Q}_${T}_hits.fasta
			echo "Done!"

		else
			echo "Please specify a BLAST query."
		fi
	fi
	echo ""

#[-g]Get READS from SRA
	if [ "$g" == "1" ]; then

		#Create fastq directory
			if [ -d "../fastq" ]; then
			    echo "Found directory '../fastq'"
				else
				mkdir ../fastq
				echo "Created directory '../fastq'"
			fi

		#Create reads variable
			reads=$(cat ../reference/reads)
			echo "Files to get: $reads"
			
		#Create and dispatch fqdump.slurm for each DATA
			if [ -d "LOG" ]; then
			    echo "Found directory 'LOG'"
				else
				mkdir LOG
				echo "Created directory 'LOG'"
			fi
			echo "#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=SAMPLE
#SBATCH --error=./LOG/get_SAMPLE.err
#SBATCH --output=./LOG/get_SAMPLE.out

module load SRAtoolkit

cd ../fastq
fastq-dump -A SAMPLE --split-files" > fqdump.slurm

			for x in $reads
				do
				sed "s/SAMPLE/$x/g" fqdump.slurm > fqdump_auto.slurm
				sbatch fqdump_auto.slurm
			done
	fi

# #rename files according to read_names
# 		echo "Renaming files..."
# 		cd ../fastq
# 		for FILE in `ls *.fastq | sort`; do
# 			R=$(sed "s/.fastq//" <<< "$FILE");
# 			#echo $R;
# 		    m=$(grep -n $R ../reference/read_names | cut -f1 -d:);
# 		    #echo $m;
# 		    n=$(($m+1));
# 		    #echo $n;
# 		    #echo Moving $FILE to `sed -n ${n}p ../reference/read_names`;
# 		    mv "${FILE}" `sed -n ${n}p ../reference/read_names`;
# 		done
# 			#add .fastq to all file names
# 			find . -type f -exec mv '{}' '{}'.fastq \;
# 		echo "Done."
# 		cd ../scripts

#[-s]run salmon
	if [ "$s" == "1" ]; then

		# #rename files according to read_names
		# echo "Renaming files..."
		# for FILE in `ls *.fastq | sort`; do
		# 	R=$(sed "s/.fastq//" <<< "$FILE");
		# 	#echo $R;
		#     m=$(grep -n $R ../reference/read_names | cut -f1 -d:);
		#     #echo $m;
		#     n=$(($m+1));
		#     #echo $n;
		#     #echo Moving $FILE to `sed -n ${n}p ../reference/read_names`;
		#     mv "${FILE}" `sed -n ${n}p ../reference/read_names`;
		# done
		# 	#add .fastq to all file names
		# 	find . -type f -exec mv '{}' '{}'.fastq \;

		#Build salmon directories
			if [ -d "../salmon" ]; then
			    rm -r ../salmon
			    mkdir ../salmon
			    mkdir ../salmon/${T}
			    echo "Overwriting directory '../salmon'"
				else
				mkdir ../salmon
				mkdir ../salmon/${T}
				echo "Created directory '../salmon'"
				echo "Created directory '../salmon/${T}'"
			fi

		#Build salmon reference
			echo "Building salmon database for ../reference/${T}..."
			module load salmon
			salmon index -t ../reference/${T}.fasta -i ../reference/${T}_index -k 31

		#Dispatch salmon runs
			echo "Dispatching salmon quantification runs..."
			reads=$(cat ../reference/reads)
			echo "#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH --mem=62GB
#SBATCH --job-name=salmon_${T}_SAMPLE
#SBATCH --error=./LOG/salmon_${T}_SAMPLE.err
#SBATCH --output=./LOG/salmon_${T}_SAMPLE.out

module load salmon
module load erne
#erne-filter --query1 ../fastq/SAMPLE_1.fastq --query2 ../fastq/SAMPLE_1.fastq --output-prefix ../fastq/SAMPLE_f --ultra-sensitive --force-standard
#salmon quant -i ../reference/${T}_index -l A -1 ../fastq/SAMPLE_f_1.fastq -2 ../fastq/SAMPLE_f_2.fastq -o ../salmon/${T}/SAMPLE
#salmon quant -i ../reference/${T}_index -l A -1 ../fastq/SAMPLE_1.fastq -2 ../fastq/SAMPLE_2.fastq -o ../salmon/${T}/SAMPLE
salmon quant -i ../reference/${T}_index -l A -r ../fastq/SAMPLE.fastq -o ../salmon/${T}/SAMPLE
" > "salmon_${T}.slurm"

		reads=$(sed -n '1~2!p' ../reference/reads)
		for x in $reads
			do
			sed "s/SAMPLE/$x/g" "salmon_${T}.slurm" > "salmon_${T}_auto.slurm"
			sbatch "salmon_${T}_auto.slurm" 
		done

		echo "Done."
	fi

#[-t]generate tree
	if [ "$t" == "1" ]; then

		#Merge salmon results, check % mapped
			reads=$(sed -n '1~2!p' ../reference/reads)
			if [ -d ../salmon/${T}/_merge ]; then
				rm -r ../salmon/${T}/_merge
				mkdir ../salmon/${T}/_merge
				echo "Overwriting ../salmon/${T}/_merge"
				else
				mkdir ../salmon/${T}/_merge
				echo "Created ../salmon/${T}/_merge"
			fi	
			for i in $reads
				do
				cp ../salmon/${T}/${i}/quant.sf ../salmon/${T}/_merge/${i}.sf
			done
			if [ -e rate ]; then
				rm rate
			fi
			touch rate
			for i in $reads
				do 
				grep "Mapping rate" LOG/salmon_${T}_${i}.err | sed 's/^.* = //' | sed 's/%//' >> rate
			done
			avg=$(awk '{x+=$1;next}END{print x/NR}' rate)
			rm rate
			echo "Average mapping rate = ${avg} %"

 		#Translate ${Q}_${T}_hits.fasta to protein
			cd ../reference
				echo "Translating ../reference/${Q}_${T}_hits.fasta to peptide..."
				if [ -d ${Q}_${T}_hits.fasta.transdecoder_dir ]; then
					rm -r ${Q}_${T}_hits.fasta.transdecoder_dir
					echo "Overwriting ../reference/${Q}_${T}_hits.fasta.transdecoder_dir"
				fi
				if [ -e ${Q}_${T}_hits.fasta.transdecoder.pep ]; then
					rm ${Q}_${T}_hits.fasta.transdecoder*
					echo "Overwriting ../reference/${Q}_${T}_hits.fasta.transdecoder*"
				fi
				module load transdecoder
				TransDecoder.LongOrfs -t ${Q}_${T}_hits.fasta
				TransDecoder.Predict --retain_long_orfs 340 --single_best_orf -t ${Q}_${T}_hits.fasta

			#Clean up ${Q}_${T}_hits.fasta.transdecoder.pep
				echo "Cleaning up ${Q}_${T}_hits.fasta.transdecoder.pep..."
					echo "Linearizing ${Q}_${T}_hits.fasta.transdecoder.pep..."
					sed -e 's/\(^>.*$\)/#\1#/' "${Q}_${T}_hits.fasta.transdecoder.pep" | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > "${Q}_${T}_l_hits.fasta.transdecoder.pep"
					rm ${Q}_${T}_hits.fasta.transdecoder.pep
					mv ${Q}_${T}_l_hits.fasta.transdecoder.pep ${Q}_${T}_hits.fasta.transdecoder.pep
				#cutting off up 'til "LOC":  (^ means beginning of line, .* is wildcard)
					sed -i 's/^.*(+) />/' ${Q}_${T}_hits.fasta.transdecoder.pep
				#cutting after ", mRNA":  ($ means end of line)
					sed -i 's/:/_/' ${Q}_${T}_hits.fasta.transdecoder.pep
					sed -i 's/_.*$//' ${Q}_${T}_hits.fasta.transdecoder.pep
				#remove * (stop codon characters) from pep.fasta
 					sed -i 's/*//g' ${Q}_${T}_hits.fasta.transdecoder.pep
					cp ${Q}_${T}_hits.fasta.transdecoder.pep ${Q}_${T}_hits_pep.fasta
			cd ../scripts

 		#concatenate the BLAST results and QUERY
			echo "Merging and linearizing ../reference/${Q}_${T}_hits_pep.fasta and ../reference/${Q}.fasta..."
			awk '{print}' "../reference/${Q}_${T}_hits_pep.fasta" "../reference/${Q}.fasta" > "../reference/${Q}_${T}_hqmerge.fasta"
			sed -e 's/\(^>.*$\)/#\1#/' "../reference/${Q}_${T}_hqmerge.fasta" | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > "../reference/${Q}_${T}_hqmerge_l.fasta"
			rm ../reference/${Q}_${T}_hqmerge.fasta
			mv ../reference/${Q}_${T}_hqmerge_l.fasta ../reference/${Q}_${T}_hqmerge.fasta

		#size filter ${Q}_${T}_hqmerge.fasta
			echo "Length filtering 400aa to 600aa..."
			awk '!/^>/ { next } { getline seq } length(seq) >= 400 { print $0 "\n" seq }' "../reference/${Q}_${T}_hqmerge.fasta" > "../reference/${Q}_${T}_hqmerge_f.fasta"
			awk '!/^>/ { next } { getline seq } length(seq) <= 600 { print $0 "\n" seq }' "../reference/${Q}_${T}_hqmerge.fasta" > "../reference/${Q}_${T}_hqmerge_ff.fasta"
			rm ../reference/${Q}_${T}_hqmerge.fasta ../reference/${Q}_${T}_hqmerge_f.fasta
			mv ../reference/${Q}_${T}_hqmerge_ff.fasta ../reference/${Q}_${T}_hqmerge.fasta

 		#Align ${Q}_${T}_hqmerge.fasta
			echo "Aligning ../reference/${Q}_${T}_hqmerge.fasta..."
			module load clustal-omega
			clustalo -i ../reference/${Q}_${T}_hqmerge.fasta -o ../reference/${Q}_${T}_hqmerge_align.fasta --force

		#Create NJ heat tree
			echo "Creating neighbor-joining heat tree..."
			cd "../salmon/${T}"
			cd _merge 
			merdir=`pwd`
			cd ..
			cd ..
			cd ../reference
			refdir=`pwd`
			cd /home/moriyama/lucasbusta/R

			module load R
			export R_LIBS="/home/moriyama/lucasbusta/R/R_libs"
			echo "library(DESeq2)
library(tximport)
library(readr)
library(tximportData)
library(ggtree)
library(ape)
library(msa)
library(phangorn)

#Create protein fasta phylo
	setwd("\""$refdir"\"")
	aa <- read.aa("\""${Q}_${T}_hqmerge_align.fasta"\"", format="\""fasta"\"")
	phylo <- NJ(dist.ml(aa))
	dim(as.data.frame(phylo))
	as.data.frame(phylo)

#Create sample Table
	setwd("\""$merdir"\"")
	files <- dir()
	names(files) <- substr(files,1,10)
	sampleTable <- data.frame(files)
	sampleTable

#Create annotation file
	setwd("\""$merdir"\"")
	txi.salmon <- tximport(files, type = "\""salmon"\"", txOut=TRUE)
	dds1 <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~1)
	dds2 <- DESeq2::fpkm(dds1, robust=FALSE) #NOTE THAT ROBUST IS AFFECTED BY N SAMPLES
	#dds2 <- log2(dds2)
	head(dds2)

	ngenes <- as.numeric(table(as.data.frame(phylo)\$isTip)["\""TRUE"\""])
	rnas <- as.data.frame(phylo)\$label[1:ngenes]
	rnas
	
	keep <- which(rownames(dds2) %in% rnas)
	annot <- dds2[keep,]
	dim(annot)
	annot

#plot protein tree
	setwd("\""$refdir"\"")
	p <- ggtree(phylo, size=1, branch.length = "\""none"\"", ladderize=TRUE)+
		geom_point(size=1)+
		theme(legend.position="\""right"\"")+
		geom_tiplab(align=TRUE, linesize=.5, size=5, hjust=0, offset=0.1)

	pdf(file="\""ggtreeP.pdf"\"",width=16,height=16)
	gheatmap(p, annot, width=3, offset=6, low = "\""#ffffff"\"", high = "\""#0091d1"\"", color="\""black"\"", colnames_position="\""top"\"", colnames_angle="\""90"\"", colnames_offset_x=0, colnames_offset_y=0, font.size=5, hjust=0)+
		theme(legend.position="\""right"\"", legend.key.size = unit(1, "\""cm"\""), legend.text=element_text(size=20))
	dev.off()" > ggtree.R

				R CMD BATCH ggtree.R
			cd "$refdir"
			echo "Done!"
	fi
		
#IFS="$OIFS"
