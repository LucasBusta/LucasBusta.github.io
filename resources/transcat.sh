#!/usr/bin/perl
#🜁air🜄water🜂fire🜃earth
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Data::Dumper;
use File::Spec::Functions qw/ splitdir rel2abs /;
use POSIX;

my $version = 0.1;

my $M = 62;
my $C = 50;
my $readLength;

my @hits_fasta;
my @reads;
my @names;
my $num;
my $i;
my $p;
my $mergedir;
my $miningdir;
my $filterIN;
my $filterMID;
my $filterOUT;
my @short;
my $anchor_name = "";
my ($IN, $IN2, $OUT, $OUT2);
my $line;
my $nReads;
my $nBases;
my $avgLength;
my @readlength;
my $kMin = 15;
my $kMax;
my $read;
my $reads;
my $file;
my %seen;

my $mainDir;
my $baseDir;

my $project;

my $T;
my $Q;

my $help = 0;
my $h = 0;

my $single;
my $paired;
my $get;
my $preanalyze;

my @assemblers = ("idba", "soap", "spades", "trinity");
my $q = 20;
my $filter;
my $normalize;
my $merge;
my $skipmergenorm;
my $mergenorm;

my $assemble;
my $idba;
my $soap;
my $spades;
my $trinity;
my $kmerge;
my $translatedenovos;
my $nonredund;
my $consensus;

my $translate;

my $anchor;
my $sindex;
my $squant;
my $blast;
my $plot;
my $quant;

my $collect;
my $subprojects;

my $height = 28;
my $width = 12;
my $ylim = 10;

GetOptions(	"help" => \$help,
			"h" => \$h,
			"version" => \$version,

			"project=s" => \$project,

			"single=s" => \$single,
			"paired=s" => \$paired,
			"get" => \$get,
			"preanalyze" => \$preanalyze,

			"M=i" => \$M,
			"T=s" => \$T,
			"Q=s" => \$Q,
			"C=s" => \$C,

			"filter" => \$filter,
			"normalize" => \$normalize,
			"merge" => \$merge,
			"skipmergenorm" => \$skipmergenorm,
			"mergenorm" => \$mergenorm,

			"assemble" => \$assemble,
			"idba" => \$idba,
			"spades" => \$spades,
			"soap" => \$soap,
			"trinity" => \$trinity,
			"kmerge" => \$kmerge,
			"translatedenovos" => \$translatedenovos,
			"nonredund" => \$nonredund,
			"consensus" => \$consensus,

			"translate=s" => \$translate,

			"anchor=s" => \$anchor,
			"sindex" => \$sindex,
			"squant" => \$squant,
			"blast=i" => \$blast,
			"plot" => \$plot,
			"quant" => \$quant,
			"height=i" => \$height,
			"width=i" => \$width,
			"ylim=i" => \$ylim,

			"collect" => \$collect,
			"subprojects=s" => \$subprojects,
			"reads=s" => \$reads,
);

print "\n 				🜁🜄  Elemental v${version}  🜃🜂\n\n";

if ($h or $help) {
	print "
	Conventions:
		-> 	Takes input files with .fa and .fq extensions.
		-> 	Assumes {name}_1.xx and {name}_2.xx nomenclature for 
			paired reads. Provide only {name} in configuration file.

	Usage: 	perl elemental.pl [functions]
 
		1. 	Create project directory with subdirectory 'scripts' 
			containing a readsfile or projectfile
				for example, a file 'ivy' with contents: 
					SRR685298
					Leaf_1
					SRR696884
					Leaf_2
					SRR696915
					Flower

		2. 	If using reads you already have, create a RawReads 
			folder containing, for example (if paired):
				SRR685298(_1).fq
				(SRR685298_2.fq)
				SRR696884(_1).fq
				(SRR696884_2.fq), etc...

			Alternatively, use -get to download SRA reads (see below)

		3. 	Read pre-processing and assembly functions:
				[provide: -single <readfile> or -paired <readfile>]:
			
				#Read pre-processesing
					-get 			Download reads from the SRA
					-filter 		Filter the rawreads 	[-M 62 can do reads at least 16GB (sum for paired), maybe more, but not 34GB]
				
				#Pre-processesing for assembly
					-normalize 		Normalize filtered reads 	[-M 62 {0-1 GB} -M 250 {1+ GB}]
					-mergenorm 		Merge specified reads and normalize them to prepare for assembly 	[-M 62 {0-5 GB} -M 250 {5-40+ GB}]
					-skipmergenorm 		(if only one set of readeads is being used then mergenorm is not necessary)
				
				#Assembly
					-assemble 		Assemble normalized reads using the specified assemblers

				#Post-assembly
					-kmerge
					-translatedenovos 	Creates four files, each containing all contigs from one assembler, then translates them
					-nonredund 		Selects nonredund AA-generating contigs from each assembler, creating 'nonredund_{assembler}_contigs.fa
					-consensus 		Selects contigs found by a minimum of <integer> assemblers

		4.	Transcriptome analyis functions:
				[provide -T <filename> (.fa) and -Q <filename> (.faa)]
			
				#Transcriptome cleaning
					truncate header 'til  (^ means beginning of line, .* is wildcard)
						sed -i 's/^.*|X/>X/' rna.fa
					truncate header after ::  (\$ means end of line)
						sed -i 's/|_.*\$//' MAH1_rna_hqmerge_align.faa
					replace all spaces:
						sed -i 's/ /_/g' <file>
					replace all |:
						sed -i 's/|/_/g' <file>
					truncation series:
						sed -i -e 's/^.*XR_/>XR_/' -e 's/^.*XM_/>XM_/' -e 's/| .*\$//' rna.fasta
					add _ to the end of the header
						sed 's/>.*/&_/' ivyConsensus.fasta > ivyConsensus_.fasta

				#Read quantification
					-sindex (-T)	Index a transcriptome to prepare for salmon-based quantification
					-squant (-T)	-single or -paired	Quantify transcripts in a reads file\n\n";
}

if ($project) {

	print "Setting up Elemental...\n";

		# print "Be sure to initialize modules using \"source modules.sh\" before running this script!\n";
			open ($OUT, ">modules.sh") or die "Can't open modules.sh";
			print $OUT "#!/bin/bash\nml R blast salmon genemarks/4.3 R clustal-omega\nexport R_LIBS=/home/moriyama/lucasbusta/R/R_libs\n";
			close $OUT;

		# print "Defining base directory and project name...\n";
			$mainDir = "/work/moriyama/lucasbusta/";
			$baseDir = "/work/moriyama/lucasbusta/${project}";

			# chdir ('..') or die "Failed to get to base directory.";
			# my $baseDir = getcwd;
			# chomp $baseDir;
			# print "	Base directory = $baseDir\n";
			# my $project = (splitdir(rel2abs($0)))[-2];
			# chomp $project;
			# print "	Project name = $project\n\n";
			# chdir ('scripts') or die "Check directory structure. Elemental should be in $baseDir/scripts\n\n";

		# print "Setting up directories...\n";
			if (-d "$baseDir/scripts/LOG") {
				print "	Found $baseDir/scripts/LOG\n";
			} else {
				mkdir ("$baseDir/scripts/LOG");
				print "	Made $baseDir/scripts/LOG\n";
			}
			if (-d "$baseDir/RawReads") {
				print "	Found $baseDir/RawReads\n";
			} else {
				mkdir ("$baseDir/RawReads");
				print "	Made $baseDir/RawReads\n";
			}
			if (-d "$baseDir/FilteredReads") {
				print "	Found $baseDir/FilteredReads\n";
			} else {
				mkdir ("$baseDir/FilteredReads");
				print "	Made $baseDir/FilteredReads\n";
			}
			if (-d "$baseDir/NormalizedReads") {
				print "	Found $baseDir/NormalizedReads\n";
			} else {
				mkdir ("$baseDir/NormalizedReads");
				print "	Made $baseDir/NormalizedReads\n";
			}
			if (-d "$baseDir/NormalizedMergedReads") {
				print "	Found $baseDir/NormalizedMergedReads\n";
			} else {
				mkdir ("$baseDir/NormalizedMergedReads");
				print "	Made $baseDir/NormalizedMergedReads\n";
			}
			if (-d "$baseDir/Assemblies") {
				print "	Found $baseDir/Assemblies\n";
			} else {
				mkdir ("$baseDir/Assemblies");
				print "	Made $baseDir/Assemblies\n";
			}
			if (-d "$baseDir/Assemblies/idba") {
				print "	Found $baseDir/Assemblies/idba\n";
			} else {
				mkdir ("$baseDir/Assemblies/idba");
				print "	Made $baseDir/Assemblies/idba\n";
			}
			if (-d "$baseDir/Assemblies/soap") {
				print "	Found $baseDir/Assemblies/soap\n";
			} else {
				mkdir ("$baseDir/Assemblies/soap");
				print "	Made $baseDir/Assemblies/soap\n";
			}
			if (-d "$baseDir/Assemblies/spades") {
				print "	Found $baseDir/Assemblies/spades\n";
			} else {
				mkdir ("$baseDir/Assemblies/spades");
				print "	Made $baseDir/Assemblies/spades\n";
			}
			if (-d "$baseDir/Assemblies/trinity") {
				print "	Found $baseDir/Assemblies/trinity\n";
			} else {
				mkdir ("$baseDir/Assemblies/trinity");
				print "	Made $baseDir/Assemblies/trinity\n";
			}
			if (-d "$baseDir/Mining") {
				print "	Found $baseDir/Mining\n";
			} else {
				mkdir ("$baseDir/Mining");
				print "	Made $baseDir/Mining\n";
			}
			if (-d "$baseDir/Mining/salmon") {
				print "	Found $baseDir/Mining/salmon\n";
			} else {
				mkdir ("$baseDir/Mining/salmon");
				print "	Made $baseDir/Mining/salmon\n";
			}
} else {
	print "Please specify a project...\n";
}


if ($single) {

	print "Setting up reads...\n";
		if (-e "$baseDir/scripts/$single") { 
			print "	Found $baseDir/scripts/$single\n"; 
			} else { 
			die "	Can't find $baseDir/scripts/$single\n";
		}

		print "	Importing read list...\n";
			my @reads = `sed -n '0~2!p' $baseDir/scripts/$single`;
			chomp @reads;
			print "	Reads are: @reads\n\n";

	if ($get) {
		print "Submitting jobs to get SRA files...\n";
		foreach my $read (@reads) {
			open ($OUT, ">$baseDir/scripts/fastqdump_${read}.slurm") or die "Can't open $baseDir/scripts/fastqdump_${read}.slurm";
			print $OUT "#!/bin/sh\n#SBATCH --time=12:00:00\n#SBATCH --mem=16GB\n#SBATCH --job-name=${project}_${read}_fqd\n#SBATCH --error=$baseDir/scripts/LOG/fqd_${read}.err\n#SBATCH --output=$baseDir/scripts/LOG/fqd_${read}.out
			cd $baseDir/RawReads
			module load SRAtoolkit
			fastq-dump ${read}
			mv $baseDir/RawReads/${read}.fastq $baseDir/RawReads/${read}.fq";
			# prefetch ${read}
			close $OUT;
			system ("sbatch $baseDir/scripts/fastqdump_${read}.slurm");
		}
	}

	if ($preanalyze) {
		print "Submitting fastqc job(s)...\n";
		foreach my $read (@reads) {
			open ($OUT, ">$baseDir/scripts/fastqc_${read}.slurm") or die "Can't open $baseDir/scripts/fastqdump_${read}.slurm";
			print $OUT "#!/bin/sh\n#SBATCH --time=12:00:00\n#SBATCH --mem=16GB\n#SBATCH --job-name=${project}_${read}_fqc\n#SBATCH --error=$baseDir/scripts/LOG/fqc_${read}.err\n#SBATCH --output=$baseDir/scripts/LOG/fqc_${read}.out
			cd $baseDir/RawReads
			module load fastqc
			fastqc ${read}.fq
			rm ${read}_fastqc.zip";
			close $OUT;
			system ("sbatch $baseDir/scripts/fastqc_${read}.slurm");
		}
	}

	if ($filter) {		
		print "Dispatching filtering job(s) for single reads at q = $q...\n";
		foreach my $read (@reads) {
			open ($OUT, ">$baseDir/scripts/filter_${read}_q${q}.slurm");
			print $OUT "#!/bin/sh\n#SBATCH --time=12:00:00\n#SBATCH --mem=${M}GB\n#SBATCH --job-name=${project}_${read}_filter_q$q\n#SBATCH --error=$baseDir/scripts/LOG/filter_${read}_q$q.err\n#SBATCH --output=$baseDir/scripts/LOG/filter_${read}_q$q.out\n
			module load cutadapt
			module load erne
			cutadapt -l 600 $baseDir/RawReads/${read}.fq > $baseDir/FilteredReads/${read}_short.fq
			erne-filter --query1 $baseDir/FilteredReads/${read}_short.fq --output-prefix $baseDir/FilteredReads/${read} --ultra-sensitive --force-standard
			mv $baseDir/FilteredReads/${read}_1.fastq $baseDir/FilteredReads/${read}.fq
			rm $baseDir/FilteredReads/${read}_short.fq";
			close $OUT;
			system("sbatch $baseDir/scripts/filter_${read}_q${q}.slurm");
		}
	}

	if ($normalize) {
		print "Dispatching normalization job(s) for single reads...\n";
		foreach my $read (@reads) {
			open ($OUT, ">$baseDir/scripts/normalize_${read}.slurm");
			print $OUT "#!/bin/sh\n#SBATCH --time=168:00:00\n#SBATCH --mem=${M}GB\n#SBATCH --job-name=${project}_${read}_normalize\n#SBATCH --error=$baseDir/scripts/LOG/normalize_${read}.err\n#SBATCH --output=$baseDir/scripts/LOG/normalize_${read}.out\n
			module load khmer
			normalize-by-median.py -o $baseDir/NormalizedReads/${read}.fq -x 1000000000 -C ${C} --n_tables 99 --force --ksize 32 $baseDir/FilteredReads/${read}.fq";
			close $OUT;
			system("sbatch $baseDir/scripts/normalize_${read}.slurm");
		}
	}

	if ($mergenorm) {
		print "Merging all reads...\n";
		open ($OUT, ">$baseDir/NormalizedReads/${single}_preclean.fq");
		foreach my $read (@reads) {
			open ($IN, "$baseDir/NormalizedReads/${read}.fq");
			while (my $line = readline($IN)) {
				print $OUT $line;
			}
		}
		close $OUT;
		close $IN;

		print "Tidying headers of merged reads file...\n"; #redo headers so that they are all >${single}.${i}/1
		open ($OUT, ">$baseDir/NormalizedReads/${single}.fq") or die "Can't open $baseDir/NormalizedReads/${single}_preclean.fq\n";
		open ($IN, "$baseDir/NormalizedReads/${single}_preclean.fq") or die "Can't open $baseDir/NormalizedReads/${single}.fq\n";
		$i = 1;
		while (my $line = readline($IN)) { 
			$line =~ s/^@.*$/\@${single}.${i}\/1/ if $i % 4 == 1;
			print $OUT $line;
			$i++;
		}
		close $OUT;
		close $IN;
		`rm $baseDir/NormalizedReads/${single}_preclean.fq`;

		print "Dispatching normalization job for merged single reads...\n";
			open ($OUT, ">$baseDir/scripts/mergenorm_${project}.slurm");
			print $OUT "#!/bin/sh\n#SBATCH --time=168:00:00\n#SBATCH --mem=${M}GB\n#SBATCH --job-name=${project}_${project}_normalize\n#SBATCH --error=$baseDir/scripts/LOG/normalize_${project}.err\n#SBATCH --output=$baseDir/scripts/LOG/normalize_${project}.out\n
			module load khmer
			normalize-by-median.py -o $baseDir/NormalizedMergedReads/${single}.fq -x 1000000000 -C ${C} --n_tables 99 --force --ksize 32 $baseDir/NormalizedReads/${single}.fq";
			close $OUT;
			system("sbatch $baseDir/scripts/mergenorm_${project}.slurm");
	}

	if ($assemble) {
		print "Calculating average read length...\n";
			my @seqs = `sed -n '2~4p' $baseDir/NormalizedMergedReads/${single}.fq`;
			foreach my $line (@seqs) {
				my $len=length($line);
				push @readlength, $len;
			}
	   		$nReads = @readlength;
	   		# print "$nReads";
			$nBases = eval join '+', @readlength;
			# print "$nBases";
	    	$avgLength = $nBases / $nReads;
	    	# $stdevLength += 
		    print "number of reads: $nReads\n";
		    print "number of bases: $nBases\n";
		    print "avg read length: $avgLength\n";

			if ($avgLength < 127) {
				$kMax = floor($avgLength);
			} else {
				$kMax = 127;
			}

			print "kMax: ${kMax}\n";

		if ($idba) {
			print "Dispatching idba run...\n";
				open ($OUT, ">$baseDir/scripts/idba_${single}.slurm") or die "Can't open $baseDir/scripts/idba_${single}.slurm";
					if ($avgLength <= 128) {
						print $OUT "#!/bin/sh\n#SBATCH --mem=${M}GB\n#SBATCH --time=168:00:00\n#SBATCH --job-name=${single}_idba\n#SBATCH --error=$baseDir/scripts/LOG/idba_${single}.err\n#SBATCH --output=$baseDir/scripts/LOG/idba_${single}.out\n
						module load idba khmer
						fastq-to-fasta.py -o $baseDir/NormalizedMergedReads/${single}.fa $baseDir/NormalizedMergedReads/${single}.fq
						idba_tran -o $baseDir/Assemblies/idba/ -r $baseDir/NormalizedMergedReads/${single}.fa --mink $kMin --maxk $kMax --step 7 --num_threads 8 --max_isoforms 50";
					}
					if ($avgLength > 128) {
						print $OUT "#!/bin/sh\n#SBATCH --mem=${M}GB\n#SBATCH --time=168:00:00\n#SBATCH --job-name=${single}_idba\n#SBATCH --error=$baseDir/scripts/LOG/idba_${single}.err\n#SBATCH --output=$baseDir/scripts/LOG/idba_${single}.out\n
						module load idba khmer
						fastq-to-fasta.py -o $baseDir/NormalizedMergedReads/${single}.fa $baseDir/NormalizedMergedReads/${single}.fq
						idba_tran -o $baseDir/Assemblies/idba/ -l $baseDir/NormalizedMergedReads/${single}.fa --mink $kMin --maxk 124 --step 7 --num_threads 8 --max_isoforms 50";
					}
				close $OUT;
				system ("sbatch $baseDir/scripts/idba_${single}.slurm");
		}

		if ($soap) {
			print "Creating soap configuration file...\n";
				open ($OUT, ">$baseDir/Assemblies/soap/SOAP_Conf.txt") or die "Can't open $baseDir/Assemblies/soap/SOAP_Conf.txt";
					print $OUT "\nmax_rd_len=600\n[LIB]\nrd_len_cutof=600\nreverse_seq=0\nasm_flags=3\nmap_len=32\nq=$baseDir/NormalizedMergedReads/${single}.fq";
				close $OUT;

			print "Dispatching soap runs...\n";
				for (my $k = $kMin; $k <= $kMax; $k += 4) {
					open ($OUT, ">$baseDir/scripts/soap_${k}_${single}.slurm") or die "Can't open $baseDir/scripts/soap_${k}_${single}.slurm";
						print $OUT "#!/bin/sh\n#SBATCH --mem=${M}GB\n#SBATCH --time=168:00:00\n#SBATCH --job-name=${single}_soap_${k}\n#SBATCH --error=$baseDir/scripts/LOG/soap_${k}_${single}.err\n#SBATCH --output=$baseDir/scripts/LOG/soap_${k}_${single}.out\n
						module load soapdenovo-trans
						SOAPdenovo-Trans-127mer all -s $baseDir/Assemblies/soap/SOAP_Conf.txt -K $k -o $baseDir/Assemblies/soap/$k -F -p 4";
					close $OUT;
					system ("sbatch $baseDir/scripts/soap_${k}_${single}.slurm");
				}
		}

		if ($spades) {
			print "Dispatching spades runs...\n";
				for (my $k = $kMin; $k <= $kMax; $k += 4) {
					open ($OUT, ">$baseDir/scripts/spades_${k}_${single}.slurm") or die "Can't open $baseDir/scripts/spades_${k}_${single}.slurm";
						print $OUT "#!/bin/sh\n#SBATCH --mem=${M}GB\n#SBATCH --time=168:00:00\n#SBATCH --job-name=${single}_spades_${k}\n#SBATCH --error=$baseDir/scripts/LOG/spades_${k}_${single}.err\n#SBATCH --output=$baseDir/scripts/LOG/spades_${k}_${single}.out\n
						module load spades
						rnaspades.py -k ${k} -o $baseDir/Assemblies/spades/${k} -s $baseDir/NormalizedMergedReads/${single}.fq -m ${M} -t 8";
					close $OUT;
					system ("sbatch $baseDir/scripts/spades_${k}_${single}.slurm");
				}
		}

		if ($trinity) {
			print "Dispatching trinity runs...\n";
				for (my $k = $kMin; $k <= 32; $k += 4) {
					open ($OUT, ">$baseDir/scripts/trinity_${k}_${single}.slurm") or die "Can't open $baseDir/scripts/trinity_${k}_${single}.slurm";
						print $OUT "#!/bin/sh\n#SBATCH --mem=${M}GB\n#SBATCH --time=168:00:00\n#SBATCH --job-name=${single}_trinity_${k}\n#SBATCH --error=$baseDir/scripts/LOG/trinity_${k}_${single}.err\n#SBATCH --output=$baseDir/scripts/LOG/trinity_${k}_${single}.out\n
						module load compiler/gcc bowtie samtools trinity
						Trinity --KMER_SIZE ${k} --no_normalize_reads --no_version_check --seqType fq --max_memory ${M}G --output $baseDir/Assemblies/trinity/${k}/trinity/ --full_cleanup --single $baseDir/NormalizedMergedReads/${single}.fq";
					close $OUT;
					system ("sbatch $baseDir/scripts/trinity_${k}_${single}.slurm");
				}
		}
	}

	if ($squant and $T) {
		print "Setting up salmon directories and dispatching salmon jobs...\n";
			open ($OUT, ">salmon_${T}.slurm") or die "Can't open salmon_${T}.slurm.\n";
			if ($single) {
				print $OUT "#!/bin/sh\n#SBATCH --time=12:00:00\n#SBATCH --mem=${M}GB\n#SBATCH --job-name=${project}_SAMPLE_salmon_${T}\n#SBATCH --error=$baseDir/scripts/LOG/salmon_${T}_SAMPLE.err\n#SBATCH --output=$baseDir/scripts/LOG/salmon_${T}_SAMPLE.out\n
				module load salmon
				salmon quant -i $baseDir/Mining/salmon/${T}_index -l A -r $baseDir/RawReads/SAMPLE.fq -o $baseDir/Mining/salmon/${T}/SAMPLE";
				close $OUT;
			}
			foreach my $x (@reads) {
				if ( -d "$baseDir/Mining/salmon/${T}/$x" ) {
			    print "Overwriting $baseDir/Mining/salmon${T}/$x\n";
			    	`rm -r $baseDir/Mining/salmon/${T}/$x`;
			    	mkdir "$baseDir/Mining/salmon/${T}/$x";
				} else {
				print "Creating directory $baseDir/Mining/salmon/${T}/$x...\n";
				mkdir "$baseDir/Mining/salmon/${T}/$x";
				}
				system ("sed \"s/SAMPLE/$x/g\" \"salmon_${T}.slurm\" > \"salmon_${T}_auto.slurm\"");
				system ("sbatch \"salmon_${T}_auto.slurm\"");
			}
	}
	print "Done.\n";
}

if ($paired) {

	print "Setting up reads...\n";
	if (-e "$baseDir/scripts/$paired") { 
		print "	Found $baseDir/scripts/$paired\n"; 
		} else { 
		die "	Can't find $baseDir/scripts/$paired\n";
	}

	print "	Importing read list...\n";
		my @reads = `sed -n '0~2!p' $baseDir/scripts/$paired`;
		chomp @reads;
		print "	Reads are: @reads\n";

	if ($get) {
		print "Submitting jobs to get SRA files...\n";
			foreach my $read (@reads) {
			open ($OUT, ">fastqdump_${read}.slurm") or die "Can't open fastqdump_${read}.slurm";
			print $OUT "#!/bin/sh\n#SBATCH --time=12:00:00\n#SBATCH --mem=16GB\n#SBATCH --job-name=${project}_${read}_fqd\n#SBATCH --error=$baseDir/scripts/LOG/fqd_${read}.err\n#SBATCH --output=$baseDir/scripts/LOG/fqd_${read}.out
			cd $baseDir/RawReads
			module load SRAtoolkit
			prefetch ${read}
			fastq-dump --defline-seq \'@\$sn[_\$rn]/\$ri\' --split-files ${read}
			mv $baseDir/RawReads/${read}_1.fastq $baseDir/RawReads/${read}_1.fq
			mv $baseDir/RawReads/${read}_2.fastq $baseDir/RawReads/${read}_2.fq";
			close $OUT;
			system("sbatch fastqdump_${read}.slurm");
			}
	}

	if ($filter) {
		print "Dispatching filtering job(s) for paired reads at q = $q...\n";
			foreach my $read (@reads) {
				open ($OUT, ">$baseDir/scripts/filter_${read}_q${q}.slurm");
				print $OUT "#!/bin/sh\n#SBATCH --time=12:00:00\n#SBATCH --mem=${M}GB\n#SBATCH --job-name=${project}_${read}_filter_q$q\n#SBATCH --error=$baseDir/scripts/LOG/filter_${read}_q$q.err\n#SBATCH --output=$baseDir/scripts/LOG/filter_${read}_q$q.out\n
				module load erne
				erne-filter --query1 $baseDir/RawReads/${read}_1.fq --query2 $baseDir/RawReads/${read}_2.fq --output-prefix $baseDir/FilteredReads/${read} --ultra-sensitive --force-standard
				mv $baseDir/FilteredReads/${read}_1.fastq $baseDir/FilteredReads/${read}_1.fq
				mv $baseDir/FilteredReads/${read}_2.fastq $baseDir/FilteredReads/${read}_2.fq
				rm $baseDir/FilteredReads/${read}_unpaired.fastq";
				close $OUT;
				system("sbatch $baseDir/scripts/filter_${read}_q${q}.slurm");
			}
	}

	if ($normalize) {
		print "Dispatching normalization job(s) for paired reads...\n";
		foreach my $read (@reads) {
			open ($OUT, ">$baseDir/scripts/normalize_${read}.slurm");
			print $OUT "#!/bin/sh\n#SBATCH --time=168:00:00\n#SBATCH --mem=${M}GB\n#SBATCH --job-name=${project}_${read}_normalize\n#SBATCH --error=$baseDir/scripts/LOG/normalize_${read}.err\n#SBATCH --output=$baseDir/scripts/LOG/normalize_${read}.out\n
			module load khmer
			interleave-reads.py  -o $baseDir/FilteredReads/${read}.fq.I --force $baseDir/FilteredReads/${read}_1.fq $baseDir/FilteredReads/${read}_2.fq
			normalize-by-median.py -o $baseDir/NormalizedReads/${read}.fq.I -x 1000000000 -C 50 --n_tables 99 --force --ksize 32 $baseDir/FilteredReads/${read}.fq.I
			extract-paired-reads.py -d $baseDir/NormalizedReads $baseDir/NormalizedReads/${read}.fq.I
			split-paired-reads.py -d $baseDir/NormalizedReads --force $baseDir/NormalizedReads/${read}.fq.I.pe
			rm $baseDir/NormalizedReads/${read}.fq.I.se
			mv $baseDir/NormalizedReads/${read}.fq.I.pe $baseDir/NormalizedReads/${read}.fq.I
			mv $baseDir/NormalizedReads/${read}.fq.I.pe.1 $baseDir/NormalizedReads/${read}_1.fq
			mv $baseDir/NormalizedReads/${read}.fq.I.pe.2 $baseDir/NormalizedReads/${read}_2.fq";
			close $OUT;
			system("sbatch $baseDir/scripts/normalize_${read}.slurm");
		}
	}

	if ($mergenorm) {
		print "Merging reads in ${paired}...\n";
		open ($OUT, ">$baseDir/NormalizedReads/${paired}_1.fq");
		foreach my $read (@reads) {
			open ($IN, "$baseDir/NormalizedReads/${read}_1.fq");
			while (my $line = readline($IN)) {
				print $OUT $line;
			}
		}
		close $OUT;
		open ($OUT, ">$baseDir/NormalizedReads/${paired}_2.fq");
		foreach my $read (@reads) {
			open ($IN, "$baseDir/NormalizedReads/${read}_2.fq");
			while (my $line = readline($IN)) {
				print $OUT $line;
			}
		}
		close $OUT;
		print "Dispatching normalization job for merged single reads...\n";
			open ($OUT, ">$baseDir/scripts/mergenorm_${paired}.slurm");
			print $OUT "#!/bin/sh\n#SBATCH --time=168:00:00\n#SBATCH --mem=${M}GB\n#SBATCH --job-name=${paired}_normalize\n#SBATCH --error=$baseDir/scripts/LOG/normalize_${paired}.err\n#SBATCH --output=$baseDir/scripts/LOG/normalize_${paired}.out\n
			module load khmer
			interleave-reads.py -o $baseDir/NormalizedReads/${paired}.fq.I --force $baseDir/NormalizedReads/${paired}_1.fq $baseDir/NormalizedReads/${paired}_2.fq
			normalize-by-median.py -o $baseDir/NormalizedMergedReads/${paired}.fq.I -x 1000000000 -C ${C} --n_tables 99 --force --ksize 32 $baseDir/NormalizedReads/${paired}.fq.I
			extract-paired-reads.py -d $baseDir/NormalizedMergedReads $baseDir/NormalizedMergedReads/${read}.fq.I
			split-paired-reads.py -d $baseDir/NormalizedMergedReads --force $baseDir/NormalizedMergedReads/${read}.fq.I.pe
			rm $baseDir/NormalizedMergedReads/${read}.fq.I.se
			mv $baseDir/NormalizedMergedReads/${read}.fq.I.pe $baseDir/NormalizedMergedReads/${read}.fq.I
			mv $baseDir/NormalizedMergedReads/${read}.fq.I.pe.1 $baseDir/NormalizedMergedReads/${read}_1.fq
			mv $baseDir/NormalizedMergedReads/${read}.fq.I.pe.2 $baseDir/NormalizedMergedReads/${read}_2.fq";
			close $OUT;
			system("sbatch $baseDir/scripts/mergenorm_${paired}.slurm");
	}

	if ($skipmergenorm) {
		print "Copying reads to NormalizedMergedReads...\n";
		`cp $baseDir/NormalizedReads/* $baseDir/NormalizedMergedReads`;
	}

	if ($assemble) {
		print "Calculating average read length...\n";
			my @seqs = `sed -n '2~4p' $baseDir/NormalizedReads/normalized_${project}.fq`;
			foreach my $line (@seqs) {
				my $len=length($line);
				push @readlength, $len;
			}
	   		$nReads = @readlength;
	   		# print "$nReads";
			$nBases = eval join '+', @readlength;
			# print "$nBases";
	    	$avgLength = $nBases / $nReads;
	    	# $stdevLength += 
		    print "number of reads: $nReads\n";
		    print "number of bases: $nBases\n";
		    print "avg read length: $avgLength\n";

			if ($avgLength < 127) {
				$kMax = floor($avgLength);
			} else {
				$kMax = 127;
			}

			print "kMax: ${kMax}\n";

		if ($idba) {
			print "Dispatching idba run...\n";
			open ($OUT, ">$baseDir/scripts/idba_${paired}.slurm") or die "Can't open $baseDir/scripts/idba_${paired}.slurm";
				if ($avgLength <= 128) {
					print $OUT "#!/bin/sh\n#SBATCH --mem=${M}GB\n#SBATCH --time=168:00:00\n#SBATCH --job-name=${paired}_idba\n#SBATCH --error=$baseDir/scripts/LOG/idba_${paired}.err\n#SBATCH --output=$baseDir/scripts/LOG/idba_${paired}.out\n
					module load idba khmer
					fastq-to-fasta.py -o $baseDir/NormalizedMergedReads/${paired}.fa $baseDir/NormalizedMergedReads/${paired}.fq
					idba_tran -o $baseDir/Assemblies/idba/ -r $baseDir/NormalizedMergedReads/${paired}.fa --mink $kMin --maxk $kMax --step 7 --num_threads 8 --max_isoforms 50";
				}
				if ($avgLength > 128) {
					print $OUT "#!/bin/sh\n#SBATCH --mem=${M}GB\n#SBATCH --time=168:00:00\n#SBATCH --job-name=${paired}_idba\n#SBATCH --error=$baseDir/scripts/LOG/idba_${paired}.err\n#SBATCH --output=$baseDir/scripts/LOG/idba_${paired}.out\n
					module load idba khmer
					fastq-to-fasta.py -o $baseDir/NormalizedMergedReads/${paired}.fa $baseDir/NormalizedMergedReads/${paired}.fq
					idba_tran -o $baseDir/Assemblies/idba/ -l $baseDir/NormalizedMergedReads/${paired}.fa --mink $kMin --maxk $kMax --step 7 --num_threads 8 --max_isoforms 50";
				}
			close $OUT;
			system ("sbatch $baseDir/scripts/idba_${paired}.slurm");
		}

		if ($soap) {
			print "Creating soap configuration file..."; ###SHOULD REVERSE_SEQ HERE BE == 1?
				open ($OUT, ">$baseDir/Assemblies/soap/SOAP_Conf.txt") or die "Can't open $baseDir/Assemblies/soap/SOAP_Conf.txt";
					print $OUT "\nmax_rd_len=600\n[LIB]\nrd_len_cutof=600\nreverse_seq=0\nasm_flags=3\nmap_len=32\nq=$baseDir/NormalizedReads/normalized_${project}.fq";
				close $OUT;

			print "Dispatching soap runs...\n";
				for (my $k = $kMin; $k <= $kMax; $k += 4) {
					open ($OUT, ">$baseDir/scripts/soap_${k}_${project}.slurm") or die "Can't open $baseDir/scripts/soap_${k}_${project}.slurm";
						print $OUT "#!/bin/sh\n#SBATCH --mem=${M}GB\n#SBATCH --time=168:00:00\n#SBATCH --job-name=${project}_soap_${k}\n#SBATCH --error=$baseDir/scripts/LOG/soap_${k}_${project}.err\n#SBATCH --output=$baseDir/scripts/LOG/soap_${k}_${project}.out\n
						module load soapdenovo-trans
						SOAPdenovo-Trans-127mer all -s $baseDir/Assemblies/soap/SOAP_Conf.txt -K $k -o $baseDir/Assemblies/soap/$k -F -p 4";
					close $OUT;
					system ("sbatch $baseDir/scripts/soap_${k}_${project}.slurm");
				}
		}

		if ($spades) { 
			print "Dispatching spades runs...\n";
				for (my $k = $kMin; $k <= $kMax; $k += 4) {
					open ($OUT, ">$baseDir/scripts/spades_${k}_${project}.slurm") or die "Can't open $baseDir/scripts/spades_${k}_${project}.slurm";
						print $OUT "#!/bin/sh\n#SBATCH --mem=${M}GB\n#SBATCH --time=168:00:00\n#SBATCH --job-name=${project}_spades_${k}\n#SBATCH --error=$baseDir/scripts/LOG/spades_${k}_${project}.err\n#SBATCH --output=$baseDir/scripts/LOG/spades_${k}_${project}.out\n
						module load spades
						rnaspades.py -k ${k} -o $baseDir/Assemblies/spades/${k} -s $baseDir/NormalizedReads/normalized_${project}.fq -m ${M} -t 8";
					close $OUT;
					system ("sbatch $baseDir/scripts/spades_${k}_${project}.slurm");		
				}
		}

		if ($trinity) {
			print "Dispatching trinity runs...\n"; #trinity has max kmer size of 32
				for (my $k = $kMin; $k <= 32; $k += 4) {
					open ($OUT, ">$baseDir/scripts/trinity_${k}_${project}.slurm") or die "Can't open $baseDir/scripts/trinity_${k}_${project}.slurm";
						print $OUT "#!/bin/sh\n#SBATCH --mem=${M}GB\n#SBATCH --time=168:00:00\n#SBATCH --job-name=${project}_trinity_${k}\n#SBATCH --error=$baseDir/scripts/LOG/trinity_${k}_${project}.err\n#SBATCH --output=$baseDir/scripts/LOG/trinity_${k}_${project}.out\n
						module load compiler/gcc bowtie samtools trinity
						Trinity --KMER_SIZE ${k} --no_normalize_reads --no_version_check --seqType fq --max_memory ${M}G --output $baseDir/Assemblies/trinity/${k}/trinity/ --full_cleanup --single $baseDir/NormalizedReads/normalized_${project}.fq";
					close $OUT;
					system ("sbatch $baseDir/scripts/trinity_${k}_${project}.slurm");		
				}
		}
	}

	if ($squant and $T) {
		print "Setting up salmon directories and dispatching salmon jobs...\n";
			open ($OUT, ">salmon_${T}.slurm") or die "Can't open salmon_${T}.slurm.\n";
			print $OUT "#!/bin/sh\n#SBATCH --time=12:00:00\n#SBATCH --mem=${M}GB\n#SBATCH --job-name=${project}_SAMPLE_salmon_${T}\n#SBATCH --error=$baseDir/scripts/LOG/salmon_${T}_SAMPLE.err\n#SBATCH --output=$baseDir/scripts/LOG/salmon_${T}_SAMPLE.out\n
				module load salmon
				salmon quant -i $baseDir/Mining/salmon/${T}_index -l A -1 $baseDir/RawReads/SAMPLE_1.fq -2 $baseDir/RawReads/SAMPLE_2.fq -o $baseDir/Mining/salmon/${T}/SAMPLE";
			close $OUT;
			foreach my $x (@reads) {
				if ( -d "$baseDir/Mining/salmon/${T}/$x" ) {
			    print "Overwriting $baseDir/Mining/salmon/${T}/$x\n";
			    	`rm -r $baseDir/Mining/salmon/${T}/$x`;
			    	mkdir "$baseDir/Mining/salmon/${T}/$x";
				} else {
				print "Creating directory $baseDir/Mining/salmon/${T}/$x...\n";
				mkdir "$baseDir/Mining/salmon/${T}/$x";
			}
				system ("sed \"s/SAMPLE/$x/g\" \"salmon_${T}.slurm\" > \"salmon_${T}_auto.slurm\"");
				system ("sbatch \"salmon_${T}_auto.slurm\"");
			}
		print "Done.\n";
	}
}

if ($kmerge) {

	print "Deleting previous kmerges...\n";
		foreach my $assembler (@assemblers) {
			if (-e "$baseDir/Assemblies/${assembler}/all_${assembler}_contigs.fa") {
				unlink "$baseDir/Assemblies/${assembler}/all_${assembler}_contigs.fa";
			}
		}

	print "Creating one file for each assembler with all the contigs it found using various kmers...\n";
		`awk 'FNR==1{print ""}1' $baseDir/Assemblies/idba/transcript-*.fa > $baseDir/Assemblies/idba/all_idba_contigs.fa`;
		`awk 'FNR==1{print ""}1' $baseDir/Assemblies/soap/*.scafSeq > $baseDir/Assemblies/soap/all_soap_contigs.fa`;
		`awk 'FNR==1{print ""}1' $baseDir/Assemblies/spades/*/transcripts.fasta > $baseDir/Assemblies/spades/all_spades_contigs.fa`;
		`awk 'FNR==1{print ""}1' $baseDir/Assemblies/trinity/*/trinity.Trinity.fasta > $baseDir/Assemblies/trinity/all_trinity_contigs.fa`;
		# system ("cat $baseDir/Assemblies/idba/transcript-*.fa > $baseDir/Assemblies/idba/all_idba_contigs.fa");
		# system ("cat $baseDir/Assemblies/soap/*.scafSeq > $baseDir/Assemblies/soap/all_soap_contigs.fa");
		# system ("cat $baseDir/Assemblies/spades/*/transcripts.fasta > $baseDir/Assemblies/spades/all_spades_contigs.fa");
		# system ("cat $baseDir/Assemblies/trinity/*/trinity.Trinity.fasta > $baseDir/Assemblies/trinity/all_trinity_contigs.fa");
	print "Done.\n";
}

if ($translatedenovos) {

	print "Deleting previous translations...\n";
	foreach my $assembler (@assemblers) {
		if (-e "$baseDir/Assemblies/${assembler}/all_${assembler}_contigs.fa.faa") {
			unlink "$baseDir/Assemblies/${assembler}/all_${assembler}_contigs.fa.faa";
		}
	}

	foreach my $assembler (@assemblers) {
	print "Submitting ${assembler} translation job...\n";
		open (my $OUT, ">$baseDir/scripts/genemarks_${assembler}.slurm") or die "Could not open $baseDir/scripts/genemarks_${assembler}.slurm";
		print $OUT "#!/bin/sh \n#SBATCH --mem=${M}GB \n#SBATCH --time=168:00:00 \n#SBATCH --job-name=${project}_gmsn_${assembler} \n#SBATCH --error=$baseDir/scripts/LOG/gmsn_${assembler}.err \n#SBATCH --output=$baseDir/scripts/LOG/gmsn_${assembler}.out
		cd $baseDir/Assemblies/${assembler}
		echo 'Loading modules...'
			ml genemarks/4.3
		echo 'Translating...'
			gmsn.pl --faa --euk all_${assembler}_contigs.fa
		echo 'Done.'\n";
		close $OUT;
		system("sbatch genemarks_${assembler}.slurm");
	}
	print "Done.\n";
}

if ($nonredund) {

	foreach my $assembler (@assemblers) {

		print "Deleting previous nonredund AA(RNA) selections...\n";
			if (-e "$baseDir/Assemblies/${assembler}/nonredund_${assembler}_contigs.fa") {
				unlink "$baseDir/Assemblies/${assembler}/nonredund_${assembler}_contigs.fa";
			}

		$file = "$baseDir/Assemblies/${assembler}/all_${assembler}_contigs.fa";
		print "Tidying and linearizing ${file}...\n";
			open ($IN, "${file}") or die "Can't open ${file}\n";
			open ($OUT, ">${file}_l") or die "Can't open ${file}_l\n";
			while (my $line = readline($IN)) {
				$line =~ s/ /_/g;
				$line =~ s/^>(.*)$/#>$1#/;
				$line =~ s/\n//g;
				$line =~ s/#/\n/g;
				print $OUT $line;
			}
			close $IN;
			close $OUT;
			system ("sed -i '1d' ${file}_l");
			unlink "${file}";
			system ("mv ${file}_l ${file}");

		$file = "$baseDir/Assemblies/${assembler}/all_${assembler}_contigs.fa.faa";
		print "Tidying and linearizing ${file}...\n";
			open ($IN, "${file}") or die "Can't open ${file}\n";
			open ($OUT, ">${file}_l") or die "Can't open ${file}_l\n";
			while (my $line = readline($IN)) {
				$line =~ s/ /_/g;
				$line =~ s/^(>.*)>/>/g;
				$line =~ s/^>(.*)$/#>$1#/;
				$line =~ s/\n//g;
				$line =~ s/#/\n/g;
				print $OUT $line;
			}
			close $IN;
			close $OUT;
			system ("sed -i '1d' ${file}_l");
			unlink "${file}";
			system ("mv ${file}_l ${file}");

		print "Preparing array of hashes containing all_${assembler}_contigs and their translations...\n";
		$file = "$baseDir/Assemblies/${assembler}/all_${assembler}_contigs.fa";
			my %rna;
			my @rna = `cat ${file}`;
			chomp @rna;
			%rna = @rna;

			my %aa;
			my @aa = `cat ${file}.faa`;
			chomp @aa;
			%aa = @aa;

			my @contig = `sed -n '0~2!p' ${file}.faa`;
			chomp @contig;

			my @results;
			foreach (@contig) {
				push @results, { rna => "$rna{\"$_\"}", aa => "$aa{\"$_\"}" };
			}

		print "Sorting array of all_${assembler}_contigs based on RNA length...\n";
			my @sorted_results;
			@sorted_results = sort { length $a->{rna} cmp length $b->{rna} } @results;

		print "Selecting and writing shortest nonredund ${assembler} proteins and their transcripts...\n";
			my @uniq_results = grep { !$seen{$_->{aa}}++ } @sorted_results;           
			# print Dumper(\@uniq_idba);
			open ($OUT, ">$baseDir/Assemblies/${assembler}/nonredund_${assembler}_contigs.fa") or die "Can't open $baseDir/Assemblies/${assembler}/nonredund_${assembler}_contigs.fa\n";
			$i = 1;
			foreach  (@uniq_results) {
				print $OUT ">${project}_${assembler}_${i}-\n";
				print $OUT "$_->{rna}\n";
				$i = $i++;
			}
			close $OUT;
			open ($OUT, ">$baseDir/Assemblies/${assembler}/nonredund_${assembler}_contigs.fa.faa") or die "Can't open $baseDir/Assemblies/${assembler}/nonredund_${assembler}_contigs.fa.faa\n";
			$i = 1;
			foreach  (@uniq_results) {
				print $OUT ">${project}_${assembler}_${i}-\n";
				print $OUT "$_->{aa}\n";
				$i = $i++;
			}
			close $OUT;
	}
	print "Done!\n";
}

if ($consensus) {

	print "Preparing all_${project}_contigs.fa...\n";
		open ($OUT, ">$baseDir/Assemblies/all_${project}_contigs.fa") or die "Can't open $baseDir/Assemblies/all_${project}_contigs.fa";
		foreach my $assembler (@assemblers) {
			open ($IN, "$baseDir/Assemblies/${assembler}/nonredund_${assembler}_contigs.fa") or die "Can't open $baseDir/Assemblies/${assembler}/nonredund_${assembler}_contigs.fa";;
				while (my $line = readline($IN)) {
					print $OUT $line;
				}
			close $IN;
		}
		close $OUT;

	print "Tidying and linearizing all_${project}_contigs.fa\n";
			$file = "$baseDir/Assemblies/all_${project}_contigs.fa";
			print "Tidying and linearizing ${file}...\n";
				open ($IN, "${file}") or die "Can't open ${file}\n";
				open ($OUT, ">${file}_l") or die "Can't open ${file}_l\n";
				while (my $line = readline($IN)) {
					$line =~ s/ /_/g;
					$line =~ s/^(>.*)>/>/g;
					$line =~ s/^>(.*)$/#>$1#/;
					$line =~ s/\n//g;
					$line =~ s/#/\n/g;
					print $OUT $line;
				}
				close $IN;
				close $OUT;
				system ("sed -i '1d' ${file}_l");
				unlink "${file}";
				system ("mv ${file}_l ${file}");

	print "Preparing all_${project}_contigs.fa.faa...\n";
		open ($OUT, ">$baseDir/Assemblies/all_${project}_contigs.fa.faa") or die "Can't open $baseDir/Assemblies/all_${project}_contigs.fa.faa";
		foreach my $assembler (@assemblers) {
			open ($IN, "$baseDir/Assemblies/${assembler}/nonredund_${assembler}_contigs.fa.faa") or die "Can't open $baseDir/Assemblies/${assembler}/nonredund_${assembler}_contigs.fa.faa";
				while (my $line = readline($IN)) {
					print $OUT $line;
				}
			close $IN;
		}
		close $OUT;

	print "Tidying and linearizing all_${project}_contigs.fa.faa\n";
			$file = "$baseDir/Assemblies/all_${project}_contigs.fa.faa";
			print "Tidying and linearizing ${file}...\n";
				open ($IN, "${file}") or die "Can't open ${file}\n";
				open ($OUT, ">${file}_l") or die "Can't open ${file}_l\n";
				while (my $line = readline($IN)) {
					$line =~ s/ /_/g;
					$line =~ s/^(>.*)>/>/g;
					$line =~ s/^>(.*)$/#>$1#/;
					$line =~ s/\n//g;
					$line =~ s/#/\n/g;
					print $OUT $line;
				}
				close $IN;
				close $OUT;
				system ("sed -i '1d' ${file}_l");
				unlink "${file}";
				system ("mv ${file}_l ${file}");

	print "Preparing array of hashes containing all_${project}_contigs and their translations...\n";
		$file = "$baseDir/Assemblies/all_${project}_contigs.fa";
		my %rna;
		my @rna = `cat "${file}"`;
		chomp @rna;
		%rna = @rna;

		my %aa;
		my @aa = `cat "${file}.faa"`;
		chomp @aa;
		%aa = @aa;

		my @contig = `sed -n '0~2!p' ${file}.faa`;
		chomp @contig;

		my @project_results;
		foreach (@contig) {
			push @project_results, { rna => "$rna{\"$_\"}", aa => "$aa{\"$_\"}", contig => "$_"};
		}

		print "Sorting array of all_${project}_contigs based on RNA length...\n";
			my @sorted_project_results;
			@sorted_project_results = sort { length $a->{rna} cmp length $b->{rna} } @project_results;

		print "Selecting AA sequences found by at least two assemblers...\n";
			my @dupl_results = grep { $seen{$_->{aa}}++ } @sorted_project_results;
			# print Dumper(\@uniq_results);

		print "Sorting array of all_${project}_contigs based on AA length...\n";
			my @resorted_results;
			@resorted_results = sort { length $a->{rna} cmp length $b->{rna} } @dupl_results;

		print "Selecting shortest RNA/nonredund AA sequences found by at least two assemblers...\n";
			my @nonredund_dupl_results = grep { !$seen{$_->{aa}}++ } @resorted_results;
			# print Dumper(\@uniq_results);

			open ($OUT, ">$baseDir/Assemblies/consensus_${project}_contigs.fa") or die "Can't open $baseDir/Assemblies/consensus_${project}_contigs.fa\n";
			$i = 1;
			my $orig;
			foreach  (@nonredund_dupl_results) {
				my $orig = substr $_->{contig}, 1;
				print $OUT ">${project}_consensus_${i}-(${orig})\n";
				print $OUT "$_->{rna}\n";
				$i = $i++;
			}
			close $OUT;
			open ($OUT, ">$baseDir/Assemblies/consensus_${project}_contigs.fa.faa") or die "Can't open $baseDir/Assemblies/consensus_${project}_contigs.fa.faa\n";
			$i = 1;
			foreach  (@nonredund_dupl_results) {
				my $orig = substr $_->{contig}, 1;
				print $OUT ">${project}_consensus_${i}-$_-($orig})\n";
				print $OUT "$_->{aa}\n";
				$i = $i++;
			}
			close $OUT;
			print "Done.\n";
}

if ($sindex and $T) {
	
	if (-d "$baseDir/Mining/salmon/${T}_index") {
		print "Removing previous index...\n";
		`rm -r $baseDir/Mining/salmon/${T}_index`
	}
	print "Building salmon database for $baseDir/Assemblies/${T}...\n";
		`salmon index -t $baseDir/Assemblies/${T}.fa -i $baseDir/Mining/salmon/${T}_index -k 31`;
	print "Done."
}

if ($blast and $T and $Q and $anchor) {

	$file = "$baseDir/Assemblies/${T}.fa";
		print "Tidying and linearizing ${file}...\n";
		open ($IN, "${file}") or die "Can't open ${file}\n";
		open ($OUT, ">${file}_l") or die "Can't open ${file}_l\n";
		while (my $line = readline($IN)) {
			$line =~ s/ /_/g;
			$line =~ s/^>(.*)$/#>$1#/;
			$line =~ s/\n//g;
			$line =~ s/#/\n/g;
			print $OUT $line;
		}
		close $IN;
		close $OUT;
		system ("sed -i '1d' ${file}_l");
		unlink "${file}";
		system ("mv ${file}_l ${file}");

	print "Making BLASTable database from $baseDir/Assemblies/${T}.fa...\n";
		`makeblastdb -in $baseDir/Assemblies/${T}.fa -parse_seqids -dbtype nucl`;

	$file = "$baseDir/Mining/${Q}.faa";
		print "Tidying and linearizing ${file}...\n";
		open ($IN, "${file}") or die "Can't open ${file}\n";
		open ($OUT, ">${file}_l") or die "Can't open ${file}_l\n";
		while (my $line = readline($IN)) {
			$line =~ s/ /_/g;
			$line =~ s/^>(.*)$/#>$1#/;
			$line =~ s/\n//g;
			$line =~ s/#/\n/g;
			print $OUT $line;
		}
		close $IN;
		close $OUT;
		system ("sed -i '1d' ${file}_l");
		unlink "${file}";
		system ("mv ${file}_l ${file}");

	print "BLASTing $baseDir/Mining/${Q}.faa against $baseDir/Assemblies/${T}.fa...\n";
		system("tblastn -db $baseDir/Assemblies/${T}.fa -query $baseDir/Mining/${Q}.faa -out $baseDir/Mining/${Q}_${T}_hits.out -outfmt '6 sacc'");

	print "Removing duplicate hits...\n";
		`sort $baseDir/Mining/${Q}_${T}_hits.out | uniq > $baseDir/Mining/${Q}_${T}_hits_u.out`;
		`rm $baseDir/Mining/${Q}_${T}_hits.out`;
		`mv $baseDir/Mining/${Q}_${T}_hits_u.out $baseDir/Mining/${Q}_${T}_hits.out`;
		my @hits;
		open (my $HITS, "<", "$baseDir/Mining/${Q}_${T}_hits.out") or die "Failed to open file: $!\n";
		while (<$HITS>) {
			chomp;
			push @hits, $_;
		}
		close $HITS;
		# print "Hits:\n";
		# print join "\n", @hits;
		# print "\n";
		print "Number of nonredund hits:\n";
		print scalar @hits;
		print "\n";

	print "Extracting hit sequences...\n";
		`blastdbcmd -db $baseDir/Assemblies/${T}.fa -entry_batch $baseDir/Mining/${Q}_${T}_hits.out -outfmt %f -out $baseDir/Mining/${Q}_${T}_hits.fa`;

	# Printing results file...\n";
	# 	foreach my $i (@hits) {
	# 		my $header_line = `grep -n $i "$baseDir/Assemblies/${T}.fa" | cut -d : -f 1`;
	# 		chomp $header_line;
	# 		# print "header_line: ${header_line}\n";
	# 		my $seq_line = $header_line+1;
	# 		chomp $seq_line;
	# 		# print "seq_line: $seq_line\n";
	# 		my $header = `sed -n ${header_line}p $baseDir/Assemblies/${T}.fa`;
	# 		chomp $header;
	# 		push @hits_fasta, "$header";
	# 		# print "Header: $header\n";
	# 		my $seq = `sed -n ${seq_line}p $baseDir/Assemblies/${T}.fa`;
	# 		chomp $seq;
	# 		push @hits_fasta, "$seq";
	# 		# print "Seq: $seq\n";			
	# 	}
	# 	open ($OUT, ">$baseDir/Mining/${Q}_${T}_hits.fa") or die "Can't open $baseDir/Mining/${Q}_${T}_hits.fa\n";
	# 		print $OUT join("\n", @hits_fasta);
	# 	close($OUT);

	print "Translating $baseDir/Mining/${Q}_${T}_hits.fa to peptide then tidying and linearizing ...\n";
		chdir ("$baseDir/Mining");
		`gmsn.pl --faa --euk --clean $baseDir/Mining/${Q}_${T}_hits.fa`;
		chdir ("$baseDir/scripts");
		
		$file = "$baseDir/Mining/${Q}_${T}_hits.fa.faa";
		print "Tidying and linearizing ${file}...\n";
			open ($IN, "${file}") or die "Can't open ${file}\n";
			open ($OUT, ">${file}_l") or die "Can't open ${file}_l\n";
			while (my $line = readline($IN)) {
				# $line =~ s/.*>${project}/>${project}/;
				$line =~ s/.*>XM/>XM/;
				$line =~ s/ /_/g;
				$line =~ s/^>(.*)$/#>$1#/;
				$line =~ s/\n//g;
				$line =~ s/#/\n/g;
				print $OUT $line;
			}
			close $IN;
			close $OUT;
			system ("sed -i '1d' ${file}_l");
			unlink "${file}";
			system ("mv ${file}_l ${file}");	

	print "Merging hits and query...\n";
		`awk '{print}' "$baseDir/Mining/${Q}_${T}_hits.fa.faa" "$baseDir/Mining/${Q}.faa" > "$baseDir/Mining/${Q}_${T}_hqmerge.faa"`;
	

	print "Merging hits, query, and anchor...\n";
		`awk '{print}' "$baseDir/Mining/${Q}_${T}_hits.fa.faa" "$baseDir/Mining/${Q}.faa" > "$baseDir/Mining/${Q}_${T}_merge.faa"`;
		`awk '{print}' "$baseDir/Mining/${Q}_${T}_merge.faa" "$baseDir/Mining/${anchor}.faa" > "$baseDir/Mining/${Q}_${T}_hqmerge.faa"`;
		
	print "Linearizing and cleaning $baseDir/Mining/${Q}_${T}_hqmerge.faa...\n";
		$file = "$baseDir/Mining/${Q}_${T}_hqmerge.faa";
		print "Tidying and linearizing ${file}...\n";
		open ($IN, "${file}") or die "Can't open ${file}\n";
		open ($OUT, ">${file}_l") or die "Can't open ${file}_l\n";
		while (my $line = readline($IN)) {
			$line =~ s/ /_/g;
			$line =~ s/^>(.*)$/#>$1#/;
			$line =~ s/\n//g;
			$line =~ s/#/\n/g;
			print $OUT $line;
		}
		close $IN;
		close $OUT;
		system ("sed -i '1d' ${file}_l");
		unlink "${file}";
		system ("mv ${file}_l ${file}");

	print "Length filtering $baseDir/Mining/${Q}_${T}_hqmerge.faa to ${blast} +/- 100 AA...\n";
		$file = "$baseDir/Mining/${Q}_${T}_hqmerge.faa";
		my %fasta;
		my @fasta = `cat ${file}`;
		chomp @fasta;
		%fasta = @fasta;
		my %long = map { $_ => $fasta{$_} } grep { length($fasta{$_}) >= ${blast}-100 } keys(%fasta);
		my %short = map { $_ => $long{$_} } grep { length($long{$_}) <= ${blast}+100 } keys(%long);
		@short = %short;
		open ($OUT, ">${file}_f") or die "Can't open ${file}_f";
		foreach (@short) {
			print $OUT "$_\n"
		}
		close ($OUT);
		unlink ${file};
		system ("mv ${file}_f ${file}");

	print "Aligning $baseDir/Mining/${Q}_${T}_hqmerge.faa...\n";
		`clustalo -i $baseDir/Mining/${Q}_${T}_hqmerge.faa -o $baseDir/Mining/${Q}_${T}_hqmerge_align.faa --force`;

	print "Linearizing '$baseDir/Mining/${Q}_${T}_hqmerge_align.faa'...\n";
	$file = "$baseDir/Mining/${Q}_${T}_hqmerge_align.faa";
		print "Tidying and linearizing ${file}...\n";
			open ($IN, "${file}") or die "Can't open ${file}\n";
			open ($OUT, ">${file}_l") or die "Can't open ${file}_l\n";
			while (my $line = readline($IN)) {
				$line =~ s/ /_/g;
				$line =~ s/^>(.*)$/#>$1#/;
				$line =~ s/\n//g;
				$line =~ s/#/\n/g;
				print $OUT $line;
			}
			close $IN;
			close $OUT;
			system ("sed -i '1d' ${file}_l");
			unlink "${file}";
			system ("mv ${file}_l ${file}");
	
	print "Done!\n";
}

if ($collect) {

	print "Parsing subproject names...\n";

		my @subprojects = `cat $baseDir/scripts/${subprojects}`;
		chomp @subprojects;
		print join("\n", @subprojects);
		print "\n";

	print "Collecting all_contigs from each subproject...\n";

		if (-e "$baseDir/Assemblies/all_${project}_contigs.fa") {
			unlink "$baseDir/Assemblies/all_${project}_contigs.fa";
			print "Deleting previous all_${project}_contigs.fa ...\n";
		}
		open ($OUT, ">>$baseDir/Assemblies/all_${project}_contigs.fa");
		foreach my $subproject (@subprojects) {
			print "Copying all_${subproject}_contigs.fa to all_${project}_contigs.fa...\n";
			open ($IN, "$mainDir/${subproject}/Assemblies/all_${subproject}_contigs.fa") or die "Can't find $mainDir/${subproject}/Assemblies/all_${subproject}_contigs.fa";
			while (my $line = readline($IN)) {
				print $OUT $line;
				print $OUT "\n";
			}
			close $IN;
		}
		close $OUT;

		if (-e "$baseDir/Assemblies/all_${project}_contigs.fa.faa") {
			unlink "$baseDir/Assemblies/all_${project}_contigs.fa.faa";
			print "Deleting previous collections of all contigs.fa.faa ...\n";
		}
		open ($OUT, ">>$baseDir/Assemblies/all_${project}_contigs.fa.faa");
		foreach my $subproject (@subprojects) {
			print "Copying all_${subproject}_contigs.fa.faa to all_${project}_contigs.fa.faa...\n";
			open ($IN, "$mainDir/${subproject}/Assemblies/all_${subproject}_contigs.fa.faa") or die "Can't find $mainDir/${subproject}/Assemblies/all_${subproject}_contigs.fa.faa";
			while (my $line = readline($IN)) {
				print $OUT $line;
				print $OUT "\n";
			}
			close $IN;
		}
		close $OUT;

	print "Collecting consensus_contigs from each subproject...\n";

		if (-e "$baseDir/Assemblies/consensus_${project}_contigs.fa") {
			unlink "$baseDir/Assemblies/consensus_${project}_contigs.fa";
			print "Deleting previous consensus_${project}_contigs.fa ...\n";
		}
		open ($OUT, ">>$baseDir/Assemblies/consensus_${project}_contigs.fa");
		foreach my $subproject (@subprojects) {
			print "Copying consensus_${subproject}_contigs.fa to consensus_${project}_contigs.fa...\n";
			open ($IN, "$mainDir/${subproject}/Assemblies/consensus_${subproject}_contigs.fa") or die "Can't find $mainDir/${subproject}/Assemblies/consensus_${subproject}_contigs.fa";
			while (my $line = readline($IN)) {
				print $OUT $line;
				print $OUT "\n";
			}
			close $IN;
		}
		close $OUT;

		if (-e "$baseDir/Assemblies/consensus_${project}_contigs.fa.faa") {
			unlink "$baseDir/Assemblies/consensus_${project}_contigs.fa.faa";
			print "Deleting previous consensus_${project}_contigs.fa.faa ...\n";
		}
		open ($OUT, ">>$baseDir/Assemblies/consensus_${project}_contigs.fa.faa");
		foreach my $subproject (@subprojects) {
			print "Copying consensus_${subproject}_contigs.fa.faa to consensus_${project}_contigs.fa.faa...\n";
			open ($IN, "$mainDir/${subproject}/Assemblies/consensus_${subproject}_contigs.fa.faa") or die "Can't find $mainDir/${subproject}/Assemblies/consensus_${subproject}_contigs.fa.faa";
			while (my $line = readline($IN)) {
				print $OUT $line;
				print $OUT "\n";
			}
			close $IN;
		}
		close $OUT;

	if ($quant) {

		if (-d "$baseDir/Mining/salmon/consensus_${project}_contigs") {
				print "Found $baseDir/Mining/salmon/consensus_${project}_contigs...\n";
			} else {
				print "Making $baseDir/Mining/salmon/consensus_${project}_contigs...\n";
				mkdir "$baseDir/Mining/salmon/consensus_${project}_contigs";
		}
		
		print "Collecting salmon quantification of consensus_subproject_contigs from each subproject...\n";
			
		foreach my $subproject (@subprojects) {
			`scp -r $mainDir/$subproject/Mining/salmon/consensus_${subproject}_contigs/* $baseDir/Mining/salmon/consensus_${project}_contigs`;
		}
	}
}
