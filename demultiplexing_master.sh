#!/bin/sh
#Created by Hugo Rody, fev 11 of 2017
#Federal University of Sao Paulo, Brazil
#Pipeline for de-multiplexing Illumina Next Seq (High output kit) GBS raw data, apply PhiX174 control, and perform Bowtie2 mapping
#usage: sh demultiplexing_master.sh

########################################################################
#                             DEPENDENCIES

#Bowtie2
#SAMTOOLS
#FASTX-TOOLKIT
#BLAST++

#VERIFY DEPENDENCIES

blastdep=`which blastn`
bowtie2dep=`which bowtie2`
makeblastdep=`which makeblastdb`
samtoolsdep=`which samtools`
fastqdep=`which fastq_to_fasta`

if [ "$blastdep" = "" ]; then
    echo "BLAST is not installed"
    break
fi

if [ "$makeblastdep" = "" ]; then
    echo "MAKEBLASTDB is not installed"
    break
fi

if [ "$samtoolsdep" = "" ]; then
    echo "SAMTOOLS is not installed"
    break
fi

if [ "$fastqdep" = "" ]; then
    echo "FASTX-toolkit is not installed"
    break
fi

if [ "$bowtie2dep" = "" ]; then
    echo "Bowtie2 is not installed"
    break
fi

########################################################################
#                    DEFINE WHICH PIPELINES TO RUN

echo -e "-------------------------------------------------------------\nDEFINE WHICH PIPELINES TO RUN\n-------------------------------------------------------------"
echo "Select 1 to perform pipeline, or 0 to jump it."
echo "De-multiplexing:"; read demultiplexing_pipeline
echo "Bowtie2 Mapping:"; read bowtie2_pipeline
echo "SAMtools statistics:"; read samtools_statistics

########################################################################
#                   DEFINE VARIABLES AND DIRECTORIES

#DE-MULTIPLEXING PIPELINE
if [ "$demultiplexing_pipeline" = "1" ]; then
    #variables
    echo -e "-------------------------------------------------------------\nDE-MULTIPLEXING PIPELINE\n-------------------------------------------------------------"
    echo -e "All paths to directories must be entered as: /home/user/directory/"
    echo -e "RAW READS MUST BE NAMED WITH EXTENTION .fastq"
    echo -e "-------------------------------------------------------------\nPath to raw reads directory:"; read rawreads_dir
    echo -e "\nATTENTION: barcodes must have size from 4 to 9pb\n"
    echo -e "\nPath to barcodes file [ID<tab>SEQUENCE]:"; read barcodes_file
	echo -e "\nHow many threads to use?"; read njobs
    echo -e "\nPath to output directory [you must have already created it]:"; read outputdir

    #directories
    mkdir "$outputdir"fqfiles/ #FASTq FILES
    mkdir "$outputdir"fafiles/ #FASTa FILES
    mkdir "$outputdir"phix_control/  #PhiX control files
    mkdir "$outputdir"fqjoined/
fi

#BOWTIE2 PIPELINE
if [ "$bowtie2_pipeline" = "1" ]; then

    #variables
    echo -e "-------------------------------------------------------------\nBOWTIE 2 PIPELINE\n-------------------------------------------------------------"
    echo "Path to reference *.fasta file [extention must be named fasta (and not fas or fa)]:"; read reference
    echo "Pick the preset Bowtie2 option: --very-fast-local, --fast-local, --sensitive-local, --very-sensitive-local:"; read bowtieoption

    #verify if outputdir has been aready created
    if [ "$outputdir" = "" ]; then
        echo "Path to output directory [you must have already created it]:"; read outputdir
    fi

    #directories
    mkdir "$outputdir"samfiles/ #SAM FILES
    mkdir "$outputdir"bamfiles/ #BAM FILES

fi

#STATISTICS PIPELINE
if [ "$samtools_statistics" = "1" ]; then

    if [ "$outputdir" = "" ]; then
        echo "Path to output directory [you must have already created it]:"; read outputdir
    fi

fi

########################################################################
#                               LOG FILE

rm "$outputdir"demultiplexing_logs.txt
touch "$outputdir"demultiplexing_logs.txt
#use this command to add logs: echo "LOG" | tee -a demultiplexing_logs.txt

########################################################################
#                       DE-MULTIPLEXING PIPELINE

if [ "$demultiplexing_pipeline" = "1" ]; then

    date | tee -a "$outputdir"demultiplexing_logs.txt
    echo "Starting de-multiplexing." | tee -a "$outputdir"demultiplexing_logs.txt

    #Split barcodes file in different files, each file merging barcodes with identical sizes
    python split_barcode_bcfiles.py "$barcodes_file" "$outputdir" "$njobs"
    echo "BarCode files created." | tee -a "$outputdir"demultiplexing_logs.txt
    
    #-----------------------------------------------------------------------
    #De-multiplexing FastX-toolkit step

    #define function
	demultiplexing()
	{
	#mismatches equal 0
	cat "$i" | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$outputdir"bar9.txt --prefix "$outputdir"fqfiles/"$idrun"_bar9_M0_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar9_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$outputdir"bar8.txt --prefix "$outputdir"fqfiles/"$idrun"_bar8_M0_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar8_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$outputdir"bar7.txt --prefix "$outputdir"fqfiles/"$idrun"_bar7_M0_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar7_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$outputdir"bar6.txt --prefix "$outputdir"fqfiles/"$idrun"_bar6_M0_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar6_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$outputdir"bar5.txt --prefix "$outputdir"fqfiles/"$idrun"_bar5_M0_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar5_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$outputdir"bar4.txt --prefix "$outputdir"fqfiles/"$idrun"_bar4_M0_ --suffix ".fq" --bol
	
	rm "$outputdir"fqfiles/"$idrun"_bar9_M0_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar8_M0_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar7_M0_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar6_M0_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar5_M0_unmatched.fq
	
	#mismatches equal 1
	cat "$outputdir"fqfiles/"$idrun"_bar4_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$outputdir"bar9.txt --prefix "$outputdir"fqfiles/"$idrun"_bar9_M1_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar9_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$outputdir"bar8.txt --prefix "$outputdir"fqfiles/"$idrun"_bar8_M1_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar8_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$outputdir"bar7.txt --prefix "$outputdir"fqfiles/"$idrun"_bar7_M1_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar7_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$outputdir"bar6.txt --prefix "$outputdir"fqfiles/"$idrun"_bar6_M1_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar6_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$outputdir"bar5.txt --prefix "$outputdir"fqfiles/"$idrun"_bar5_M1_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar5_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$outputdir"bar4.txt --prefix "$outputdir"fqfiles/"$idrun"_bar4_M1_ --suffix ".fq" --bol
	
	rm "$outputdir"fqfiles/"$idrun"_bar4_M0_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar9_M1_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar8_M1_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar7_M1_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar6_M1_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar5_M1_unmatched.fq
	
	#mismatches equal 2
	cat "$outputdir"fqfiles/"$idrun"_bar4_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$outputdir"bar9.txt --prefix "$outputdir"fqfiles/"$idrun"_bar9_M2_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar9_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$outputdir"bar8.txt --prefix "$outputdir"fqfiles/"$idrun"_bar8_M2_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar8_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$outputdir"bar7.txt --prefix "$outputdir"fqfiles/"$idrun"_bar7_M2_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar7_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$outputdir"bar6.txt --prefix "$outputdir"fqfiles/"$idrun"_bar6_M2_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar6_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$outputdir"bar5.txt --prefix "$outputdir"fqfiles/"$idrun"_bar5_M2_ --suffix ".fq" --bol
	cat "$outputdir"fqfiles/"$idrun"_bar5_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$outputdir"bar4.txt --prefix "$outputdir"fqfiles/"$idrun"_bar4_M2_ --suffix ".fq" --bol
	
	rm "$outputdir"fqfiles/"$idrun"_bar4_M1_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar9_M2_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar8_M2_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar7_M2_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar6_M2_unmatched.fq
	rm "$outputdir"fqfiles/"$idrun"_bar5_M2_unmatched.fq	
	}

    #run demultiplexing function in parallel for each, normally four, raw reads files
    idrun=0
    for i in "$rawreads_dir"*.fastq
    do
        idrun=$(($idrun+1))
	    echo "Demultplexing run $idrun started." | tee -a "$outputdir"demultiplexing_logs.txt
        demultiplexing &
    done

    wait #Script will wait for all above subprocess

    date | tee -a "$outputdir"demultiplexing_logs.txt
    echo "De-multiplexing done!" | tee -a "$outputdir"demultiplexing_logs.txt

    #Join files to get individual complete files
	date
	echo "Merging individuals FASTQ files..." | tee -a "$outputdir"demultiplexing_logs.txt
    python merge_fastq_demultiplexing.py "$outputdir"list_individuals.txt "$outputdir"fqfiles/ #individuals final files
    echo "FASTq files joined." | tee -a "$outputdir"demultiplexing_logs.txt
	date | tee -a "$outputdir"demultiplexing_logs.txt

    #rm "$outputdir"fqfiles/*bar*       #remove all temporary fq files

	date | tee -a "$outputdir"demultiplexing_logs.txt
    echo "All temporary files were removed." | tee -a "$outputdir"demultiplexing_logs.txt

    #-----------------------------------------------------------------------
    #PHIX CONTROL

    echo "Starting PhiX control." | tee -a "$outputdir"demultiplexing_logs.txt

    #Create Local Blast DB for PhiX genome
    makeblastdb -in phix_genome.fasta -dbtype nucl -out "$outputdir"phix_control/phix.blastdb
    echo "BLAST database for PhiX genome created." | tee -a "$outputdir"demultiplexing_logs.txt

    #Convert FASTq to FASTA to perform BLAST

	fastqtofasta()
	{
	cat "$i" | while read line
	do
		echo "Converting FASTq $line to FASTA." | tee -a "$outputdir"demultiplexing_logs.txt
		fastq_to_fasta -Q33 -i "$outputdir"fqfiles/joined_"$line".fastq -o "$outputdir"fafiles/joined_"$line".fasta  #convert fastq to fasta using phred33 quality
	done
	}
	
	for i in "$outputdir"*.tmpp
	do
		fastqtofasta &  #works in parallel for njobs set
	done
	
	wait
	
	date | tee -a "$outputdir"demultiplexing_logs.txt
	echo "All FASTq files converted to FASTA." | tee -a "$outputdir"demultiplexing_logs.txt
	
	#BLAST searches
	
	blastsearch()
	{
	cat "$i" | while read line
	do
		echo "Blasting $line against PhiX." | tee -a "$outputdir"demultiplexing_logs.txt
		blastn -query "$outputdir"fafiles/joined_"$line".fasta -db "$outputdir"phix_control/phix.blastdb -evalue 0.001 -outfmt 6 -out "$outputdir"phix_control/"$line".blastn   #do blastn search
	done
	}
	
	for i in "$outputdir"*.tmpp
	do
		blastsearch &   #works in parallel for njobs set
	done
	
	wait
	
	date | tee -a "$outputdir"demultiplexing_logs.txt
	echo "Blast searches has finished." | tee -a "$outputdir"demultiplexing_logs.txt

	#Filter out PhiX reads from individual FASTq files
	phixfilter()
	{
	cat "$i" | while read line
	do
		echo "Filtering PhiX reads from $line." | tee -a "$outputdir"demultiplexing_logs.txt
		python phix_filter.py "$outputdir"fqfiles/joined_"$line".fastq "$outputdir"phix_control/"$line".blastn > "$outputdir"fqfiles/nophix_"$line".fastq
	done
	}
	
	for i in "$outputdir"*.tmpp
	do
		phixfilter &   #works in parallel for njobs set
	done

    wait

    date | tee -a "$outputdir"demultiplexing_logs.txt
    echo "All PhiX filtered." | tee -a "$outputdir"demultiplexing_logs.txt

    echo "Moving all extra joined FASTq files to fqjoined directory. Keeping no-phix FASTq only." | tee -a "$outputdir"demultiplexing_logs.txt
    mv "$outputdir"fqfiles/joined*.fastq "$outputdir"fqjoined/

fi

########################################################################
#                       BOWTIE2 MAPPING PIPELINE

if [ "$bowtie2_pipeline" = "1" ]; then

	#SAMTOOLS: INDEX REFERENCE FASTA FILE
	
	echo "Starting mapping process." | tee -a "$outputdir"demultiplexing_logs.txt
	echo "Indexing reference FASTA file." | tee -a "$outputdir"demultiplexing_logs.txt
	nameindex=`echo "$reference" | sed "s/.fasta/_index/g"`
	bowtie2-build "$reference" "$nameindex"
	
	#SAMTOOLS: THE ALIGNMENT
	
	bowtie2mapping()
	{
	cat "$i" | while read line
	do
		echo "Mapping $line to reference." | tee -a "$outputdir"demultiplexing_logs.txt
		bowtie2 "$bowtieoption" --phred33 --trim5 4 --trim3 10 -x "$nameindex" -U "$outputdir"fqfiles/nophix_"$line".fastq -S "$outputdir"samfiles/"$line".sam
	done
	}
	
	for i in "$outputdir"*.tmpp
	do
		bowtie2mapping &   #works in parallel for njobs set
	done
	
	wait
	
	date | tee -a "$outputdir"demultiplexing_logs.txt
	echo "Mapping process done." | tee -a "$outputdir"demultiplexing_logs.txt

fi
########################################################################
#                         SAMTOOLS STATISTCS

if [ "$samtools_statistics" = "1" ]; then

	echo "Starting SAMtools calculations." | tee -a "$outputdir"demultiplexing_logs.txt
	
	#SAMTOOLS CALCULATIONS
	
	samtoolscalc()
	{
	cat "$i" | while read line
	do
		echo "SAMtools calculations for $line" | tee -a "$outputdir"demultiplexing_logs.txt
		
		samtools view -bS "$outputdir"samfiles/"$line".sam > "$outputdir"bamfiles/"$line".bam   #convert SAM to BAM
		samtools sort "$outputdir"bamfiles/"$line".bam > "$outputdir"bamfiles/"$line".bam_sorted      #sort BAM
		samtools flagstat "$outputdir"bamfiles/"$line".bam > "$outputdir"bamfiles/"$line".map   #mapping percentage
		samtools depth "$outputdir"bamfiles/"$line".bam_sorted > "$outputdir"bamfiles/"$line"_sorted.coverage #depth ave
	done
	}
	
	for i in "$outputdir"*.tmpp
	do
		samtoolscalc &   #works in parallel for njobs set
	done
	
	wait
	
	date | tee -a "$outputdir"demultiplexing_logs.txt
	echo "All SAMtools calculation is done." | tee -a "$outputdir"demultiplexing_logs.txt
	
	echo "Creating tables." | tee -a "$outputdir"demultiplexing_logs.txt
	
	#CREATE TABLES
	cat "$outputdir"list_individuals.txt | while read line
	do
	
		#mapped reads table
		number=`cat "$outputdir"bamfiles/"$line".map | grep "mapped" | head -1 | sed "s/^.*(//g" | sed "s/\% :.*//g"`
		
		echo "$line $number" | tee -a "$outputdir"mappedreads.txt
		
		#depth average table
		coveragefile=`echo "$outputdir"bamfiles/"$line"_sorted.coverage | sed "s/\//+/g"`
		
		sed -i 's/myave <- read.table (".*")/myave <- read.table ("'"$coveragefile"'")/g' calc_depth_average.R
		sed -i "s/+/\//g" calc_depth_average.R
		
		ave=`Rscript calc_depth_average.R`
		
		echo "$line $ave" | tee -a "$outputdir"coverage.txt
		
		#date | tee -a "$outputdir"demultiplexing_logs.txt
		echo "Tables have been created." | tee -a "$outputdir"demultiplexing_logs.txt
	
	done
	
	#date | tee -a "$outputdir"demultiplexing_logs.txt
	echo "De-multiplexing Master has finished." | tee -a "$outputdir"demultiplexing_logs.txt
fi
