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

echo "DEFINE WHICH PIPELINES TO RUN"
echo "Select 1 to perform pipeline, or 0 to jump it."
echo "De-multiplexing:"; read demultiplexing_pipeline
echo "Bowtie2 Mapping:"; read bowtie2_pipeline

########################################################################
#                   DEFINE VARIABLES AND DIRECTORIES

#DE-MULTIPLEXING PIPELINE
if [ "$demultiplexing_pipeline" = "1" ]; then

#variables
echo "DE-MULTIPLEXING PIPELINE"
echo "All paths to directories must be entered as: /home/user/directory/"
echo "RAW READS MUST BE NAMED WITH EXTENTION .fastq"
echo "Path to raw reads directory:"; read rawreads_dir
echo "ATTENTION: barcodes must have size from 4 to 9pb"
echo "Path to barcodes file [ID<tab>SEQUENCE]:"; read barcodes_file
echo "How many threads to use?"; read njobs

#directories
mkdir "$rawreads_dir"fqfiles/ #FASTq FILES
mkdir "$rawreads_dir"fafiles/ #FASTa FILES
mkdir "$rawreads_dir"phix_control/  #PhiX control files

fi

#BOWTIE2 PIPELINE
if [ "$bowtie2_pipeline" = "1" ]; then

#variables
echo "Path to reference *.fasta file [extention must be named fasta (and not fas or fa)]:"; read reference
echo "Pick the preset Bowtie2 option: --very-fast-local, --fast-local, --sensitive-local, --very-sensitive-local:"; read bowtieoption

#directories
mkdir "$rawreads_dir"samfiles/ #SAM FILES
mkdir "$rawreads_dir"bamfiles/ #BAM FILES

fi
########################################################################
#                               LOG FILE

rm demultiplexing_logs.txt
touch demultiplexing_logs.txt
#use this command to add logs: echo "LOG" | tee -a demultiplexing_logs.txt

########################################################################
#                       DE-MULTIPLEXING PIPELINE

if [ "$demultiplexing_pipeline" = "1" ]; then

date | tee -a demultiplexing_logs.txt
echo "Starting de-multiplexing." | tee -a demultiplexing_logs.txt

#Split barcodes file in different files, each file merging barcodes with identical sizes
python split_barcode_bcfiles.py "$barcodes_file" "$rawreads_dir" "$njobs"
echo "BCfiles created." | tee -a demultiplexing_logs.txt

#De-multiplexing FastX-toolkit step

#define function
demultiplexing()
{
#mismatches equal 0
cat "$i" | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar9.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar9_M0_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar9_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar8.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar8_M0_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar8_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar7.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar7_M0_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar7_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar6.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar6_M0_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar6_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar5.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar5_M0_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar5_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar4.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar4_M0_ --suffix ".fq" --bol

rm "$rawreads_dir"fqfiles/"$idrun"_bar9_M0_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar8_M0_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar7_M0_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar6_M0_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar5_M0_unmatched.fq

#mismatches equal 1
cat "$rawreads_dir"fqfiles/"$idrun"_bar4_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar9.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar9_M1_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar9_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar8.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar8_M1_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar8_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar7.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar7_M1_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar7_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar6.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar6_M1_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar6_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar5.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar5_M1_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar5_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar4.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar4_M1_ --suffix ".fq" --bol

rm "$rawreads_dir"fqfiles/"$idrun"_bar4_M0_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar9_M1_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar8_M1_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar7_M1_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar6_M1_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar5_M1_unmatched.fq

#mismatches equal 2
cat "$rawreads_dir"fqfiles/"$idrun"_bar4_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar9.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar9_M2_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar9_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar8.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar8_M2_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar8_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar7.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar7_M2_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar7_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar6.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar6_M2_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar6_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar5.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar5_M2_ --suffix ".fq" --bol
cat "$rawreads_dir"fqfiles/"$idrun"_bar5_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar4.txt --prefix "$rawreads_dir"fqfiles/"$idrun"_bar4_M2_ --suffix ".fq" --bol

rm "$rawreads_dir"fqfiles/"$idrun"_bar4_M1_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar9_M2_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar8_M2_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar7_M2_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar6_M2_unmatched.fq
rm "$rawreads_dir"fqfiles/"$idrun"_bar5_M2_unmatched.fq	
}

#run demultiplexing function in parallel for each, normally four, raw reads files
idrun=0
for i in "$rawreads_dir"*.fastq
do
idrun=$(($idrun+1))
demultiplexing &
done

wait #Script will wait for all above subprocess

date | tee -a demultiplexing_logs.txt
echo "De-multiplexing done!" | tee -a demultiplexing_logs.txt

#Join files to get individual complete files
cat "$rawreads_dir"list_individuals.txt | while read line
do
cat "$rawreads_dir"fqfiles/*"$line".fq > "$rawreads_dir"fqfiles/joined_"$line".fastq  #individuals final files
done

echo "FASTq files joined." | tee -a demultiplexing_logs.txt

rm "$rawreads_dir"fqfiles/*bar*       #remove all temporary fastq files

echo "All temporary files were removed." | tee -a demultiplexing_logs.txt

#-----------------------------------------------------------------------
#PHIX CONTROL

echo "Starting PhiX control." | tee -a demultiplexing_logs.txt

#Create Local Blast DB for PhiX genome
makeblastdb -in phix_genome.fasta -dbtype nucl -out "$rawreads_dir"phix_control/phix.blastdb
echo "BLAST database for PhiX genome created." | tee -a demultiplexing_logs.txt

#FASTq to FASTA

fastqtofasta()
{
cat "$i" | while read line
do
echo "Converting FASTq $line to FASTA." | tee -a demultiplexing_logs.txt
fastq_to_fasta -Q33 -i "$rawreads_dir"fqfiles/joined_"$line".fastq -o "$rawreads_dir"fafiles/joined_"$line".fasta  #convert fastq to fasta using phred33 quality
done
}

for i in "$rawreads_dir"*.tmpp
do
fastqtofasta &  #works in parallel for njobs set
done

wait

date | tee -a demultiplexing_logs.txt
echo "All FASTq files converted to FASTA." | tee -a demultiplexing_logs.txt

#BLAST searches

blastsearch()
{
cat "$i" | while read line
do
echo "Blasting $line against PhiX." | tee -a demultiplexing_logs.txt
blastn -query "$rawreads_dir"fafiles/joined_"$line".fasta -db "$rawreads_dir"phix_control/phix.blastdb -evalue 0.001 -outfmt 6 -out "$rawreads_dir"phix_control/"$line".blastn   #do blastn search
done
}

for i in "$rawreads_dir"*.tmpp
do
blastsearch &   #works in parallel for njobs set
done

wait

date | tee -a demultiplexing_logs.txt
echo "Blast searches has finished." | tee -a demultiplexing_logs.txt

#Filter out PhiX reads from individual FASTq files
phixfilter()
{
cat "$i" | while read line
do
echo "Filtering PhiX reads from $line." | tee -a demultiplexing_logs.txt
python phix_filter.py "$rawreads_dir"fqfiles/joined_"$line".fastq "$rawreads_dir"phix_control/"$line".blastn > "$rawreads_dir"fqfiles/nophix_"$line".fastq
done
}

for i in "$rawreads_dir"*.tmpp
do
phixfilter &   #works in parallel for njobs set
done

wait

date | tee -a demultiplexing_logs.txt
echo "All PhiX filtered." | tee -a demultiplexing_logs.txt

echo "Removing all extra joined FASTq files. Keeping no-phix FASTq only." | tee -a demultiplexing_logs.txt
rm "$rawreads_dir"fqfiles/joined*.fastq

#rm "$rawreads_dir"*.tmpp
#rm "$rawreads_dir"list_individuals.txt

fi

########################################################################
#                       BOWTIE2 MAPPING PIPELINE

if [ "$bowtie2_pipeline" = "1" ]; then

#SAMTOOLS: INDEX REFERENCE FASTA FILE

echo "Starting mapping process." | tee -a demultiplexing_logs.txt
echo "Indexing reference FASTA file." | tee -a demultiplexing_logs.txt
nameindex=`echo "$reference" | sed "s/.fasta/_index/g"`
bowtie2-build "$reference" "$nameindex"

#SAMTOOLS: THE ALIGNMENT

bowtie2mapping()
{
cat "$i" | while read line
do
echo "Mapping $line to reference." | tee -a demultiplexing_logs.txt
bowtie2 "$bowtieoption" --phred33 --trim5 4 --trim3 10 -x "$nameindex" -U "$rawreads_dir"fqfiles/nophix_"$line".fastq -S "$rawreads_dir"samfiles/"$line".sam
done
}

for i in "$rawreads_dir"*.tmpp
do
bowtie2mapping &   #works in parallel for njobs set
done

wait

date | tee -a demultiplexing_logs.txt
echo "Mapping process done." | tee -a demultiplexing_logs.txt

#-----------------------------------------------------------------------
#STATISTCS

echo "Starting SAMtools calculations." | tee -a demultiplexing_logs.txt

#SAMTOOLS CALCULATIONS

samtoolscalc()
{
cat "$i" | while read line
do
echo "SAMtools calculations for $line" | tee -a demultiplexing_logs.txt

samtools view -bS "$rawreads_dir"samfiles/"$line".sam > "$rawreads_dir"bamfiles/"$line".bam   #convert SAM to BAM
samtools sort "$rawreads_dir"bamfiles/"$line".bam "$rawreads_dir"bamfiles/"$line"_sorted      #sort BAM
samtools flagstat "$rawreads_dir"bamfiles/"$line".bam > "$rawreads_dir"bamfiles/"$line".map   #mapping percentage
samtools depth "$rawreads_dir"bamfiles/"$line"_sorted.bam > "$rawreads_dir"bamfiles/"$line"_sorted.coverage #depth ave
done
}

for i in "$rawreads_dir"*.tmpp
do
samtoolscalc &   #works in parallel for njobs set
done

wait

date | tee -a demultiplexing_logs.txt
echo "All SAMtools calculation is done." | tee -a demultiplexing_logs.txt

echo "Creating tables." | tee -a demultiplexing_logs.txt

#CREATE TABLES
cat "$rawreads_dir"list_individuals.txt | while read line
do

#mapped reads table
number=`cat "$rawreads_dir"bamfiles/"$line".map | grep mapped | head -1 | sed 's/^.*(//g' | sed 's/%:-.*//g'`
echo "$line $number" | tee -a "$rawreads_dir"mappedreads.txt

#depth average table
coveragefile=`"$rawreads_dir"bamfiles/"$line"_sorted.coverage`
sed -i "s/myave <- read.table (\".*\")/myave <- read.table \(\"$coveragefile\"\)/g" calc_depth_average.R
ave=`Rscript calc_depth_average.R`
echo "$i $ave" | tee -a "$rawreads_dir"coverage.txt

date | tee -a demultiplexing_logs.txt
echo "Tables have been created." | tee -a demultiplexing_logs.txt

done

fi

date | tee -a demultiplexing_logs.txt
echo "De-multiplexing Master has finished." | tee -a demultiplexing_logs.txt
