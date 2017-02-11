#!/bin/bash
#Created by Hugo Rody, fev 11 of 2017
#Federal University of Sao Paulo, Brazil
#Pipeline for de-multiplexing Illumina Next Seq GBS raw data and High output kit, applying PhiX174 control

########################################################################
#DEPENDENCIES

#SAMTOOLS
#FASTX-TOOLKIT
#BLAST++

#VERIFY DEPENDENCIES

blast=`which blastn`
makeblast=`which makeblastdb`
samtools=`which samtools`
fastq=`which fastq_to_fasta`

if [ "$blast" = "" ]; then
echo "BLAST is not installed"
fi

if [ "$makeblast" = "" ]; then
echo "MAKEBLASTDB is not installed"
fi

if [ "$samtools" = "" ]; then
echo "SAMTOOLS is not installed"
fi

if [ "$fastq" = "" ]; then
echo "FASTX-toolkit is not installed"
fi


########################################################################
#DEFINE VARIABLES

#demultiplexing
echo "DE-MULTIPLEXING PIPELINE"
echo "All paths to directories must be entered as: /home/user/directory/"
echo "RAW READS MUST BE NAMED WITH EXTENTION .fastq"
echo "Path to raw reads directory:"; read rawreads_dir
echo "ATTENTION: barcodes must have size from 4 to 9pb"
echo "Path to barcodes file [ID<tab>SEQUENCE]:"; read barcodes_file
echo "Set an ID for your run:"; read idrun

#alignment
echo "Path to reference *.fasta file [extention must be named fasta (and not fas or fa)]:"; read reference
echo "Pick the preset Bowtie2 option:\n--very-fast-local\n--fast-local\n--sensitive-local\n--very-sensitive-local\n"; read bowtieoption


########################################################################
#CREATES DIRECTORIES

mkdir "$rawreads_dir"fqfiles/ #FASTq FILES
mkdir "$rawreads_dir"fafiles/ #FASTa FILES
mkdir "$rawreads_dir"samfiles/ #SAM FILES
mkdir "$rawreads_dir"bamfiles/ #BAM FILES
mkdir "$rawreads_dir"phix_control/  #PhiX control files


########################################################################

#FASTX-TOOLKIT: DEMULTIPLEXING

#Split barcodes file into different files, each file merging barcodes with identical sizes
python split_barcode_bcfiles.py "$barcodes_file" "$rawreads_dir"

#De-multiplexing step
for i in "$rawreads_dir"*.fastq
do

#mismatches equal 0
cat "$i" | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar9.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar9_M0_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar9_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar8.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar8_M0_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar8_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar7.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar7_M0_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar7_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar6.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar6_M0_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar6_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar5.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar5_M0_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar5_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 0 --bcfile "$rawreads_dir"bar4.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar4_M0_ --suffix ".fq" --bol

rm "$rawreads_dir"/fqfiles/"$idrun"_bar9_M0_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar8_M0_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar7_M0_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar6_M0_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar5_M0_unmatched.fq

#mismatches equal 1
cat "$rawreads_dir"/fqfiles/"$idrun"_bar4_M0_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar9.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar9_M1_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar9_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar8.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar8_M1_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar8_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar7.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar7_M1_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar7_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar6.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar6_M1_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar6_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar5.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar5_M1_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar5_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 1 --bcfile "$rawreads_dir"bar4.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar4_M1_ --suffix ".fq" --bol

rm "$rawreads_dir"/fqfiles/"$idrun"_bar4_M0_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar9_M1_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar8_M1_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar7_M1_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar6_M1_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar5_M1_unmatched.fq

#mismatches equal 2
cat "$rawreads_dir"/fqfiles/"$idrun"_bar4_M1_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar9.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar9_M2_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar9_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar8.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar8_M2_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar8_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar7.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar7_M2_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar7_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar6.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar6_M2_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar6_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar5.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar5_M2_ --suffix ".fq" --bol
cat "$rawreads_dir"/fqfiles/"$idrun"_bar5_M2_unmatched.fq | perl fastx_barcode_splitter.pl --mismatches 2 --bcfile "$rawreads_dir"bar4.txt --prefix "$rawreads_dir"/fqfiles/"$idrun"_bar4_M2_ --suffix ".fq" --bol

rm "$rawreads_dir"/fqfiles/"$idrun"_bar4_M1_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar9_M2_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar8_M2_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar7_M2_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar6_M2_unmatched.fq
rm "$rawreads_dir"/fqfiles/"$idrun"_bar5_M2_unmatched.fq

done

#Join files to get individual complete files

cat "$rawreads_dir"list_individuals.txt | sed "s/^[a-z]//g" | sort | uniq | while read line
do

cat "$rawreads_dir"/fqfiles/*"$line".fq > "$rawreads_dir"/fqfiles/joined_"$line".fastq  #individuals final files

done

rm "$rawreads_dir"/fqfiles/*bar*       #remove all temporary fastq files

########################################################################
#PHIX CONTROL

#Create Local Blast DB
makeblastdb -in phix_genome.fasta -dbtype nucl -out "$rawreads_dir"phix_control/phix.blastdb

#Blast searches
cat "$rawreads_dir"list_individuals.txt | sed "s/^[a-z]//g" | sort | uniq | while read line
do

fastq_to_fasta -Q33 -i "$rawreads_dir"/fqfiles/joined_"$line".fastq -o "$rawreads_dir"/fafiles/joined_"$line".fasta  #convert fastq to fasta using phred33 quality

blastn -query "$rawreads_dir"/fafiles/joined_"$line".fasta -db "$rawreads_dir"phix_control/phix.blastdb -evalue 0.001 -outfmt 6 -out "$rawreads_dir"phix_control/"$line".blastn

done


#Filter out PhiX from individual FASTq files
cat "$rawreads_dir"list_individuals.txt | sed "s/^[a-z]//g" | sort | uniq | while read line
do

python phix_filter.py "$rawreads_dir"/fqfiles/joined_"$line".fastq "$rawreads_dir"phix_control/"$line".blastn > "$rawreads_dir"/fqfiles/nophix_"$line".fastq

done


########################################################################
#SAMTOOLS: INDEX REFERENCE FASTA FILE

nameindex=`echo "$reference" | sed "s/.fasta/_index/g"`

bowtie2-build "$reference" "$nameindex"



#SAMTOOLS: DO THE ALIGNMENT

cat "$rawreads_dir"list_individuals.txt | sed "s/^[a-z]//g" | sort | uniq | while read line
do

bowtie2 "$bowtieoption" --phred33 --trim5 4 --trim3 10 -x "$nameindex" -U "$rawreads_dir"/fqfiles/nophix_"$line".fastq -S "$rawreads_dir"samfiles/"$line".sam

done

########################################################################
#Statistics

#SAMTOOLS CALCULATIONS
cat "$rawreads_dir"list_individuals.txt | sed "s/^[a-z]//g" | sort | uniq | while read line
do

name=`echo "$i" | sed "s/.sam//g"`

echo "$name"

samtools view -bS "$rawreads_dir"samfiles/"$line".sam > "$rawreads_dir"bamfiles/"$line".bam   #convert SAM to BAM
samtools sort "$rawreads_dir"bamfiles/"$line".bam "$rawreads_dir"bamfiles/"$line"_sorted.bam  #sort BAM

samtools flagstat "$rawreads_dir"bamfiles/"$line".bam > "$rawreads_dir"bamfiles/"$line".map   #mapping percentage
samtools depth "$rawreads_dir"bamfiles/"$line"_sorted.bam > "$rawreads_dir"bamfiles/"$line"_sorted.coverage #depth ave

done


#GENERATE TABLES
cat "$rawreads_dir"list_individuals.txt | sed "s/^[a-z]//g" | sort | uniq | while read line
do

#mapped reads
number=`cat "$rawreads_dir"bamfiles/"$line".map | grep mapped | head -1 | sed 's/^.*(//g' | sed 's/%:-.*//g'`
echo "$line $number" | tee -a "$rawreads_dir"mappedreads.txt

#depth ave

coveragefile=`"$rawreads_dir"bamfiles/"$line"_sorted.coverage`

sed -i "s/myave <- read.table (\".*\")/myave <- read.table \(\"$coveragefile\"\)/g" calc_depth_average.R

ave=`Rscript calc_depth_average.R`

echo "$i $ave" | tee -a "$rawreads_dir"coverage.txt


done
