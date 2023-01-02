#This code is adapted from what Tyler sent me, plus I played around with the settings a bit
#The allelic series reads live here: /2/scratch/dworkin_temp/allelicSeries/rawReads

#name=/allelicSeries/rawReads/SAM_sdETX4_r1_ACAGTG_L005_R2_001.fastq



raw_dir=/2/scratch/amandaN/GBE/unzipped_reads
trim_dir=/2/scratch/amandaN/GBE/bbduk_trimmed

files=${raw_dir}/*.fastq

#First, it seems that I cannot use zipped files. So unzip those (delete later)

#For loop over every file
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .fastq`
base=${base%??????}
/usr/local/BBmap/bbduk.sh \
in1=${raw_dir}/${base}R1_001.fastq \
in2=${raw_dir}/${base}R2_001.fastq \
out1=${trim_dir}/${base}R1_001_trim-bbduk.fastq \
out2=${trim_dir}/${base}R2_001_trim-bbduk.fastq \
ref=/2/scratch/amandaN/cornell_seqdata/BBMap_adapters.fa \
threads=8 ktrim=r k=23 mink=10 hdist=1 tpe tbo \
qtrim=rl trimq=15 minlength=36
done
