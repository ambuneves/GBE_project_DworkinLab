#####################
#First set up some shell variables that we will want for directories, and files that do not change.
genome=/2/scratch/amandaN/cornell_seqdata/STAR_index

# Note that (I assume for speed) you can keep the reference genome loaded
# into memory with STAR. So we want to do this before we run the loop,
# and then turn it off after we finish mapping everything.

#  load genome for shared use.

genome=/2/scratch/amandaN/cornell_seqdata/STAR_index

STAR --genomeLoad LoadAndExit \
  --genomeDir ${genome}

file_dir=/2/scratch/amandaN/GBE/bbduk_trimmed/unfinished
mapped_dir=/2/scratch/amandaN/GBE/star_quant
gtf_file=/2/scratch/amandaN/index_files/dmel-all-r6.38.gtf


files=${file_dir}/*.fastq

# Aliging
for file in ${files[@]}
do
  name=${file}
  base=`basename ${name} .fastq`
  base=${base%?????????????????}
  STAR --runThreadN 16 \
  --genomeDir ${genome} \
  --readFilesIn ${file_dir}/${base}R1_001_trim-bbduk.fastq ${file_dir}/${base}R2_001_trim-bbduk.fastq \
  --outFileNamePrefix ${mapped_dir}/${base}_bbduktrim \
  --sjdbGTFfile ${gtf_file} \
  --quantMode TranscriptomeSAM GeneCounts
done
