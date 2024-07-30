#!/bin/bash
## RNA quantification
cd ~/OneDrop/OD_20231010/lane_03
for f in `ls trimmed_*read1*`
do
  newname=${f#trimmed_}
  newname="lane03_$newname"
  echo $newname
  kallisto quant --single -l 200 -s 30 -t 24 --bias -i ~/reference/humanRNA/kallisto_ref/gcv29_rRNA_ebv_ercc.idx -o ~/OneDrop/OD_20231010/mapping/kallisto/${newname%_read1.fastq.gz}_human $f
  kallisto quant --single -l 200 -s 30 -t 24 --bias -i ~/reference/mourse/kallisto_mm27/K_GC_vM27_ERCC.idx -o ~/OneDrop/OD_20231010/mapping/kallisto/${newname%_read1.fastq.gz}_mouse $f

  STAR --runThreadN 18 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --genomeDir ~/reference/humanRNA/star_from10x_ercc --readFilesCommand zcat --readFilesIn $f --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/OneDrop/OD_20231010/mapping/star/${newname%_read1.fastq.gz}/

  bwa mem -t 24 -M ~/reference/Ecoli/bwa_ecoli/GCA_000019425.1_ASM1942v1_DH10B.fna $f > temp.sam
  samtools view -@ 24 -S -b -q 40 temp.sam > ~/OneDrop/OD_20231010/mapping/bwa/${newname%_read1.fastq.gz}_ecoli.bam
  rm temp.sam

  bwa mem -t 24 -M  ~/reference/mourse/BWA_mm10/GRCm38.primary_assembly.genome.fa $f > temp.sam
  samtools view -@ 24 -S -b -q 40 temp.sam > ~/OneDrop/OD_20231010/mapping/bwa/${newname%_read1.fastq.gz}_mouse.bam
  rm temp.sam

  bwa mem -t 24 -M  ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna $f > temp.sam
  samtools view -@ 24 -S -b -q 40 temp.sam > ~/OneDrop/OD_20231010/mapping/bwa/${newname%_read1.fastq.gz}_human.bam
  rm temp.sam
done
