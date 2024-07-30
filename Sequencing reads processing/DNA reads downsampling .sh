#!/bin/bash
for i in `ls *bam`
do
  samtools view -@ 48 -b -q 59 $i > temp.bam
  samtools sort -@ 48 temp.bam -o ${i%.bam}_sorted.bam
  ~/software/gatk-4.2.0.0/gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=24" MarkDuplicates I=${i%.bam}_sorted.bam O=${i%.bam}_dedup_sorted.bam M=ASXS.metrix REMOVE_DUPLICATES=true ASSUME_SORTED=true
  samtools index -@ 48 ${i%.bam}_dedup_sorted.bam
  echo ${i%.bam}_dedup_sorted.bam >> bwa.txt
  samtools flagstat $i >> bwa.txt
done

for i in `ls *_dedup_sorted.bam`
do
  #samtools index -@ 24 $i
  fraction=$(samtools idxstats $i | cut -f3 | awk -v ct=5000000 'BEGIN {total=0} {total += $1} END {print ct/total}')
  st=$(echo "$fraction < 1.0" | bc)
  if [ $st -eq 1 ]; then
    samtools view -@ 48 -b -s ${fraction} $i > ./downsampling/${i%_dedup_sorted.bam}_5000000_sampled.bam
  fi
done

for i in `ls *5000000_sampled.bam`
do
  samtools index $i
done
for i in `ls *bam`
do
  bamToBed -i $i | pigz > ${i%.bam}.bed.gz
done
