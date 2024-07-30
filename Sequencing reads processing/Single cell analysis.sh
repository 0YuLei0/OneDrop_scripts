############################# RNA PART
#!/bin/bash
for i in `ls lane_0*/*.fastq.gz`
do
  a=${i:16}
  a=${a%_CB.fastq.gz}
  STAR --runThreadN 18 --genomeDir ~/reference/humanRNA/star_from10x_ercc --readFilesCommand zcat \
  --soloCBwhitelist [cell barcode file] --soloMultiMappers Rescue --soloType CB_UMI_Complex \
  --readFilesIn read1.fastq.gz ---soloAdapterMismatchesNmax 3 --soloCBmatchWLtype EditDist_2 \
  --soloCBposition [cell barcode posisiton] --soloUMIposition [UMI posisiton] --soloFeatures Gene GeneFull Velocyto \
  --outSAMtype BAM SortedByCoordinate --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --outFileNamePrefix ./mapping/scRNA/
done
############################# DNA PART
for i in `ls lane_0*/*CHIP*DNA*Beads*CB.fastq.gz`
do
  a=${i:16}
  a=${a%_CB.fastq.gz}
  mkdir -p ./mapping/scDNA/$a

  umi_tools whitelist -I $i -S OUT_TSV.gz --plot-prefix=plots_ --error-correct-threshold=2 \
  --method=reads --extract-method=regex --bc-pattern= [cell barcode parrtern] \
  --read2-in=${i%CB.fastq.gz}read1.fastq.gz -L ./mapping/scDNA/$a/extract.log -S ./mapping/scDNA/$a/denovo.whitelist

  awk '{print $1,$1}' ./mapping/scDNA/$a/denovo.whitelist > ./mapping/scDNA/$a/cells.txt

  umi_tools extract -I $i -S ./mapping/scDNA/${a}.CB.fastq.gz \
  --whitelist=denovo.whitelist --bc-pattern=[cell barcode parrtern] \
  --read2-in=${i%CB.fastq.gz}read1.fastq.gz --read2-out=./mapping/scDNA/${a}.read1.fastq.gz -L ./mapping/scDNA/$a/extract.log --umi-separator="_"

  bwa mem -t 24 -M  ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna CHIP-Mock_DNA-Beads_read1.fastq.gz \
  | samtools view -@ 48 -b -q 59 > ./mapping/scDNA/$a/temp.bam

  samtools sort -@ 48 ./mapping/scDNA/$a/temp.bam -o ./mapping/scDNA/$a/sorted.bam

  ~/software/gatk-4.2.0.0/gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=24" MarkDuplicates I=./mapping/scDNA/$a/sorted.bam O=./mapping/scDNA/$a/dedup_sorted.bam M=ASXS.metrix REMOVE_DUPLICATES=true ASSUME_SORTED=true

  cd ./mapping/scDNA/$a

  samtools index dedup_sorted.bam
  split -n l/30 cells.txt

  for i in `ls x*`
  do
    sinto filterbarcodes -b dedup_sorted.bam -c $i -p 108 --barcode_regex '(?<=_)[^_]+(?=_)' --barcodetag "^V35" --outdir sinto_split
  done

  cd ~/OneDrop/OD_20231010
done
