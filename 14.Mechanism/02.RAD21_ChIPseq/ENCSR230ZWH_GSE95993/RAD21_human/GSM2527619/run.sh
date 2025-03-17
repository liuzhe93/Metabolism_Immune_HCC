fastq-dump SRR5331454.sra
fastqc --noextract SRR5331454.fastq
bowtie2 -p 6 -3 5 --local -x /mnt/zhangzheng_group/liuz-53/data/Human/Bowtie2_index/human_ref -U SRR5331454.fastq -S SRR5331454.sam
samtools view -bS SRR5331454.sam > SRR5331454.bam
samtools sort -@ 10 SRR5331454.bam -o SRR5331454.sorted.bam
samtools index SRR5331454.sorted.bam
macs2 callpeak -t SRR5331454.sorted.bam -q 0.05 -f BAM -g hs -n rep1
Rscript rep1_model.r
python Py_preoperate_peaks.py rep1_peaks.narrowPeak rep1_peaks.narrowPeak_clean
bamCoverage -p 4 -b SRR5331454.sorted.bam -o SRR5331454.bw



