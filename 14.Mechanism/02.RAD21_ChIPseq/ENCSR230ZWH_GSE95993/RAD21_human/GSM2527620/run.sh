fastq-dump SRR5331455.sra
fastqc --noextract SRR5331455.fastq
bowtie2 -p 6 -3 5 --local -x /mnt/zhangzheng_group/liuz-53/data/Human/Bowtie2_index/human_ref -U SRR5331455.fastq -S SRR5331455.sam
samtools view -bS SRR5331455.sam > SRR5331455.bam
samtools sort -@ 10 SRR5331455.bam -o SRR5331455.sorted.bam
samtools index SRR5331455.sorted.bam
macs2 callpeak -t SRR5331455.sorted.bam -q 0.05 -f BAM -g hs -n rep2
Rscript rep2_model.r
python Py_preoperate_peaks.py rep2_peaks.narrowPeak rep2_peaks.narrowPeak_clean
bamCoverage -p 4 -b SRR5331455.sorted.bam -o SRR5331455.bw

