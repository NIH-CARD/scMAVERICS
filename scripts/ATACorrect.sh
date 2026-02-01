bed=data/celltypes/Astro/Astro_fragments.bed
genome=/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa
peaks=data/celltypes/Astro/Astro_peaks.bed

bedToBam -i $bed -g /fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2024-A/star/chrNameLength.txt > data/celltypes/Astro/Astro_fragments.bam
echo 'Done converting to Bam file'
TOBIAS ATACorrect --bam $bam --genome $genome --peaks $peaks --blacklist input/hg38-blacklist.bed --outdir ATACorrect_test --cores 32