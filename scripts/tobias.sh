
motifs=input/jaspar_2024_hsapiens.meme
signal1=data/celltypes/DaN/DaN_normalized_bigwig_test.bw
signal2=data/celltypes/ExN/ExN_control_normalized_bigwig.bw
signal3=data/celltypes/InN/InN_control_normalized_bigwig.bw
genome=/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa
peaks=data/consensus_regions.bed
out=data/neuron_TFs
cond1=DaN
cond2=ExN
cond3=InN


TOBIAS BINDetect --motifs $motifs --signals $signal1 --genome $genome --peaks $peaks --outdir $out --cond_names $cond1 --cores 32