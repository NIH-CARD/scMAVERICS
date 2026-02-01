# Define files 

# TF motifs
motifs=input/jaspar_2024_hsapiens.meme

# Normalized bigwig files
signal1=data/celltypes/Astro/Astro_control_normalized_bigwig.bw
signal2=data/celltypes/DaN/DaN_control_normalized_bigwig.bw
signal3=data/celltypes/EC/EC_control_normalized_bigwig.bw
signal4=data/celltypes/EpC/EpC_control_normalized_bigwig.bw
signal5=data/celltypes/ExN/ExN_control_normalized_bigwig.bw
signal6=data/celltypes/FB/FB_control_normalized_bigwig.bw
signal7=data/celltypes/InN/InN_control_normalized_bigwig.bw
signal8=data/celltypes/MG/MG_control_normalized_bigwig.bw
signal9=data/celltypes/OPC/OPC_control_normalized_bigwig.bw
signal10=data/celltypes/Oligo/Oligo_control_normalized_bigwig.bw
signal11=data/celltypes/PC/PC_control_normalized_bigwig.bw
signal12=data/celltypes/TC/TC_control_normalized_bigwig.bw

# Reference genome
genome=/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa
# Pan-cell type peaks
peaks=data/consensus_regions.bed
# Output folder
out=data/neuron_TFs

# Conditions to compare for enrichment
cond1=Astro
cond2=DaN
cond3=EC
cond4=EpC
cond5=ExN
cond6=FB
cond7=InN
cond8=MG
cond9=OPC
cond10=Oligo
cond11=PC
cond12=TC

# Command
TOBIAS BINDetect --motifs $motifs --signals $signal1 $signal2 $signal3 $signal4 $signal5 $signal6 $signal7 $signal8 $signal9 $signal10 $signal11 $signal12 --genome $genome --peaks $peaks --outdir $out --cond_names $cond1 $cond2 $cond3 $cond4 $cond5 $cond6 $cond7 $cond8 $cond9 $cond10 $cond11 $cond12 --cores 32