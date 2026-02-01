tfbs=data/neuron_TFs/MA0080.4.SPI1_MA0080.4/beds/MA0080.4.SPI1_MA0080.4_all.bed
# Normalized bigwig files
signal1=data/celltypes/Astro/Astro_control_normalized_bigwig.bw
signal2=data/celltypes/DaN/DaN_control_normalized_bigwig.bw
signal3=data/celltypes/EC/EC_control_normalized_bigwig.bw
signal4=data/celltypes/EpC/EpC_control_normalized_bigwig.bw
signal5=data/celltypes/ExN/ExN_control_normalized_bigwig.bw
signal6=data/celltypes/FB/FB_control_normalized_bigwig.bw
signal7=data/celltypes/InN/InN_control_normalized_bigwig.bw
signal8=data/MG_control_normalized_bigwig_footprints.bw
signal9=data/celltypes/OPC/OPC_control_normalized_bigwig.bw
signal10=data/celltypes/Oligo/Oligo_control_normalized_bigwig.bw
signal11=data/celltypes/PC/PC_control_normalized_bigwig.bw
signal12=data/celltypes/TC/TC_control_normalized_bigwig.bw
regions=data/consensus_regions.bed

#TOBIAS FootprintScores --score footprint --flank-max 0 --window 200 --signal $signal1 $signal2 $signal3 $signal4 $signal5 $signal6 $signal7 $signal8 $signal9 $signal10 $signal11 $signal12 --regions $regions --output Astro_test_footprints.bw --cores 8
TOBIAS PlotAggregate --TFBS $tfbs --signals $signal1 $signal2 $signal3 $signal4 $signal5 $signal6 $signal7 $signal8 $signal9 $signal10 $signal11 $signal12 --output SPIC_test.pdf --share_y both --plot_boundaries --signal-on-x