tfbs=BINDectect_output/MA1633.1.BACH1_MA1633.1/beds/MA1633.1.BACH1_MA1633.1_all.bed
bw_1=data/celltypes/Astro/Astro_control_normalized_bigwig.bw
bw_2=data/celltypes/Astro/Astro_PD_normalized_bigwig.bw
regions=data/celltypes/Astro/Astro_peaks.bed

TOBIAS FootprintScores --score footprint --flank-max 0 --window 200 --signal $bw_1 --regions $regions --output Astro_test_footprints.bw --cores 8
#TOBIAS PlotAggregate --TFBS $tfbs  --signals $bw_1 $bw_2 --output BACH1_test.pdf --share_y both --plot_boundaries --signal-on-x