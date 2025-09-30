import pandas as pd
import os

PD_GWAS_df = pd.read_csv('/data/CARD_singlecell/SN_atlas/data/DLB_GWAS_scrape.tsv', delimiter='\t')
RSIDs = [x.split('-')[0].strip(' ') for x in PD_GWAS_df['STRONGEST SNP-RISK ALLELE']]

for RSID in RSIDs:
    template = f"curl -k -X GET 'https://ldlink.nih.gov/LDlinkRest/ldproxy?var={RSID}&pop=ALL&r2_d=r2&window=500000&genome_build=grch38&token=e748e680c922' -o /data/CARD_singlecell/SN_atlas/data/snp_LD_files/{RSID}.csv"
    os.system(template)
    print(RSID)