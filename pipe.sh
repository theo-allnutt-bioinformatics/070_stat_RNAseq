#/bin/bash

#Cliff RNAseq pipeline
#trim reads
trimmo.py 'reads/*' 24 ~/db/adaptors.fa s y

#reference genome:
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Danio_rerio/reference/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_cds_from_genomic.fna.gz -o cds.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Danio_rerio/reference/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.fna.gz -o genome.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Danio_rerio/reference/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_protein.faa.gz -o protein.fa.gz

gunzip *.gz

#map to cds
bbmap_general.py 'clip/514_CC0U1ANXX_TTAGGC_L001_R1.fastq' cds.fa ./ 1
bbmap_general.py 'clip/515_CC0U1ANXX_TGACCA_L001_R1.fastq' cds.fa ./ 1
bbmap_general.py 'clip/516_CC0U1ANXX_ACAGTG_L001_R1.fastq' cds.fa ./ 1 
bbmap_general.py 'clip/524_CC0U1ANXX_GCCAAT_L001_R1.fastq' cds.fa ./ 1 
bbmap_general.py 'clip/525_CC0U1ANXX_CAGATC_L001_R1.fastq' cds.fa ./ 1 
bbmap_general.py 'clip/526_CC0U1ANXX_ACTTGA_L001_R1.fastq' cds.fa ./ 1 
bbmap_general.py 'clip/WT4_CC0U1ANXX_ATTCCT_L001_R1.fastq' cds.fa ./ 1 
bbmap_general.py 'clip/WT5_CC0U1ANXX_ATCACG_L001_R1.fastq' cds.fa ./ 1 
bbmap_general.py 'clip/WT6_CC0U1ANXX_CGATGT_L001_R1.fastq' cds.fa ./ 1

#tabulate results and trim
python ~/s/maps2tbl.py 'mapped/*.map' counts.txt

cut -f1,2,3,4,8,9,10 counts.tab > stat1.tab
cut -f1,5,6,7,8,9,10 counts.tab > stat2.tab
cut -f1,2,3,4,5,6,7 counts.tab > stat1v2.tab

R --file=stat.R