#!/bin/bash
# version: 2

: '
分析時會遇到 genome 版本問題，或是要選擇物種，
因此將物種版本的選擇交給 scripts
但是 scripts 裡的路徑寫死就會很難更改，

因此另外用此 genomePATH.sh 將 genome 位置定義成變數，
直接 source 給環境，後續 scripts 像是 STAR.sh 就能引用這些變數。
'

export staridx_mm10=/staging/biology/ls807terra/0_genomes/star_index/mm10_star_2.7.9
export staridx_chm13=/staging/biology/ls807terra/0_genomes/star_index/CHM13_human
export staridx_hg38=/staging/biology/ls807terra/0_genomes/star_index/GRCh38_human/UCSC

export gtf_mm10=/staging/biology/ls807terra/0_genomes/genome_gtf/mm10/mm10.refGene.gtf
export gtf_chm13=/staging/biology/ls807terra/0_genomes/genome_gtf/CHM13/CHM13_v2.0.gtf
export gtf_hg38=/staging/biology/ls807terra/0_genomes/genome_gtf/hg38/UCSC_GRCh38_refGene.gtf

export bowtie2idx_mm10=/staging/biology/ls807terra/0_genomes/bowtie2_index/mm10/mm10
export bowtie2idx_chm13=/staging/biology/ls807terra/0_genomes/bowtie2_index/bowtie2_index_CHM13/CHM13
export bowtie2idx_hg38=/staging/biology/ls807terra/0_genomes/bowtie2_index/UCSC_GRCh38/GRCh38.p14

