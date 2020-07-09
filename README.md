# mosdepth_region_vis


### Run Mosdepth
```
mosdepth \
    --by 500 \
    -t 4 \
    -f human_g1k_v38_decoy_phix.fasta \
    sample-01 \
    sample-01.cram


mosdepth \
    --by 500 \
    -t 4 \
    -f human_g1k_v38_decoy_phix.fasta \
    sample-02 \
    sample-02.cram


mosdepth \
    --by 500 \
    -t 4 \
    -f human_g1k_v38_decoy_phix.fasta \
    sample-03 \
    sample-03.cram

```

### Tabix each sample file

```
tabix sample-01.per-base.bed.gz 
tabix sample-02.per-base.bed.gz 
tabix sample-03.per-base.bed.gz 

```

### Run region vis script
```
python plot_region_info.py \
    --gtf gencode.v34.annotation.gtf.gz \
    --gene-file gene_file.txt \
    --input sample-01.per-base.bed.gz sample-02.per-base.bed.gz sample-03.per-base.bed.gz \
    --combine \
    --sample-colors green blue orange
```
