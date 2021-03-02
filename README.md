# mosdepth_region_vis

[Mosdepth](https://github.com/brentp/mosdepth) is a fast method to get per-base coverage from sequencing data. However, 
Mosdepth does not provide an interactive approach to visualize the coverage data it creates. 

Visualizing coverage data for a set of genomic regions or genes of interest is a vital QC step in order to determine there is 
sufficient coverage support for each individual in a study. *mosdepth_region_vis* is an approach to take data coverage data 
from tools like mosdepth and provide a visual and interactive approach to quickly evaluate coverage information for region and
genes of interests. 


## SeqCover

Seqcover is a new approach to interactively visualize sequencing coverage at genes of interests. Seqcover replaces mosdepth_region_vis where
mosdepth_region_vis is no longer maintained.

Seqcover can be found here: https://github.com/brentp/seqcover



## How to use mosdepth_region_vis

Below are the common steps that should be taken to generate an interactive html file to evaluate coverage.

### Step 0: Software Requirements

Install the software requirements using the `requirements.txt` file. I recommend setting up a conda environment to house 
software requirements for this tool. 

Installing with conda:
```
conda install --file https://raw.githubusercontent.com/mikecormier/mosdepth_region_vis/master/requirements.txt 
```

### Step 1: Per-base Coverage

To get the per-base coverage information run `mosdepth` on each sample. 

Here is an example of running `modepth` for 3 samples:
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

Next, the per-base bed files need to be tabixed for quick interval lookup. 
Make sure that you tabix the *per-base* bed file. 
```
tabix sample-01.per-base.bed.gz 
tabix sample-02.per-base.bed.gz 
tabix sample-03.per-base.bed.gz 

```

### Step 2: Download Scripts

With the tabixed per-base files, use the `plot_regions_info.py` script and the `tmpl.html` file to generate an interactive html file.

First you need to download the `plot_regions_info.py` script and the `tmpl.html` file from this repository. 
```
wget https://raw.githubusercontent.com/mikecormier/mosdepth_region_vis/master/plot_region_info.py

wget https://raw.githubusercontent.com/mikecormier/mosdepth_region_vis/master/tmpl.htm
```

The `plot_regions_info.py` file will use the `tmpl.html` to generate the final html file.

The input parameters include: 
```
$ python plot_region_info.py -h

    usage: plot_region_info.py [-h] --gtf GTF File
                               [--region Genomic Region [Genomic Region ...]]
                               [--region-file File of Genomic Region]
                               [--gene Gene Symbol [Gene Symbol ...]]
                               [--gene-file File of Gene Symbols] [-o OUTPUT]
                               [--tmpl Template HTML] [--combine]
                               [--lch-cutoff Low Coverage Highlight Cutoff]
                               [--sample-colors Sample Colors [Sample Colors ...]]
                               --input Input Coverage Files)
                               [Input Coverage File(s ...]

    Creates html plots from mosdepth results based on user defined region(s)

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT, --output OUTPUT
                            Path and/or name of output file. Directories must
                            exist. (Default = 'region_coverage.html')
      --tmpl Template HTML  Path and/or name of the template html file. This is
                            the template html file that is distributed with the
                            script. (Default = 'tmpl.html')
      --combine             True or False, whether or not to combine all regions
                            into a single html file or not. If '--combine' is
                            added, all regions will be combined into a single html
                            file. If '--combine' is not added, each region will be
                            written to a separate html file. (NOTE: when --combine
                            is set, the size of the html file will increase.
                            Depending on the number of regions and the size of
                            each region, the html file may be to large to load)
      --lch-cutoff Low Coverage Highlight Cutoff
                            A coverage value cutoff, where any coverage at or
                            bellow the cutoff will be highlighted in per base
                            region plot. (Default = 10)
      --sample-colors Sample Colors [Sample Colors ...]
                            A space separated list of colors to use for the
                            samples while plotting. If --sample-colors is not used
                            or the number of colors provided by --sample-colors
                            does not match the number of samples then colors will
                            be chosen at random. (Example --sample-colors green
                            blue orange) (Default = random color per sample)

    Required Arguments:
      --gtf GTF File        (Required) Path to a gtf file with gene features
      --input Input Coverage File(s) [Input Coverage File(s) ...]
                            One ore more coverage bed files from mosdepth to get
                            coverage info for. (3 sample example: --input
                            sample1.per-base.bed.gz sample2.per-base.bed.gz
                            sample3.per-base.bed.gz)

    Required Argument - ONLY ONE:
      --region Genomic Region [Genomic Region ...]
                            A tabix styled region. (Example: --region
                            chr11:1234567-1234578) Multiple regions can be added,
                            separated by a space. (Example: --region
                            chr11:1234567-1234578 chr1:987654321:987655432).
                            Either --region, --region-file, --gene, or --gene-file
                            is required.
      --region-file File of Genomic Region
                            A file of tabix styled regions. One region per line.
                            Either --region, --region-file, --gene, or --gene-file
                            is required.
      --gene Gene Symbol [Gene Symbol ...]
                            A Gene symbol to get info for. (Example: --gene
                            APEX1).Multiple genes can be added, separated by a
                            space. (Example: --gene APEX1 ABCA1) Either --region,
                            --region-file, --gene, or --gene-file is required.
      --gene-file File of Gene Symbols
                            A file of gene symbols to get info for. One gene
                            symbol per line. Either --region, --region-file,
                            --gene, or --gene-file is required.

```

NOTE: if you use the `--combine` parameters the output file will be quite large. Without the `--combine` parameter a separate html file will 
be generated for each region or gene.


### Step 3: Generate interactive html file

The `plot_region_info.py` script generates a few plots for each gene/region. Every sample included in the parameters will be plotted for each gene.

The script will generate the following figures for each gene:

    - A table of descriptive statistic for each sample
    - A Kernel Density Plot histogram of coverage for each sample
    - A z-score adjusted Kernel Density Plot histogram of coverage for each sample
    - A plot for each sample showing the proportion of bases for a region covered at a certain coverage cutoff
    - Sample vs Sample per base coverage value plot. Used to detect any outlier coverage 
    - A region/gene plot with the per-base across the region. A gene track is included along with a per-base coverage track for each sample. 


Each plot is interactive with hover information, zooming, scaling, and other features. 

The most popular plot is the gene + coverage track plot. This shows where within the genic region the coverage changes. Any bases that fall 
below a user defined cutoff will be highlighted. Zooming along the genomic position is coordinate across the different tracks. 

NOTE:
    The `tmpl.html` file needs to be in the same directory where the `plot_region_info.py` file is. If it is not then the script won't be able to find it.

An example of generate a combined html file:
```
python plot_region_info.py \
    --gtf <gtf-file> \
    --gene-file <gene_file.txt> \
    --input sample-01.per-base.bed.gz sample-02.per-base.bed.gz sample-03.per-base.bed.gz \
    --combine \
    --sample-colors green blue orange
```

NOTE: 
    
    - The <gtf-file> variable above represents a gtf-file to use for annotations. (This needs to be the same genome build used to align the files)
    - The <gene-file> variable above represents a file of gene names to look at.  



### Step 4: Interactive Visualization 

With the new file created you can now start exploring the per sample coverage across different genes or regions of interest. 

Here is an example of the html output. This output includes 3 random samples from the 1000G project and a few random genes. 



