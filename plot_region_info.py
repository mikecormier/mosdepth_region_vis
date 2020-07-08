import sys
import os
import string
import scipy
import copy
from scipy import stats
import collections
from argparse import ArgumentParser
import pyranges as pr
import pandas as pd
import numpy as np
import pysam 
import random 
import matplotlib
import matplotlib.colors as pltc
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots


#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def get_args():
    parser = ArgumentParser(description="Creates html plots from mosdepth results based on user defined region(s)")

    req = parser.add_argument_group("Required Arguments")
    one_req = parser.add_argument_group("Required Argument - ONLY ONE")

    req.add_argument(
        "--gtf",
        metavar = "GTF File",
        required = True,
        help = "(Required) Path to a gtf file with gene features"
    )

    one_req.add_argument(
        "--region",
        metavar = "Genomic Region",
        nargs='+',
        default = None,
        help = "A tabix styled region. (Example: --region chr11:1234567-1234578) Multiple regions can be added, seperated by a space. (Example: --region chr11:1234567-1234578 chr1:987654321:987655432). Either --region, --region-file, --gene, or --gene-file is required."
    )

    one_req.add_argument(
        "--region-file",
        metavar = "File of Genomic Region",
        default = None,
        help = "A file of tabix styled regions. One region per line. Either --region, --region-file, --gene, or --gene-file is required."
    )

    one_req.add_argument(
        "--gene",
        metavar = "Gene Symbol",
        nargs='+',
        default = None,
        help = "A Gene sybmol to get info for. (Example: --gene APEX1).Multiple genes can be added, sperated by a space. (Example: --gene APEX1 ABCA1) Either --region, --region-file, --gene, or --gene-file is required."
    )

    one_req.add_argument(
        "--gene-file",
        metavar = "File of Gene Symbols",
        default = None,
        help = "A file of gene symbols to get info for. One gene sybmol per line. Either --region, --region-file, --gene, or --gene-file is required."
    )

    parser.add_argument(
        "-o", 
        "--output",
        default="region_coverage.html",
        help="Path and/or name of output file. Directories must exist. (Default = 'region_coverage.html')"
    )

    parser.add_argument(
        "--tmpl",
        metavar = "Template HTML",
        default = "tmpl.html",
        help="Path and/or name of the template html file. This is the template html file that is distributed with the script. (Default = 'tmpl.html')"
    )

    parser.add_argument(
        "--combine",
        action="store_true",
        help="True or False, whether or not to combine all regions into a single html file or not. If '--combine' is added, all regions will be combined into a single html file. If '--combine' is not added, each region will be written to a seperate html file. (NOTE: when --combine is set, the size of the html file will increase. Depending on the number of regions and the size of each region, the html file may be to large to load)"
    )

    parser.add_argument(
        "--lch-cutoff",
        metavar = "Low Coverage Highlight Cutoff",
        default = 10,
        help = "A coverage value cutoff, where any coverage at or bellow the cutoff will be highlighted in per base region plot. (Default = 10)"
    )

    parser.add_argument(
        "--sample-colors",
        metavar = "Sample Colors",
        nargs='+',
        help = "A space seperated list of colors to use for the samples while plotting. If --sample-colors is not used or the number of colors provided by --sample-colors does not match the number of samples then colors will be chosen at random. (Example --sample-colors green blue orange) (Default = random color per sample)"
    )

    req.add_argument(
        "--input",
        metavar = "Input Coverage File(s)",
        nargs='+',
        required = True,
        help="One ore more coverage bed files from mosdepth to get coverage info for. (3 sample example: --input sample1.per-base.bed.gz sample2.per-base.bed.gz sample3.per-base.bed.gz)"
    )

    return parser.parse_args()


#---------------------------------------------------------------------------------------------------------------------------------
## Functions/Methods
#---------------------------------------------------------------------------------------------------------------------------------

def get_dict_from_region_string(region_string):
    """
    get_dict_from_region_string
    ===========================
    Get a diciontary from a tabix styled region string. Keys will include "chrom", "start", and "end"

    Parameters: 
    -----------
    1) regoin_string: (str) A tabix styled region string. (Example: chr14:123456-124456) 

    Returns:
    ++++++++
    1) (dict) A dictionary with keys: 'chrom', 'start', 'end' and genomic region info as values
    """
    
    region_list = region_string.strip().replace(",","").replace("-",":").split(":")
    region_dict = {"chrom":str(region_list[0]), "start":int(region_list[1]), "end":int(region_list[2])}

    return(region_dict)



def get_sample_region_coverage(sample_tabix_object, 
                               region_dict):
    """
    get_sample_region_coverage
    ==========================
    Get the coverage values for a specific sample at a specific region interval. 

    A list containing per base coverage values will be returend. 

    Parameters:
    -----------
    1) sample_tabix_object: (psyam Object) A pysam tabix object representing the bed file of a sample 
    2) region_dict: (dict) A dictionary with region information. Required keys include "chrom", "start", and "end"

    Returns:
    +++++++
    1) (list) A list of coverage values representing a per base coverage score
    """

    coverage = []
    ## Iterate over the overlaping coverate intervals
    for i,coverage_line in enumerate(sample_tabix_object.fetch(region_dict["chrom"], int(region_dict["start"]), int(region_dict["end"]))):
    
        ## Create a dict for the current line 
        cov_dict = dict(zip(["chrom","start","end","coverage_value"],coverage_line.strip().split("\t")))

        ## Add the coverage value for each position defined by the current interval 
        coverage.extend([int(cov_dict["coverage_value"]) for x in range(int(cov_dict["end"]) - int(cov_dict["start"]))])

    return(coverage)


def plot_kde_per_sample(per_base_coverage_by_sample, 
                        sample_labels, 
                        sample_color_dict, 
                        region_dict):
    """
    plot_kde_per_sample
    ===================
    Plot the kernal density estimator plot for each sample. Use values from a histrogram to generate the kernal desnity estimator predictions,
     and plot the density predictions. 

    Get the distribution of coverage per sample at a specific region.

    Parameters:
    -----------
    1) per_base_coverage_by_sample: (2d list) A 2d list, with each list representing the per base coverage values for a specific sample for the current query region 
    2) sample_labels: (list) A list of sample ids/labels, with each index associated with the index in per_base_coverage_by_sample. (That is, index 1 in 
                              per_base_coverage_by_sample should represent the coverage for the sample at index 1 in sample_labels)
    3) sample_color_dict: (dict) A dictionary with keys as sample ids/labels and value as a color. (Used for plotting colors)
    4) region_dict: (dict) A dictionary with region finromation. Required keys include "chrom", "start", and "end". (See the 'get_dict_from_region_string' function)

    Returns:
    +++++++
    1) (str) A plotly embedded html string for the per sample kde plot.
    """
    
    fig = ff.create_distplot(per_base_coverage_by_sample, 
                             sample_labels, 
                             #bin_size = .1, 
                             #curve_type = "normal",
                             colors = [sample_color_dict[x] for x in sample_labels],
                             show_rug = False,
                             show_hist = False,
                             )

    fig.update_layout(title_text = "Density Coverage Distribution for {}:{}-{}".format(region_dict["chrom"], region_dict["start"], region_dict["end"]))
    fig.update_layout(legend_title_text = "Samples")
    fig.update_layout(xaxis_title = "Coverage")
    fig.update_layout(yaxis_title = "Density")
    return(fig.to_html(full_html = False, include_plotlyjs = "cdn"))


def plot_z_score_distribution(per_base_coverage_by_sample, 
                              sample_labels, 
                              sample_color_dict, 
                              region_dict):
    """
    plot_z_score_distriubtion 
    =========================
    Plot the distribution of per base coverage scores by sample converted to z-scores. The distribution will be converted to desnity scores using a kernal density estimator. 

    The z-scores are relative to each samples mean. That is, each sample's z-scores are determined by that samples coverage distrubtion. 

    Parameters:
    -----------
    1) per_base_coverage_by_sample: (2d list) A 2d list, with each list representing the per base coverage values for a specific sample for the current query region 
    2) sample_labels: (list) A list of sample ids/labels, with each index associated with the index in per_base_coverage_by_sample. (That is, index 1 in 
                              per_base_coverage_by_sample should represent the coverage for the sample at index 1 in sample_labels)
    3) sample_color_dict: (dict) A dictionary with keys as sample ids/labels and value as a color. (Used for plotting colors)
    4) region_dict: (dict) A dictionary with region finromation. Required keys include "chrom", "start", and "end". (See the 'get_dict_from_region_string' function)

    Returns:
    +++++++
    1) (str) A plotly embedded html string for the per sample z-score kde plot.
    2) (2d list) A 2d list representing the relative coverage corrected z-scores per sample
    """

    ## get z scores
    z_scores = [stats.zscore(x) for x in per_base_coverage_by_sample]
    
    fig = ff.create_distplot(z_scores, 
                             sample_labels, 
                             colors = [sample_color_dict[x] for x in sample_labels],
                             show_rug = False,
                             show_hist = False,
                             )

    fig.update_layout(title_text = "Density Z-Score Coverage Distribution for {}:{}-{}".format(region_dict["chrom"], region_dict["start"], region_dict["end"]))
    fig.update_layout(legend_title_text = "Samples")
    fig.update_layout(xaxis_title = "Relative Z-Score for Sample Coverage")
    fig.update_layout(yaxis_title = "Density")
    return(fig.to_html(full_html = False, include_plotlyjs = "cdn"), z_scores)


def get_plotly_table(per_base_coverage_by_sample, 
                     sample_labels,
                     low_coverage_cutoff):
    """
    get_plotly_table 
    ================
    Get a table that provides descriptive statistics on coverage for each sample at a specific region.

    Parameters:
    -----------
    1) per_base_coverage_by_sample: (2d list) A 2d list, with each list representing the per base coverage values for a specific sample for the current query region 
    2) sample_labels: (list) A list of sample ids/labels, with each index associated with the index in per_base_coverage_by_sample. (That is, index 1 in 
                              per_base_coverage_by_sample should represent the coverage for the sample at index 1 in sample_labels)
    3) low_coverage_cutoff: (int) An int that is used to identify positions at or below low coverage 

    Returns:
    +++++++
    1) (str) A plotly embedded html string for per sample descriptive statistic table.
    2) (int) The max converage value across all samples. 
    3) (bool) True of False, whether or not any of the samples has a position at or below the low coverage cutoff
    """
    
    header = ["<b>Samples</b>","<b>Mean</b>","<b>Median</b>","<b>SD</b>","<b>Min</b>","<b>Max</b>","<b>Range</b>","<b>Low Coverage Bases</b><br>   coverage <= {}".format(low_coverage_cutoff)]

    samples_cell = []
    mean_cell = []
    median_cell = []
    min_cell = []
    max_cell = []
    sd_cell = []
    range_cell = []
    low_coverage_cell = []

    ## get descriptive statistics
    for i,cov_list in enumerate(per_base_coverage_by_sample):
        
        samples_cell.append("<b>{}</b>".format(sample_labels[i]))
        mean_cell.append("{:.3f}".format(np.mean(cov_list)))
        median_cell.append("{:.3f}".format(np.median(cov_list)))
        min_cell.append(np.amin(cov_list))
        max_cell.append(np.amax(cov_list))
        sd_cell.append("{:.3f}".format(np.std(cov_list)))
        range_cell.append(np.ptp(cov_list))
        low_coverage_cell.append((np.array(cov_list) <= low_coverage_cutoff).sum())
        

    fig = go.Figure(data = [go.Table(
        columnwidth = [100,70,80,60,60,60,70,200],
        header = dict(values = header,
                      line_color = "darkslategray",
                      fill_color = "royalblue",
                      align=["left","center"],
                      font = dict(color = "white", size = 15),
                      height = 40
                      ),
        cells = dict(values = [samples_cell, mean_cell, median_cell, sd_cell, min_cell, max_cell, range_cell, low_coverage_cell],
                     line_color = "darkslategray",
                     fill= dict(color = ["royalblue","white","white","white","white","white","white","white"]),
                     align=["left","center"],
                     font = dict(color = ["white","black","black","black","black","black","black", "black"],
                                 size = [15,12,12,12,12,12,12,12]
                                ),
                    height = 30
                    )
        )])
                     
    fig.update_layout(title_text = "Descriptive Statistics Table")

    return(fig.to_html(full_html = False, include_plotlyjs = "cdn"), max(max_cell), bool( max(low_coverage_cell) > 0))


def plot_sample_vs_sample_per_base_coverage(per_base_coverage_by_sample, 
                                            sample_labels, 
                                            sample_color_dict, 
                                            region_dict):
    """
    plot_sample_vs_sample_per_base_covearge 
    =======================================
    Plot each sample's coverage vs anothe sample's coverage. Currenlty, this is limited to a max of 3 samples

    Parameters:
    -----------
    1) per_base_coverage_by_sample: (2d list) A 2d list, with each list representing the per base coverage values for a specific sample for the current query region 
    2) sample_labels: (list) A list of sample ids/labels, with each index associated with the index in per_base_coverage_by_sample. (That is, index 1 in 
                              per_base_coverage_by_sample should represent the coverage for the sample at index 1 in sample_labels)
    3) sample_color_dict: (dict) A dictionary with keys as sample ids/labels and value as a color. (Used for plotting colors)
    4) region_dict: (dict) A dictionary with region finromation. Required keys include "chrom", "start", and "end". (See the 'get_dict_from_region_string' function)

    Returns:
    +++++++
    1) (str) A plotly embedded html string for the sample vs sample coverate plots.
    """

    if len(per_base_coverage_by_sample) > 3: 
        print("\n!!WARNING!! There are more then 3 samples. Unable to plot sample vs sample")
        return()

    fig = make_subplots(rows=3)

    row_index = 0

    ## Subplot of per base coverage
    for i,coverage_list1 in enumerate(per_base_coverage_by_sample):
        for j,coverage_list2 in enumerate(per_base_coverage_by_sample):

            ## Skip if sample info already used
            if j <= i:
                continue
            row_index += 1

            ## Get a set of uniqe x,y pairs
            unique_pairs = set()
            for cov1, cov2 in zip(coverage_list1, coverage_list2):
                unique_pairs.add((cov1,cov2))

            ## Seperate the pairs into x values and y values
            x_values = []
            y_values = []   
            for pair in unique_pairs:
                x_values.append(pair[0])
                y_values.append(pair[1])
                

            fig.add_trace(go.Scatter(x = x_values, 
                                     y = y_values, 
                                     showlegend = False,
                                     mode = "markers"),
                          row = row_index, 
                          col = 1)

            fig.update_xaxes(title_text = "{} Coverage".format(sample_labels[i]),row=row_index, col=1)
            fig.update_yaxes(title_text = "{} Coverage".format(sample_labels[j]),row=row_index, col=1)

    fig.update_layout(title_text = "Sample vs Sample per Base Coverage for region {}:{}-{}".format(region_dict["chrom"], region_dict["start"], region_dict["end"]))

    return(fig.to_html(full_html = False, include_plotlyjs = "cdn"))
        

def plot_porportion_coverage(per_base_coverage_by_sample, 
                             sample_labels, 
                             sample_color_dict, 
                             region_dict):
    """
    plot_proportion_coverage 
    ========================
    Plot the proportion of bases covered in a region per each coverage cutoff for each sample. 

    Parameters:
    -----------
    1) per_base_coverage_by_sample: (2d list) A 2d list, with each list representing the per base coverage values for a specific sample for the current query region 
    2) sample_labels: (list) A list of sample ids/labels, with each index associated with the index in per_base_coverage_by_sample. (That is, index 1 in 
                              per_base_coverage_by_sample should represent the coverage for the sample at index 1 in sample_labels)
    3) sample_color_dict: (dict) A dictionary with keys as sample ids/labels and value as a color. (Used for plotting colors)
    4) region_dict: (dict) A dictionary with region finromation. Required keys include "chrom", "start", and "end". (See the 'get_dict_from_region_string' function)

    Returns:
    +++++++
    1) (str) A plotly embedded html string for the proportion covearge plot.
    """

    fig = go.Figure()

    for i,coverage_list in enumerate(per_base_coverage_by_sample):

        x_values, y_values = get_region_coverage_proportion(coverage_list)


        fig.add_trace(go.Scatter(x = x_values, 
                                 y = y_values, 
                                 mode = "lines", 
                                 name = sample_labels[i],
                                 line = {"color":sample_color_dict[sample_labels[i]]})) 

    fig.update_layout(title_text = "Proportion of bases covered at coverage cutoff for region {}:{}-{}".format(region_dict["chrom"], region_dict["start"], region_dict["end"]))
    fig.update_xaxes(title_text = "Coverage")
    fig.update_yaxes(title_text = "Proportion of bases")
    fig.update_layout(legend_title_text = "Samples")
    return(fig.to_html(full_html = False, include_plotlyjs = "cdn"))


def get_region_coverage_proportion(coverage_list):
    """
    get_region_coverage_proportion 
    ===============================
    Method to get the proportion of bases covered for each coverage cutoff

    Parameters:
    -----------
    1) coverage_list: (list) A list of coverage values 

    Returns:
    +++++++
    1) x_values: (list) A list of x axis values. (Coverage cutoffs)
    2) y_values: (list) A list of y axis values. (Proportion of bases)
    NOTE: indicies between the two list represent an x,y pair
    """
    
    ## Get coverage count by coverage value
    coverage_dict = collections.defaultdict(int)
    for cov in coverage_list:
        coverage_dict[cov] += 1
    
    ## Convert to df
    rows = []
    for key in sorted(coverage_dict.keys()):
        rows.append([key,coverage_dict[key]])
    df = pd.DataFrame(rows, columns = ["cov_value","cov_count"])
    
    ## Get proportion of bases covered per coverage value
    x_values = []
    y_values = []
    for key in sorted(coverage_dict.keys(), reverse = True):
        
        ## Get proportion of bases per cov value
        ## Get the total number of bases equal to or above the coverage value divided by the total number of bases 
        proportion = ( df.loc[df.cov_value >= key].cov_count.sum() ) / df.cov_count.sum() 
        x_values.append(key)
        y_values.append(proportion)

    return(x_values, y_values)


def get_plotly_shapes_from_gtf_dict(gtf_info_dict, 
                                    xref = "x3", 
                                    yref = "y3"):
    """
    get_plotly_shapes_from_gtf_dict
    ===============================
    Method to create plotly shapes for features in a gtf file. These shapes will be used to plot a "gene track" for a 
     region of interests. 

    Current features used include:
        transcripts 
        exons

    Parameters:
    ----------
    1) gtf_info_dict: (dict) A dictionary of gene feature info for the current region of interests. The dictionary should be formated like the output from 
                              the 'get_gtf_region_position_info' function. (It is best to use the 'get_gtf_region_position_info' to get the correct dict) 
    2) xref: (str) A string represeting the plotly subplot column x position to add the output to
    2) yref: (str) A string represeting the plotly subplot column y position to add the output to

    Returns:
    ++++++++
    1) shapes: (list of dicts) A list of shape dictionaries, with a shape reprenting a transcript or an exon 
    2) exon_position_df: (pandas DataFrame) A pandas DataFrame representing the positions of each exon in the "shapes" object
    3) text_position_df: (pandas DataFrame) A pandas DataFrame representing transcript text, (strand), for each transcript in the "shapes" object
    4) y_tick_labels: (list) A list of y_tick_labels to use for a gene track. The tick lables represent transcript ids 
    5) y_tick_ref_vals: (list) A list of y tick reference values, or coordiantes, associated with each y_tick_label. The combination of these 
                                coordinates with the tick labels allows for each transcript feature to be plotted with a transcript id as a tick label.
    NOTE: If no gene feature info available, None is returned 
    """
    ## Check for an empty dict
    if not gtf_info_dict:
        return(None, None, None, None, None)

    start = 10
    end = 9
    offset = 0.25
    shapes = []
    exon_point_header = ["x_position","y_position","gene_id","transcript_id","gene_symbol","strand","gene_biotype","exon_start","exon_end","exon_number"]
    exon_point_info = []

    text_point_header = ["x_position","y_position","text"]
    text_point_info = []

    y_tick_labels = []
    y_tick_ref_vals = []

    for gene_id, transcript_dict in gtf_info_dict.items():
        
        ## Get transcript line shape info
        for transcript_id, transcript_info in transcript_dict.items():
            shapes.append({"type": "line", 
                           "x0": transcript_info["start"], 
                           "x1": transcript_info["end"], 
                           "y0": end, 
                           "y1": end, 
                           "line_color":"black", 
                           "xref":xref, 
                           "yref":yref}) 

            y_tick_labels.append(transcript_id)
            y_tick_ref_vals.append(end)

            ## Add strand text
            for pos in (np.arange(transcript_info["end"] - transcript_info["start"]) + transcript_info["start"]):
                text_point_info.append([pos, end + 0.05, ">" if transcript_info["strand"] == "+" else "<" if transcript_info["strand"] == "-" else "" ])

            ## Get Exon rectangle shape info
            ## Get exon scatter point info
            for exon_dict in transcript_info["exons"]:
                shapes.append({"type":"rect", 
                               "x0": exon_dict["start"], 
                               "x1": exon_dict["end"], 
                               "y0": end - offset, 
                               "y1": end + offset, 
                               "line_color": "black",
                               "fillcolor":"black",
                               "xref":xref, 
                               "yref":yref})
                exon_point_info.append([(((exon_dict["end"] - exon_dict["start"]) / 2) + exon_dict["start"]), 
                                        end, 
                                        gene_id, 
                                        transcript_id, 
                                        transcript_info["gene_symbol"],
                                        transcript_info["strand"],
                                        transcript_info["biotype"],
                                        exon_dict["start"],
                                        exon_dict["end"],
                                        exon_dict["exon_number"]
                                      ])
            end -= 1


    exon_position_df = pd.DataFrame(exon_point_info, columns = exon_point_header)
    text_position_df = pd.DataFrame(text_point_info, columns = text_point_header)

    return(shapes, exon_position_df, text_position_df, y_tick_labels, y_tick_ref_vals)


def get_shapes_for_low_coverage_cutoff(per_base_coverage_by_sample,
                                       region_position_list,
                                       low_coverage_cutoff,
                                       max_coverage_value,
                                       max_z_score,
                                       min_z_score,
                                       has_low_coverage):

    """
    get_shapes_for_low_coverage_cutoff
    ==================================
    Method to make rectangle plotly shapes that represent regions of low coverage based on a low coverage cutoff.

    Any single or consecutive regions of low coverage will be assigned a low coverage shape.

    Parameters:
    -----------
    1) per_base_coverage_by_sample: (2d list) A 2d list, with each list representing the per base coverage values for a specific sample for the current query region 
    2) region_position_list: (list) A list of positions that correlate with the per_base_coveage for each sample
    3) low_coverage_cutoff: (int) An int representing the coverage at or below that should be flagged as "low coverage"
    4) max_coverage_value: (int) The max coverage value between all samples for the current region 
    5) max_z_score: (int) The max z-score coverage value between all samples for the current region
    6) min_z_score: (int) The min z-score coverage value between all samples for the current region
    7) has_low_coverage: (bool) True or False, whether or not at least on sample has a position of low coverage based on the low coverage cutoff

    Returns:
    ++++++++
    1) shapes: (list of dicts) A list of dictionaries, where each dict represents a low coverage shape to add to the plot.
    """

    ## Identify low coverage areas
    shapes = []
    if has_low_coverage:
        for cov_list in per_base_coverage_by_sample:
            prev_pos = -1
            low_cov_start = -1
            for cov,pos in zip(cov_list,region_position_list):

                ## If low covearge 
                if cov <= low_coverage_cutoff:
                    ## If consecutive positions
                    if prev_pos + 1 == pos:
                        prev_pos = pos
                    ## If new region
                    else:
                        if prev_pos != -1:
                            ## Add shapes
                            shapes.append({"type":"rect",
                                           "x0": low_cov_start - 0.1,
                                           "x1": prev_pos + 0.1 if prev_pos != low_cov_start else low_cov_start + 0.1,
                                           "y0": 0,
                                           "y1": max_coverage_value,
                                           "line_color": "grey",
                                           "fillcolor":"red",
                                           "opacity":0.5,
                                           "layer":"below",
                                           "xref":"x1",
                                           "yref":"y1",
                                           "visible":True})
                                
                            shapes.append({"type":"rect",
                                           "x0": low_cov_start - 0.1,
                                           "x1": prev_pos + 0.1 if prev_pos != low_cov_start else low_cov_start + 0.1,
                                           "y0": min_z_score,
                                           "y1": max_z_score,
                                           "line_color": "grey",
                                           "fillcolor":"red",
                                           "opacity":0.5,
                                           "layer":"below",
                                           "xref":"x2",
                                           "yref":"y2",
                                           "visible":True})

                        ## set position trackers to current pos
                        low_cov_start = pos
                        prev_pos = pos

            ## Add shapes
            if prev_pos != -1:
                shapes.append({"type":"rect",
                               "x0": low_cov_start - 0.1,
                               "x1": prev_pos + 0.1 if prev_pos != low_cov_start else low_cov_start + 0.1,
                               "y0": 0,
                               "y1": max_coverage_value,
                               "line_color": "grey",
                               "fillcolor":"red",
                               "opacity":0.5,
                               "layer":"below",
                               "xref":"x1",
                               "yref":"y1",
                               "visible":True})
                    
                shapes.append({"type":"rect",
                               "x0": low_cov_start - 0.1,
                               "x1": prev_pos + 0.1 if prev_pos != low_cov_start else low_cov_start + 0.1,
                               "y0": min_z_score,
                               "y1": max_z_score,
                               "line_color": "grey",
                               "fillcolor":"red",
                               "opacity":0.5,
                               "layer":"below",
                               "xref":"x2",
                               "yref":"y2",
                               "visible":True})

    return(shapes)
                

def plot_per_base_coverage(per_base_coverage_by_sample, 
                           sample_labels, 
                           sample_color_dict, 
                           region_dict, 
                           gtf_info_dict, 
                           coverage_highlight_cutoff,
                           max_coverage_value,
                           has_low_coverage):
    """
    plot_per_base_coverage
    ======================
    plot the coverage profile for each sample at each base for the current query region.

    For each sample, subplots include:
        - Coverage profile at each base for the current region
        - Z-Score converted coverage profile at each base for the current region 
        - Gene track for the current region 

    The size of the coverage profiles are reduced based on consecutive bases with the same coverage value. 

    Parameters:
    -----------
    1) per_base_coverage_by_sample: (2d list) A 2d list, with each list representing the per base coverage values for a specific sample for the current query region 
    2) sample_labels: (list) A list of sample ids/labels, with each index associated with the index in per_base_coverage_by_sample. (That is, index 1 in 
                              per_base_coverage_by_sample should represent the coverage for the sample at index 1 in sample_labels)
    3) sample_color_dict: (dict) A dictionary with keys as sample ids/labels and value as a color. (Used for plotting colors)
    4) region_dict: (dict) A dictionary with region finromation. Required keys include "chrom", "start", and "end". (See the 'get_dict_from_region_string' function)
    5) gtf_info_dict: (dict) A dictionary of gene features from a gtf file for the current region generated using the 'get_gtf_region_position_info' function. 
    6) coverage_highlight_cutoff: (int) A int representing coverage at or below that should be highlighted  
    7) max_coverage_value: (int) The max coverage value between the samples
    8) has_low_coverage (bool) True or False, whether or not at least on sample has coverage at or below the coverage cutoff

    Returns:
    +++++++
    1) (str) A plotly embedded html string for the proportion covearge plot.
    2) (int) The number of low coverage shapes added to the plot
    """
    ## Get a list of positions 
    region_position_list = np.arange(region_dict["end"] - region_dict["start"]) + region_dict["start"]
    
    ## Reduce the number of bases looked at by combining bases with same coverage score
    compressed_per_base_coverage_by_sample = [] 
    compressed_region_position_list = []
    for cov_list in per_base_coverage_by_sample:
        new_cov_list = []
        new_pos_list = []
        prev_cov = -1
        for base_cov, pos, in zip(cov_list, region_position_list):
            
            ## if coverage differs
            if base_cov != prev_cov:
                new_cov_list.append(base_cov)
                new_pos_list.append(pos)
                prev_cov = base_cov 
        
        ## If the last position is not added, add it 
        if new_pos_list[-1] != region_position_list[-1]:
            new_pos_list.append(region_position_list[-1])
            new_cov_list.append(cov_list[-1])

        compressed_per_base_coverage_by_sample.append(new_cov_list)
        compressed_region_position_list.append(new_pos_list)

    ## get z scores
    z_scores = [stats.zscore(x) for x in compressed_per_base_coverage_by_sample]

    ## get max and min z scores
    max_z_score = 0
    min_z_score = 0
    for z_score_list in z_scores:
        sample_max_z_score = np.amax(z_score_list)
        sample_min_z_score = np.amin(z_score_list)
        if sample_max_z_score > max_z_score:
            max_z_score = sample_max_z_score
        if sample_min_z_score < min_z_score:
            min_z_score = sample_min_z_score

    ## get gene track shapes
    shapes, exon_pos_df, text_pos_df, y_tick_labels, y_tick_ref_vals = get_plotly_shapes_from_gtf_dict(gtf_info_dict)

    ## get low coverage shapes
    low_cov_shapes = get_shapes_for_low_coverage_cutoff(per_base_coverage_by_sample, 
                                                        region_position_list, 
                                                        coverage_highlight_cutoff,
                                                        max_coverage_value,
                                                        max_z_score,
                                                        min_z_score,
                                                        has_low_coverage)

    ## Subplots
    if shapes is not None:
        fig = make_subplots(rows=3, shared_xaxes=True, vertical_spacing=0.005)
        max_rows = 3

        ## get final list of shapes
        final_shapes = low_cov_shapes + shapes

    else:
        fig = make_subplots(rows=2, shared_xaxes=True, vertical_spacing=0.005)
        max_rows = 2

        ## get final list of shapes
        final_shapes = low_cov_shapes 

    for i,coverage_list in enumerate(compressed_per_base_coverage_by_sample):
        fig.add_trace(go.Scatter(x = compressed_region_position_list[i], 
                                 y = coverage_list, 
                                 mode = "lines", 
                                 legendgroup = "samples{}".format(sample_labels[i]),
                                 name = sample_labels[i],
                                 line = {"color":sample_color_dict[sample_labels[i]]}), 
                      row = 1, 
                      col = 1)

    for i,z_score_list in enumerate(z_scores): 
        fig.add_trace(go.Scatter(x = compressed_region_position_list[i], 
                                 y = z_score_list, 
                                 mode = "lines", 
                                 legendgroup = "samples{}".format(sample_labels[i]),
                                 showlegend = False,
                                 name = sample_labels[i], 
                                 line = {"color":sample_color_dict[sample_labels[i]]}), 
                      row = 2, 
                      col = 1)

    ## If gene features found
    if shapes is not None:

        ## create gene track hover template
        hover_temp = []
        for row in exon_pos_df.itertuples():
             hover_temp.append( 
             ("<br><b>Gene<b> : {}<b>"
             "<br>Gene ID: {}"
             "<br>Transcript ID: {}"
             "<br>Biotype: {}"
             "<br>strand: {}"
             "<br>Exon Start: {}"
             "<br>Exon End: {}"
             "<br>Exon Number: {}").format(row.gene_symbol,
                                       row.gene_id,
                                       row.transcript_id,
                                       row.gene_biotype,
                                       row.strand,
                                       row.exon_start,
                                       row.exon_end,
                                       row.exon_number))
                                       
        ## Subplot of gene track
        fig.add_trace(go.Scatter(x = exon_pos_df.x_position, 
                                 y = exon_pos_df.y_position, 
                                 mode = "markers",
                                 hovertemplate = hover_temp,
                                 showlegend = False,
                                 marker = {"color":"white"},
                                 hoverlabel={"bgcolor":"blue"}),
                      row = 3, 
                      col = 1)

        ## Add gene feature shapes
        #fig["layout"].update(shapes = shapes)

        ## Add stranded text
        strand_text = []
        strand_indicator_interval = round((region_dict["end"] - region_dict["start"]) / 20)
        for i,row in enumerate(text_pos_df.itertuples()):
            if i % strand_indicator_interval == 0 :
                strand_text.append({"x":row.x_position, "y":row.y_position, "text":row.text, "font":{"size":12}, "showarrow":False, "xref":"x3","yref":"y3"})

        fig["layout"].update(annotations = strand_text)
        fig.update_layout(yaxis3 = {"tickmode":"array", "ticktext": y_tick_labels, "tickvals": y_tick_ref_vals})

    fig["layout"].update(shapes = final_shapes)
    fig.update_layout(title_text = "Per Base Coverage Track for region {}:{}-{}".format(region_dict["chrom"], region_dict["start"], region_dict["end"]))
    fig.update_layout(legend_title_text = "Samples")
    fig.update_xaxes(title_text = "{} Region Position".format(region_dict["chrom"]), tickformat = "d", row=max_rows, col=1)
    fig.update_yaxes(title_text= "by Region, Sample Relative Z-Score", row =2, col = 1 )
    fig.update_yaxes(title_text= "per Base Coverage",row =1, col = 1 )
    return(fig.to_html(full_html = False, include_plotlyjs = "cdn"), len(low_cov_shapes))


def intersect_gtf_pr(gtf_pr, 
                     region_dict):
    """
    intersect_gtf_pr
    ================
    Method to get a gene features for a specific region. 

    Parameters:
    ----------
    1) gtf_pr: (PyRanges Object) A gtf/gff read into a pyranges object
    2) region_dict: (dict) A dictionary representing the region to interesct the pyranges object with. Required keys incldue: "chrom", "start", and "end"

    Returns:
    +++++++
    1) (PyRanges Object) The pyranges gtf object that intersected the query region
    """
    
    region_gtf_pr = gtf_pr.intersect(pr.from_dict({"Chromosome":[region_dict["chrom"]],
                                                   "Start":[region_dict["start"]],
                                                   "End":[region_dict["end"]]}))

    return(region_gtf_pr)


def get_gtf_region_position_info(region_gtf_pr):
    """
    get_gtf_region_position_info
    ============================
    Method to convert a PyRanges GTF/GFF object into a dict of gene features

    Parameters:
    -----------
    1) region_gtf_pr: (PyRanges Object) A pyranges gtf object to extract gene features from

    Returns:
    +++++++
    1) gene_info_dict: (dict) A dictionary of gene features from the pyranges gtf object. The dict is formated as such:
                            {key = gene_id: 
                                {key = transcript_id:
                                    {key = chrom: value = Current Chrom},
                                    {key = start: value = feature start},
                                    {key = end: value = feature end},
                                    {key = gene_symbol: value = feature gene symbol},
                                    {key = biotype: value = feature gene biotype},
                                    {key = strand: value = feature strand},
                                    {key = exons:
                                        [{key = start: value = feature start, key = end: value = feature end, key = exon_number: value = feature exon number}]
                                    }
                                }
                            }

    NOTE: If no features are available for the current region, an empty dict will be returned
    """
    ## Check if region has no features
    if region_gtf_pr.empty:
        return(dict())
        

    gene_info_dict = dict()
    for name, group in region_gtf_pr.df.groupby("transcript_id"):
        for row in group.itertuples():

            ## Add Transcript Info
            if row.Feature == "transcript":

                if row.gene_id not in gene_info_dict:
                    gene_info_dict[row.gene_id] = dict()
                if row.transcript_id not in gene_info_dict[row.gene_id]:
                    gene_info_dict[row.gene_id][row.transcript_id] = {"exons":[]}

                gene_info_dict[row.gene_id][row.transcript_id]["chrom"] = row.Chromosome
                gene_info_dict[row.gene_id][row.transcript_id]["start"] = row.Start
                gene_info_dict[row.gene_id][row.transcript_id]["end"] = row.End
                gene_info_dict[row.gene_id][row.transcript_id]["gene_symbol"] = row.gene_name
                gene_info_dict[row.gene_id][row.transcript_id]["biotype"] = row.gene_type 
                gene_info_dict[row.gene_id][row.transcript_id]["strand"] = row.Strand

            ## Add exon feature info 
            elif row.Feature == "exon":

                if row.gene_id not in gene_info_dict:
                    gene_info_dict[row.gene_id] = dict()
                if row.transcript_id not in gene_info_dict[row.gene_id]:
                    gene_info_dict[row.gene_id][row.transcript_id] = {"exons":[]}

                gene_info_dict[row.gene_id][row.transcript_id]["exons"].append({"start":row.Start,"end":row.End,"exon_number":row.exon_number}) 

    return(gene_info_dict)


def extract_plotly_plot_from_html(html_string,
                                  div_id):
    """
    extract_plotly_plot_from_html
    =============================
    Method to extract the Plotly javascript string from an html string create from the python plotly .to_html() function 

    Only the Plotly javascript string will be extract and returend. The rest of the html string will be discarded. 

    The div id for the Plotly plot will be updated based on the "div_id" provided

    Parameters:
    -----------
    1) html_string: (str) An html string create for a plotly plot using the python plotly function '.to_html'
    2) div_id: (str) A div id to use for the extracted plot

    Returns:
    +++++++
    1) (str) A string that represents the plotly javascript script for a plot with an updated div id 
    """
    
    plot_string = "Plotly.newPlot("
    end_brace = ")"
    ## Extract plotly js 
    plotly_plot_string = html_string.strip().split(plot_string)[1].split(end_brace)[-2]
    ## Replace the div id
    plotly_plot_string = plotly_plot_string.replace(plotly_plot_string[0:plotly_plot_string.find(",")], "'"+div_id+"'")
        
    ## Return the plotly js string
    return(plot_string + plotly_plot_string + end_brace) 


def get_regions_from_genes(gene_list, 
                           gtf_pr):
    """
    get_regions_from_genes
    ======================
    Method to identify a region based on a gene symbol 

    Parameters:
    -----------
    1) gene_list: (list) A list of of gene symbols to use to extract region information for  
    2) gtf_pr: (PyRanges GTF Object) A gtf/gff file loaded into a pyragnes object

    Returns:
    ++++++++
    1) regions: (list) A list of a tabix styled regions based on the gene sybmol coordiantes in the gtf file
    """

    regions = []
    
    bad_genes = []
    for gene in gene_list:

        gene_pr = gtf_pr[gtf_pr.gene_name == gene]

        if gene_pr.empty:
            bad_genes.append(gene)
            continue

        chrom = gene_pr.df.Chromosome.to_list()[0]
        start = gene_pr.df.Start.min() - 100 
        end = gene_pr.df.End.max() + 100

        regions.append("{}:{}-{}".format(chrom,start,end))

    if bad_genes:
        print("\n!!ERROR!! At least one gene from the list was not found in the gtf file. Please make sure the gene symbol provided is correct and in the gtf file. If the symbol is a correct symbol, check for alternative gene symbols in the gtf file.")
        print("Bad Gene(s):\n\t- {}\n".format("\n\t- ".join(bad_genes))) 
        sys.exit(1)

    return(regions)
        

#---------------------------------------------------------------------------------------------------------------------------------
## Templates
#---------------------------------------------------------------------------------------------------------------------------------
    
## Template for regions added to an html file (This populations the region selector button)
available_regions_template = """
                      <a class="dropdown-item" href="#{region}">{region}</a>"""

## Template for each regions divs
div_template = """    <div class="container-fluid p-3">
        <div class="row pt-3">
            <div class="col-6">
                <h4 id="{reg}"> REGION: {reg}</h4>
                <div id="{reg}_1"></div>
            </div>
            <div class="col-sm-6">
                <div id="{reg}_2"></div>
            </div>
        </div>
        <div class="row pt-3">
            <div class="col-6">
                <h4>      </h4>
                <div id="{reg}_3"></div>
            </div>
            <div class="col-6">
                <div id="{reg}_4"></div>
            </div>
        </div>
        <div class="row pt-3">
            <div class="col-12">
                <h4>  </h4>
                <div class=big_div id="{reg}_5"></div>
            </div>
        </div>
        {toggle_div}
        <div class="row">
            <div class="col-lg-12">
                <h4>  </h4>
                <div class=big_div id="{reg}_6"></div>
            </div>
        </div>
        <div class="row h-10">
            <div class=row-boardered></div>
        </div>
    </div>
    <a href="#top" class="btn btn-outline-primary">Back to Top</a>
    <hr style="color:blue;background-color:black;border-width:2.5">

"""

## Template for low coverage toggle div
low_covearge_toggle_event_div = """

        <div class="row">
            <div class="col-10">
            </div>
            <div class="col-2">
                <h5><u> Low Coverage Bars </u></h5>
                <input id="{toggle_reg}_toggle-event" type="checkbox" checked data-toggle="toggle" data-on="Displaying" data-off="NOT Dispalying"  data-size="small" data-onstyle="danger", data-offstyle="info">
            </div>
        </div>
"""

## Template turning on and off low coverage bars
low_coverage_toggle_js_tempalte = """

$$("#${toggle_reg}_toggle-event").change(function() {

    toggled_on = $$(this).prop('checked')

    var div = document.getElementById("${reg}_6");

    var update = {}
    for (i = 0; i < ${shape_num}; i++){
        update[`shapes[$${i}].visible`] = toggled_on;
    }
    Plotly.relayout(div,update)
})
"""


#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def main():
    args = get_args()

    ## Check for region or gene argument
    if not args.region and not args.region_file and not args.gene and not args.gene_file:
        print("\n!!ERROR!! Missing region or gene parameter. Please add one of '--region', '--gene', '--region-file', '--gene-file'")
        sys.exit(1)

    ## Check for only one region/gene argument
    if sum(x is not None for x in [args.region, args.region_file, args.gene, args.gene_file]) > 1:
        print("\n!!ERROR!! To many region or gene parameters provided."
              "\n\t Arguments Provided: \n\t\t{}"
              "\n\tPlease choose only one of '--region', '--gene', '--region-file', '--gene-file'".format("\n\t\t".join([x for x in ["--region" if args.region else args.region, 
                                                                                                                                     "--region-file" if args.region_file else args.region_file, 
                                                                                                                                     "--gene" if args.gene else None, 
                                                                                                                                     "--gene-file" if args.gene_file else args.gene_file] 
                                                                                                                         if x is not None
                                                                                                                        ])
                                                                                                        )
             )  
        sys.exit(1)

    ## Set up html output
    if args.combine:
        print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print("**WARNING** The '--combine' argument has been set. All regions will be written to a single html file. This may make the html file to large to load")
        print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

        output_file_path = args.output if args.output.endswith(".html") else args.output + ".html" 

        print("\ncombined html file will be written to '{}'".format(output_file_path))

    else:
        output_dir = ""
        temp_dir = os.path.join(os.getcwd(),"html_output")
        dirname = os.path.dirname(args.output)
        if os.path.exists(dirname):
            if os.path.isdir(dirname):
                output_dir = dirname
            else:
                print("\nThe directory provided does not exists: '{}'".format(dirname))
                print("A temporary directory will be used: '{}'".format(temp_dir))
                output_dir = temp_dir
        else:
            print("\nThe directory provided does not exists: '{}'".format(dirname if dirname else "(No directory provided)"))
            print("A temporary directory will be used: '{}'".format(temp_dir))
            output_dir = temp_dir

            if not os.path.exists(output_dir):
                try:
                    os.mkdir(output_dir)
                except OSError as e:
                    print("\n!!ERROR!! Problem creating outputdir")
                    print(str(e))
                    sys.exit(1)

        print("\nOutput files will be written to '{}'".format(output_dir))
        basename = os.path.basename(args.output) if args.output.endswith(".html") else os.path.basename(args.output) + ".html"
        print("Each region html file will be written with the following convention: '<region>_{}'".format(basename))

    ## Load gtf/gff file 
    print("\nLoading the GTF/GFF file")
    gtf_pr = pr.read_gff(args.gtf)

    ## Gather regions to get coverate for 
    gene_list = []
    print("\nGathering Regions")
    if args.region:
        regions = args.region
        gene_list = [""] * len(regions)
    elif args.gene:
        regions = get_regions_from_genes(args.gene, gtf_pr)
        gene_list = args.gene
    elif args.region_file:
        try:
            with open(args.region_file) as rf:
                regions = [x.strip() for x in rf]
                gene_list = [""] * len(regions)
        except IOError as e:
            print("\n!!ERROR!! Unable to open the region file '{}'. Please make sure it exists".format(args.region_file))
    elif args.gene_file:
        try:
            with open(args.gene_file) as rf:
                genes = [x.strip() for x in rf]
                gene_list = genes
            regions = get_regions_from_genes(genes, gtf_pr)
        except IOError as e:
            print("\n!!ERROR!! Unable to open the region file '{}'. Please make sure it exists".format(args.region_file))

    ## Get samples 
    samples = [x.strip().split(".")[0] for x in args.input]
    print("\nProcessing region coverage for samples:\n\t- {}\n".format("\n\t- ".join(samples)))
    
    ## Get colors
    random_colors = True
    sample_colors = dict()
    if args.sample_colors:
        if len(args.sample_colors) != len(samples):
            print("\n**WARNING** The number of colors provided does not match the number of samples. Random colors will be chosen")
            print("\t- Number of Samples: {}".format(len(samples)))
            print("\t- Number of Colors:  {}".format(len(args.sample_colors)))

        else:
            sample_colors = dict(zip(samples, args.sample_colors))
            random_colors = False
        
    if random_colors:
        ## Get a random color for each sample 
        all_colors = [x for k,x in pltc.cnames.items()]
        for s in samples:
            sample_colors[s] = all_colors.pop(random.randint(0,len(all_colors) -1))


    print("\nAdding low coverage annotation based on coverage <= {}".format( args.lch_cutoff))

    ## Temnplate strings
    available_regions = ""
    region_divs = ""
    plotly_plots = ""

    ## Iterate over each region
    for i,region in enumerate(regions):

        region_dict = get_dict_from_region_string(region)
        print("\nProcessing Region: {} {}:{}-{}".format(gene_list[i], region_dict["chrom"], region_dict["start"], region_dict["end"]))

        ## Intersect gtf with region
        region_gtf_pr = intersect_gtf_pr(gtf_pr, region_dict)

        ## get a dict of region position info from gtf file
        gtf_region_dict = get_gtf_region_position_info(region_gtf_pr)

        ## Get per base coverage per sample
        by_sample_coverage = []
        sample_list = []
        for sample_file in args.input:

            sample_object = pysam.TabixFile(sample_file)
            
            sample_list.append(sample_file.strip().split(".")[0])
            by_sample_coverage.append(get_sample_region_coverage(sample_object, region_dict))
            
            get_region_coverage_proportion(get_sample_region_coverage(sample_object, region_dict))

            sample_object.close()

        base_div_id = "{}{}:{}-{}".format(gene_list[i] + ": " if gene_list[i] != "" else gene_list[i], region_dict["chrom"], region_dict["start"], region_dict["end"])
        ## Plot region coverage info per sample 
        print("\t- Plotting KDE coverage distribution")
        kde_html_string = plot_kde_per_sample(by_sample_coverage, 
                                              sample_list, 
                                              sample_colors, 
                                              region_dict)
        kde_plotly_js = extract_plotly_plot_from_html(kde_html_string, base_div_id + "_1")
        plotly_plots += kde_plotly_js + "\n\n" 

        print("\t- Plotting KDE z-score coverage distribution")
        kde_z_score_html_string, z_scores_per_base_per_sample = plot_z_score_distribution(by_sample_coverage, 
                                                                                          sample_list, 
                                                                                          sample_colors, 
                                                                                          region_dict)
        kde_z_score_plotly_js = extract_plotly_plot_from_html(kde_z_score_html_string, base_div_id + "_3")
        plotly_plots += kde_z_score_plotly_js + "\n\n" 

        print("\t- Generating General Stats Table")
        table_html_string, max_coverage_value, has_low_coverage = get_plotly_table(by_sample_coverage, 
                                                                                   sample_list, 
                                                                                   args.lch_cutoff)
        table_js = extract_plotly_plot_from_html(table_html_string, base_div_id + "_2")
        plotly_plots += table_js + "\n\n" 


        print("\t- Plotting proportion of bases covered")
        proportion_html_string = plot_porportion_coverage(by_sample_coverage, 
                                                          sample_list, 
                                                          sample_colors, 
                                                          region_dict)
        proportion_plotly_js = extract_plotly_plot_from_html(proportion_html_string, base_div_id + "_4")
        plotly_plots += proportion_plotly_js + "\n\n" 

        print("\t- Plotting sample vs sample per base coverage")
        sample_vs_sample_html_string = plot_sample_vs_sample_per_base_coverage(by_sample_coverage, 
                                                                               sample_list, 
                                                                               sample_colors, 
                                                                               region_dict)
        sample_vs_sample_plotly_js = extract_plotly_plot_from_html(sample_vs_sample_html_string, base_div_id + "_5")
        plotly_plots += sample_vs_sample_plotly_js + "\n\n" 

        print("\t- Plotting regional coverage by base per sample")
        region_html_string, count_low_cov_shapes  = plot_per_base_coverage(by_sample_coverage, 
                                                                           sample_list, 
                                                                           sample_colors, 
                                                                           region_dict, 
                                                                           gtf_region_dict, 
                                                                           args.lch_cutoff, 
                                                                           max_coverage_value,
                                                                           has_low_coverage)
        region_plotly_js = extract_plotly_plot_from_html(region_html_string, base_div_id + "_6")
        plotly_plots += region_plotly_js + "\n\n" 


        ## Add template info
        print("\t- Updating html template")
        available_regions += available_regions_template.format(region = base_div_id)

        ## Add low coverage toggle
        if has_low_coverage:
            ## Add toggle 
            plotly_plots += string.Template(low_coverage_toggle_js_tempalte).substitute(toggle_reg = base_div_id.strip().replace(" ","").replace(":","-"), 
                                                                                        reg = base_div_id, 
                                                                                        shape_num = count_low_cov_shapes) + "\n\n"

            region_divs += div_template.format(reg = base_div_id, toggle_div = low_covearge_toggle_event_div.format(toggle_reg = base_div_id.strip().replace(" ","").replace(":","-"), reg = base_div_id))
        else:
            region_divs += div_template.format(reg = base_div_id, toggle_div = "") 

        ## If not combining html files
        if not args.combine:

            print("\n\tCreating html file")
            tmpl_dict = dict()

            tmpl_dict["region_references"] = available_regions
            tmpl_dict["region_div"] = region_divs
            tmpl_dict["plotly_plots"] = plotly_plots

            try:
                tmpl = string.Template(open(args.tmpl).read())
            except OSError as e:
                print("!!ERROR!! Problem loading the template html file. Please make sure it exists. Template file = '{}'".format(args.tmpl))
                print(str(e))
                sys.exit(1)

            html_output = args.output if args.output.endswith(".html") else args.output + ".html"

            prefix = gene_list[i] if gene_list[i] != "" else "{}-{}-{}".format(region_dict["chrom"], region_dict["start"], region_dict["end"])
            output_file_path = os.path.join(output_dir,"{}_{}".format(prefix,os.path.basename(args.output)))
            print("\n\tWriting output to '{}'".format(output_file_path))
            with open(output_file_path, "w") as html_out:
                html_out.write(tmpl.safe_substitute(**tmpl_dict))
                    
            ## Reset  temnplate strings
            available_regions = ""
            region_divs = ""
            plotly_plots = ""
            

    ## If combining html plots into a single html file
    if args.combine:
        print("\nCreating combined html file")
        tmpl_dict = dict()

        tmpl_dict["region_references"] = available_regions
        tmpl_dict["region_div"] = region_divs
        tmpl_dict["plotly_plots"] = plotly_plots

        try:
            tmpl = string.Template(open(args.tmpl).read())
        except OSError as e:
            print("!!ERROR!! Problem loading the template html file. Please make sure it exists. Template file = '{}'".format(args.tmpl))
            print(str(e))
            sys.exit(1)

        print("\nWriting output to '{}'".format(output_file_path))
        with open(output_file_path, "w") as html_out:
            html_out.write(tmpl.safe_substitute(**tmpl_dict))
        
    print("\nDONE")
            

if __name__ == '__main__':
    main()

