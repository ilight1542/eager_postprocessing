import argparse
import csv
import pysam
import matplotlib.pyplot as plt

## eventually take in a list of bams (single file with pointers?)

parser = argparse.ArgumentParser(description='create edit distance plot from bam imputs')
parser.add_argument('-bams', '--bams', metavar=".bam file(s)", required=True, nargs='+', help="bam inputs")
parser.add_argument('--relative', metavar="reference --> lengths", nargs='?', help="compute edit distance relative to genome length, reference --> lengths, if by_chrom just input T")
parser.add_argument('--output', action='store_true', help="outputs plot as png file")
parser.add_argument('--file_name', nargs=1, help="file name output if printing")
parser.add_argument('--name_loc', metavar="pos", nargs='?', default=4, help="position of sample name, (num of / from end)")
parser.add_argument('--by_chrom', action='store_true',help="parse edit distances by chromosome")
parser.add_argument('--single_reference', metavar='name', nargs=1, help="name of reference all bams were mapped to")
args=parser.parse_args()

def get_edit_dist(bam, by_chrom=False):
    """parses bam files to get edit distances for reads mapped to reference"""
    edit_dists=[]
    list_max=0
    if by_chrom:
        references=bam.references
        by_chrom_edit={ref:[] for ref in references}
        for ref in references:
            for i in bam.fetch(ref):
                current_edit_dist=i.get_tags()[6][1]
                if current_edit_dist > list_max:
                    list_max=current_edit_dist
                by_chrom_edit[ref].append(current_edit_dist)
        return (by_chrom_edit, list_max)
    for i in bam.fetch():
        current_edit_dist=i.get_tags()[6][1]
        if current_edit_dist > list_max:
            list_max=current_edit_dist
        edit_dists.append(current_edit_dist)
    return (edit_dists, list_max)

def bar_plot(ax, data, colors=None, total_width=0.8, single_width=1, legend=True):
    """
    adapted from  https://stackoverflow.com/questions/14270391/python-matplotlib-multiple-bars user pascscha
    Draws a bar plot with multiple bars per data point.

    Parameters
    ----------
    ax : matplotlib.pyplot.axis
        The axis we want to draw our plot on.

    data: dictionary
        A dictionary containing the data we want to plot. Keys are the names of the
        data, the items is a list of the values.

        Example:
        data = {
            "x":[1,2,3],
            "y":[1,2,3],
            "z":[1,2,3],
        }

    colors : array-like, optional
        A list of colors which are used for the bars. If None, the colors
        will be the standard matplotlib color cyle. (default: None)

    total_width : float, optional, default: 0.8
        The width of a bar group. 0.8 means that 80% of the x-axis is covered
        by bars and 20% will be spaces between the bars.

    single_width: float, optional, default: 1
        The relative width of a single bar within a group. 1 means the bars
        will touch eachother within a group, values less than 1 will make
        these bars thinner.

    legend: bool, optional, default: True
        If this is set to true, a legend will be added to the axis.
    """

    # Check if colors where provided, otherwhise use the default color cycle
    if colors is None:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Number of bars per group
    n_bars = len(data)

    # The width of a single bar
    bar_width = total_width / n_bars

    # List containing handles for the drawn bars, used for the legend
    bars = []

    # Iterate over all data
    for i, (name, values) in enumerate(data.items()):
        # The offset in x direction of that bar
        x_offset = (i - n_bars / 2) * bar_width + bar_width / 2

        # Draw a bar for every value of that type
        for x, y in enumerate(values):
            bar = ax.bar(x + x_offset, y, width=bar_width * single_width, color=colors[i % len(colors)])

        # Add a handle to the last drawn bar, which we'll need for the legend
        bars.append(bar[0])

    # Draw legend if we need
    if legend:
        ax.legend(bars, data.keys())

def reference_sizes(ref_size_file):
    with open(ref_size_file, newline="\n") as output:
        csv_list=list(csv.reader(output, delimiter=" "))
    name_to_len_dict={}
    for line in csv_list:
        name_to_len_dict[line[0]] = int(line[1])
    return name_to_len_dict

def output_plot(input, loc, single_reference=None, reference_size_file=None, relative=False, output=False, filename=None, by_chrom=False):
    edit_dist_list={}
    max=0
    edit_dist_summary={}
    for bams in input:
        alignment_file=pysam.AlignmentFile(str(bams),"rb")
        if by_chrom:
            edit_dist_output=get_edit_dist(alignment_file, True)
            edit_dist_list=edit_dist_output[0]
            name=str(bams.split("/")[-loc])
            name+="_by_chromosome"
        else:
            name=str(bams.split("/")[-loc]) ## this requires the naming scheme to be the same across inputs
            edit_dist_output=get_edit_dist(alignment_file)
            edit_dist_list[name]=edit_dist_output[0]
        if edit_dist_output[1] > max:
            max = edit_dist_output[1]
    for bams in edit_dist_list.keys():
        list_of_edit_dist=[]
        for i in range(max):
            if relative:
                if by_chrom:
                    reference_sizes_list=alignment_file.lengths
                    reference_names=alignment_file.references
                    reference_sizes_dict={}
                    for refs in range(len(reference_names)):
                        reference_sizes_dict[reference_names[refs]]=reference_sizes_list[refs]
                    list_of_edit_dist.append(edit_dist_list[bams].count(i)/reference_sizes_dict[bams])
                else:
                    reference_sizes_dict=reference_sizes(reference_size_file)
                    list_of_edit_dist.append(edit_dist_list[bams].count(i)/reference_sizes_dict[bams])
            else:
                list_of_edit_dist.append(edit_dist_list[bams].count(i))
        edit_dist_summary[bams]=list_of_edit_dist
    fig, ax = plt.subplots()
    bar_plot(ax, edit_dist_summary)
    # Check if single reference used
    if single_reference is not None:
        header=str(single_reference)
        fig.suptitle("Edit Distances for "+header)
    else:
        fig.suptitle("Edit Distances per Reference")
    plt.xticks(range(max))
    plt.xlabel("Edit Distance")
    if relative:
        plt.ylabel("Reads relative to Reference size")
    else:
        plt.ylabel("Reads")
    if not output:
        plt.show()
    else:
        if filename != None:
            plt.savefig(fname=str(filename[0]))
        else:
            plt.savefig(fname=(name+".jpg"))

if __name__ == '__main__':
    if args.single_reference != None:
        single_ref=args.single_reference[0]
    else: 
        single_ref=None
    if args.relative != None:
        output_plot(args.bams, loc=int(args.name_loc),single_reference=single_ref, reference_size_file=args.relative, relative=True, output=args.output, filename=args.file_name, by_chrom=args.by_chrom)
    else:
        output_plot(args.bams,loc=int(args.name_loc),single_reference=single_ref, output=args.output, filename=args.file_name, by_chrom=args.by_chrom)