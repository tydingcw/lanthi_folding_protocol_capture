#!/usr/bin/env python
import os
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from math import exp
from argparse import ArgumentParser

def main():

    # Parse command line
    args = cmdlineparse()

    # Load terms to be plotted against one another
    file_headers = np.genfromtxt(args.input, skip_header=0, dtype='string', max_rows=1, delimiter=args.delimiter)
    index = int(0)

    min_disc_bins = int(args.funnel_disc_min_bin_size)

    # First we need to get all of the receptor-ligand gas-phase interaction terms
    for col in file_headers:
        if (col == str(args.x_str)):
            x_col_num = index
        elif(col == str(args.y_str)):
            y_col_num = index
        index = index + 1

    index = int(0)
    if args.z_str is not None:
        for col in file_headers:
            if(col == str(args.z_str)):
                z_col_num = index
            index = index + 1
 
   # Collect the specified columns' data
    col_nums = []
    col_nums.append(x_col_num)
    col_nums.append(y_col_num)
    if args.z_str is not None:
        col_nums.append(z_col_num)
    data = np.genfromtxt(args.input, skip_header=1, dtype='float', usecols=col_nums, delimiter=args.delimiter)
    data = np.transpose(data)
    x = data[0]
    y = data[1]
    if args.z_str is not None:
        z = data[2]

    ofile = open(args.output, "w+")

    # Specify axis cutoffs and identify data that extend beyond each axis
    if args.x_range:
        xmin = int(args.x_range[0])
        xmax = int(args.x_range[1])
        x_bad_index = np.array([int(row) for row in range(0,len(x)) if x[row] < xmin or x[row] > xmax])
    else:
        x_bad_index = np.array([])
    if args.y_range:
        ymin = int(args.y_range[0])
        ymax = int(args.y_range[1])
        y_bad_index = np.array([int(row) for row in range(0,len(y)) if y[row] < ymin or y[row] > ymax])
    elif args.y_max:
        ymin = float(np.amin(y))
        ymax = int(args.y_max)
        y_bad_index = np.array([int(row) for row in range(0,len(y)) if y[row] < ymin or y[row] > ymax])
    elif args.y_standard:
        ymin = float(np.amin(y))
        ymax = int(0)
        y_bad_index = np.array([int(row) for row in range(0,len(y)) if y[row] < ymin or y[row] > ymax])
    else:
        y_bad_index = np.array([])

    # Exclude data that extends beyond axes
    if (args.x_range or args.y_range or args.y_max or args.y_standard) and args.exclude_axis_limits:
        bad_indices = np.unique(np.sort(np.concatenate((x_bad_index, y_bad_index))))
        x = np.delete(x, bad_indices)
        y = np.delete(y, bad_indices)
        z = np.delete(z, bad_indices)

    # Set global plot options
    font = {'family' : 'sans-serif',
    'weight' : 'normal',
    'size'   : args.font_size}
    plt.rc('font', **font)

    # Plot data
    if args.z_str is not None:
        sc = plt.scatter(x, y, marker='.', c=z, edgecolors='black', linewidths=0.10, cmap=str(args.cmap), alpha=1.0)
        cb = plt.colorbar(sc)
        cb.locator = MaxNLocator(integer=True)
        cb.update_ticks()
        cb.set_label('Energy (REU)')
    else:
        plt.scatter(x, y, marker='o', color=str(args.color), edgecolors='black', alpha=0.5)


    # Set axes options
    plt.xlabel(args.x_label, fontsize=args.font_size)
    plt.ylabel(args.y_label, fontsize=args.font_size)
    if args.x_range:
        plt.xlim([xmin - 1, xmax])
    if args.y_range:
        plt.ylim([ymin - 1, ymax])
    elif args.y_max:
        plt.ylim([ymin - 1, ymax])
    elif args.y_standard:
        plt.ylim([ymin - 1, ymax])

    # Funnel discrimination metric
    if args.funnel_score:
        funnel_discrimination = ComputeFunnelScore(x, y, min_disc_bins)
        odatafile = open(args.output+".dat", "w")
        odatafile.write("Funnel discrimination score: \n")
        odatafile.write('%0.2f' % np.mean(funnel_discrimination) + '\n')
        odatafile.close()
    if args.pnear:
        pnear = CalculatePNear(y, x, lambda_val=float(args.pnear_lambda), kbt=float(args.pnear_kbt))
        odatafile = open(args.output+".dat", "a")
        odatafile.write("Pnear value: \n")
        odatafile.write('%0.2f' % pnear + '\n')
        odatafile.close()

    # Setup plot
    if args.funnel_score and args.pnear:
        plt.legend(['Funnel score: %0.2f' % np.mean(funnel_discrimination) + "\n" 'Pnear value: %0.2f' % pnear],
                   loc='upper left')
    elif args.funnel_score:
        plt.legend(['Funnel score: %0.2f' % np.mean(funnel_discrimination)], loc='upper left')
    elif args.pnear:
        plt.legend(['Pnear value: %0.2f' % pnear], loc='upper left')

    # Save the figure to output and be done
    plt.savefig(ofile, bbox_inches='tight', dpi=600)
    ofile.close()
    return 0

# Funnel score
def ComputeFunnelScore(RMS, SCORE, MIN_DISC_BIN_SIZE):
    # Adapted from scripts by Sevy et al. 2020, DOI: 10.1016/j.str.2020.04.005
    # Sevy et al. funnel score is modified from Conway et al. 2013, DOI: 10.1002/pro.2389
    rms, score = NormalizeXY(RMS, SCORE)
#     rms_bins = [0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0,5.0,6.0]
    rms_bins = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0]

    bin_scores = []
    for bin in rms_bins:

        below_bin = [sc for rm, sc in zip(rms, score) if rm < bin]
        above_bin = [sc for rm, sc in zip(rms, score) if rm > bin]

        if len(below_bin) < MIN_DISC_BIN_SIZE or \
                len(above_bin) < MIN_DISC_BIN_SIZE:
            bin_scores.append(0.0)
            continue

        min_below_bin = min(below_bin)
        min_above_bin = min(above_bin)

        gap = min_below_bin - min_above_bin
        bin_scores.append(gap)

    return bin_scores


def NormalizeXY(X, Y):
    x, y = zip(*sorted(zip(X, Y), key=lambda a: a[0]))
    p_95 = np.percentile(y, 95)
    p_5 = np.percentile(y, 5)
    min_value = min(y)
    y = [float(z - min_value) / (p_95 - p_5) for z in y]
    return x, y

# NOTE
# The following function for calculation of PNear is adapted from code written by Vikram Mulligan:
#
# AUTHOR
# Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org).
#
# DATE
# Created 10 October 2019.
#
# DESCRIPTION (paraphrased)
# PNear ranges from 0 to 1, with 0 representing no propensity to favour the desired conformation and 1 representing
# 100% propensity to favour that conformation.
# The PNear computation takes three inputs: data (scores and RMSDs), a value "lambda", and a value "kbt".
# The lambda value,  in Angstroms, indicates the extent to which
# a structure may deviate from the native state and still be considered "native-like" (e.g. 1.5 to 2.0 A for peptides,
# and 2.0 to 4.0 A for proteins). The kbt value, in kcal/mol, determines how large the energy
# gap between native and nonnative must be in order for a design to be considered a good folder.
# Typically, 0.5961 (300K) to 0.62 (310K) are used.

# PNear value
def CalculatePNear(scores, rmsds, lambda_val=1.5, kbt=0.62):
    nscores = len(scores)
    assert nscores == len(rmsds), "Error in calculate_pnear(): The scores and rmsds lists must be of the same length."
    assert nscores > 0, "Error in calculate_pnear(): At least one score/rmsd pair must be provided."
    assert kbt > 1e-15, "Error in calculate_pnear(): kbt must be greater than zero!"
    assert lambda_val > 1e-15, "Error in calculate_pnear(): lambda must be greater than zero!"
    minscore = min(scores)
    weighted_sum = 0.0
    Z = 0.0
    lambdasq = lambda_val * lambda_val
    for i in range(nscores):
        val1 = exp(-(rmsds[i] * rmsds[i]) / lambdasq)
        val2 = exp(-(scores[i] - minscore) / kbt)
        weighted_sum += val1 * val2
        Z += val2
    assert Z > 1e-15, "Math error in calculate_pnear()!  This shouldn't happen."
    return weighted_sum / Z


# Argument parser
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments")
    parser.add_argument("-input", dest="input", required=True, help="Input file; 2 column file where the first column is RMSD and the second column is score", metavar="<input file>")
    parser.add_argument("-output", dest="output", required=True, help="Name of output file", metavar="<output file>")
    parser.add_argument("-delimiter", dest="delimiter", required=False,default=" ", help="Delimiter separating columns in input file", metavar="<delimiter>")
    parser.add_argument("-font_size", dest="font_size", required=False,default=14, help="Font size", metavar="<font size>")
    parser.add_argument("-x_str", dest="x_str", required=False,help="Specify x score term from input file", default="ligand_rms_no_super_X", metavar="<x_str>")
    parser.add_argument("-y_str", dest="y_str", required=False,help="Specify y score term from input file", default="interface_delta_X", metavar="<y_str>")
    parser.add_argument("-z_str", dest="z_str", required=False,help="Specify values used to color plot points", default="", metavar="<z_str>")
    parser.add_argument("-x_range", dest="x_range", required=False,help="Set x-axis range", nargs="+", metavar="<x range>")
    parser.add_argument("-y_range", dest="y_range", required=False,help="Set y-axis range", nargs="+", metavar="<y range>")
    parser.add_argument("-y_standard", dest="y_standard", default=False, required=False, help="Set y-axis range to Min=min(y),Max=0", action="store_true")
    parser.add_argument("-y_max", dest="y_max", default=False, required=False, help="Set y-axis range to Min=min(y),Max", metavar="<y_max>")
    parser.add_argument("-x_label", dest="x_label", required=False,default="RMSD ($\AA$)", help="Set x-axis label", metavar="<x label>")
    parser.add_argument("-y_label", dest="y_label", required=False,default="Interaction Energy (REU)", help="Set x-axis label", metavar="<y label>")
    parser.add_argument("-color", dest="color", required=False,default="lime", help="Facecolor of markers in scatter plot", metavar="<color>")
    parser.add_argument("-cmap", dest="cmap", required=False,default="ocean_r", help="Facecolor of markers in scatter plot if z-value specified", metavar="<cmap>")
    parser.add_argument("-exclude_axis_limits", dest="exclude_axis_limits", required=False, default=False, action="store_true", help="Do not plot points extending beyond axes")
    parser.add_argument("-funnel_score", dest="funnel_score", required=False, default=False, action="store_true",
                        help="Compute funnel score for each funnel")
    parser.add_argument("-funnel_disc_min_bin_size", dest="funnel_disc_min_bin_size", required=False, default=1,
                        help="Minimum number of counts per bin to have score below 0.0")
    parser.add_argument("-pnear", dest="pnear", required=False, default=False, action="store_true",
                        help="Compute pnear for each funnel")
    parser.add_argument("-pnear_lambda", dest="pnear_lambda", required=False, default=1.5,
                        help="indicates the extent to which a structure may deviate (in Angstroms) from the native state "
                             "and still be considered native-like (e.g. 1.5 to 2.0 A for peptides and 2.0 to 4.0 A for proteins)."
                        )
    parser.add_argument("-pnear_kbt", dest="pnear_kbt", required=False, default=0.62,
                        help="The kbt value, in kcal/mol, determines how large the energy gap between native and "
                             "nonnative must be in order for a design to be considered a good folder. "
                             "Typically, 0.5961 to 0.62 are used for room and physiological temperature, "
                             "respectively."
                        )
    args=parser.parse_args()
    return args

# Run main() from command line
if __name__ == '__main__':
    main()
