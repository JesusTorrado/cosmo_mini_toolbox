# Common imports
import sys
import argparse
import json
import matplotlib.pyplot as plt

# Local imports
sys.path.append("../src")
from Chain import Chain
from plot_lik import plot_lik_2D

# Parsing input
parser = argparse.ArgumentParser(
    description="Plot the likelihood sampling of a chain in a 2D projection.")
parser.add_argument("-k", "--keywords", type=str, dest="kwargs",
    metavar=("Optional list of keyword arguments of the 'plot_lik_2D' function, "+
             """as a dictionary; e.g. '{"format":"-loglik"}'. """ +
             "The order of the quotes is important!"))
parser.add_argument("-p", "--params", type=str, dest="params", nargs=2,
    required=True,
    metavar="List of 2 parameters to plot")
parser.add_argument("-o", "--output", type=str, dest="output", nargs=1,
    default=None,
    metavar="Optional file in which to store the plot.")
parser.add_argument("folders", type=str, nargs="+",
    metavar="Folders of the chains to be plotted.")
args = parser.parse_args()

# Loading the chains
chains = [Chain(folder) for folder in args.folders]

# Plotting
try:
    kwargs = json.loads(args.kwargs)
except ValueError:
    raise ValueError("The dictionary after '-k'/'--keywords' could not be parsed!" +
                     " Check the input.")
fig, axarr = plt.subplots(2,2)
axes_locations = {                    "profile":  axarr[0,1],
                  "mean": axarr[1,0], "marginal": axarr[1,1]}
for mode, axes in axes_locations.items():
    ax, options = plot_lik_2D(mode, chains, params=args.params, #format = "-loglik", dpi=200, fontsize_labels=14, fontsize_ticks=8,
                save=0, axes=axes, **kwargs
                )
    axes.set_title(mode.title(), fontdict={'fontsize':16})

# Text:
text = ("Chain:\n %s\n\n"%[c.name() for c in chains] +
        "Parameters:\n%s"%(args.params))
axarr[0,0].set_axis_off()
axarr[0,0].text(0, 1, text, weight="bold", verticalalignment="top")

# Plot
plt.tight_layout()
if args.output:
    plt.savefig(args.output[0], **options)
else:
    plt.show()
plt.close()
