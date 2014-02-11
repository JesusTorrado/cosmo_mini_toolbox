###########################################################
# Tools for analyzing, plotting and doing stuff to chains #
###########################################################

import os
import numpy as np
import matplotlib.pyplot as plt
from math import floor

class Chain():
    """
    Class to load the parameters of a chain in memory (incl. the chain files,
    but not their content.

    """
    def __init__(self, folder=None, prefix=None, code="MontePython"):
        codes = ["montepython", "cosmomc"]
        assert code.lower() in codes, "Code not known. Known codes are "+str(codes)
        self._code = code.lower()
        # MontePython case
        if self._code == "montepython":
            self._load_params_montepython(folder)
            self._chains = [os.path.join(folder, a) for a in os.listdir(folder)
                            if a[-4:]==".txt"] 
        # CosmoMC case
        if self._code == "cosmomc":
            self._load_params_cosmomc(folder)
            self._chains = [os.path.join(folder, a) for a in os.listdir(folder)
                            if a[-4:]==".txt"] 
        # Points
        individual_chains = []
        for chain in self._chains:
            individual_chains.append(np.loadtxt(chain))
        self._points = np.concatenate(individual_chains)
        # Scaling -- MontePython
        for i, param in enumerate(self.varying_params() + self.derived_params()):
            self._points[:, i+2] *= float(self._params["parameter"][param][4])
        # Finding best fit(s) -- faster if done noe (only once!)
        maxloglik = self.get_min("mloglik")
        bf = [i for i, mloglik in enumerate(self.points("mloglik"))
              if mloglik == maxloglik]
        self._bf_points = [self.points()[i] for i in bf]

    # Load parameters
    def _load_params_montepython(self, chain_folder) :
        assert os.path.exists(chain_folder) ,\
            "The chain folder does not exist: "+chain_folder
        self._name = chain_folder.strip("/").split("/")[-1]
        logparam = open(os.path.join(chain_folder, "log.param"))
        cosmo_arguments = {}
        parameters      = {}
        path            = {}
        self._sorted_varying_params = [] # list of the (sorted) varying parameters
        self._sorted_derived_params = [] # list of the (sorted) derived parameters
        for line in logparam:
            if not line.strip():
                continue
            if line.strip()[0] == "#":
                continue
            if "data.cosmo_arguments" in line:
                if ".update" in line:
                    continue
                left, right = line.split("=")
                key = left[1+left.find("["):left.find("]")].strip("'").strip('"')
                cosmo_arguments[key] = right.split(";")[0].strip().strip("'").strip('"')
            if "data.parameters" in line:
                left, right = line.split("=")
                key = left[1+left.find("["):left.find("]")].strip("'").strip('"')
                parameters[key] = [a.strip() for a in
                    right.split(";")[0].strip().lstrip("[").rstrip("]").split(",")]
                param_type = parameters[key][5].strip('"').strip("'")
                if param_type in ["cosmo", "nuisance"]:
                    self._sorted_varying_params.append(key)
                if param_type == "derived":
                    self._sorted_derived_params.append(key)

            if "data.path" in line:
                left, right = line.split("=")
                if not "data.path" in left:
                    continue
                key = left[1+left.find("["):left.find("]")].strip("'").strip('"')
                path[key] = right.split(";")[0].strip()
        logparam.close()
        params = {"cosmo_argument": cosmo_arguments,
                  "parameter":      parameters,
                  "path":           path}
        self._params = params
    def _load_params_cosmomc(self, chain_folder) :
        assert os.path.exists(chain_folder) ,\
            "The chain folder does not exist: "+chain_folder
        self._name = chain_folder.strip("/").split("/")[-1]
        logparam = open(os.path.join(chain_folder, chain_folder.split("/")[-1]+".inputparams"))
        cosmo_arguments = {}
        parameters      = {}
        for line in logparam :
            if line[0] == "#" :
                continue
            # Varying parameters
            if "param[" in line and "=" in line :
                left, right = line.split("=")
                paramname = left.split("[")[-1].strip()[:-1]
                parameters[paramname] = [a.strip() for a in right.strip().split()]
            # Cosmo code arguments (i.e. fixed)
            if not("param[") in line and "=" in line :
                left, right = [a.strip() for a in line.split("=")]
                paramname = left
                cosmo_arguments[paramname] = right
        logparam.close()
        params = {"cosmo_arguments": cosmo_arguments,
                  "parameters":      parameters}
        self._params = params

    # Get chain data in a code independent way
    def name(self):
        return self._name
    def varying_params(self):
        return self._sorted_varying_params
    def derived_params(self):
        return self._sorted_derived_params
    def chain_files(self):
        return self._chains
    def index_of_param(self, param):
        """
        Returns the index of the given parameter in a chain point.
        """
        if param == "#":
            return 0
        elif param == "mloglik":
            return 1
        elif param in self.varying_params():
            return 2 + self.varying_params().index(param)
        elif param in self.derived_params():
            return (2 + len(self.varying_params()) +
                    self.derived_params().index(param))
        else:
            raise ValueError("Unrecognized parameter: '"+str(param)+"'.")
    def points(self, param=None):
        """
        Possibilities:
        
        * param == None (or not defined) : all chain points, as rows
        * param == "#" : number of steps / (non-normalizaed) prob. of the sample
        * param == "mloglik" : minus log-likelihood
        * param == <param_name> : chain points for said parameter

        """
        if not param:
            return self._points
        else:
            return self._points[:, self.index_of_param(param)]
    def get_min(self, param):
        """
        Gets the minimum value of the given parameter that the chain has reached.
        """
        return self.points(param=param).min()
    def get_max(self, param):
        """
        Gets the maximum value of the given parameter that the chain has reached.
        """
        return self.points(param=param).max()
    def get_limits(self, param):
        """
        Gets the limits *imposed* on the search for the given parameter
        (different from 'get_min', 'get_max').

        Output: [lower, upper], with 'None' where no limit was imposed.
        """
        if self._code == "montepython":
            if param in self.varying_params():
                limits = [None, None]
                for i in [1, 2]:
                   lim = self._params["parameter"][param][i].strip()
                   scale = float(self._params["parameter"][param][4])
                   limits[i-1] = float(lim)*scale if lim != '-1' else None
                return limits
            elif param in self.derived_params():
                raise ValueError("The parameter '%s' has no limits: "%param +
                                 "it is a derived parameter.")
            else:
                raise ValueError("The parameter '%s' is not recognised."%param)
        elif self._code == "cosmomc":
            raise NotImplementedError()
    def best_fit(self, param=None, more_than_one=False):
        """
        Returns the best fit point(s) of the chain, as a list of at least one point.
        (in case there is only more than one best fit point).

        A parameter can be specified to get only its best fit value.
        """
        if not param:
            return_value = self._bf_points
        else:
            return_value = [a[self.index_of_param(param)] for a in self._bf_points]
        if not(more_than_one) and len(return_value)>1:
            print ("WARNING: Only one best fit required, but more than one present"+
                   " with the same value of the likelihood.\n"+
                   "Use the keyword 'more_than_one=True' to get them all.")
            return return_value[0]
        elif not(more_than_one) and len(return_value)==1:
            return return_value[0]
        elif more_than_one:
            return return_value


### Plot of 2D profile likelihoods
def plot_2Dprofile(chains, params,
                   # Main parameters
                   labels=None, format="-loglik", central_mloglik=None,
                   limits=None, n_grid=100, aspect=1,
                   color_map="jet_r", black_and_white=False,
                   cb_orientation="vertical",
                   bf_show=1, regions_show=True, save_file=None,
                   # Fine tuning
                   fontsize_labels=18, fontsize_ticks=12,
                   cb_ticks_formatter=None, cb_shrink=float(1),
                   padding = 0.02, dpi=150, transparent=True, not_yet=False,
                   bf_alpha=1, bf_radius=1, bf_thickness=1,
                   bf_color_in="white", bf_color_out="black",
                   regions_color="0.5", regions_thickness=1, regions_style="--"
                  ):
    """

    Plots the profile likelihood of the given chains with respect to the given
    parameters, on a grid.

    The plotting region, unless stated otherwise (see keyword 'limits'), includes
    every point reached by the chain (in general a smaller region than that
    defined by the prior) plus a small padding (see keyword 'padding').

    A caveat: when mixing chains, the user is left to his/her own responsibility
    with respect to the correctness of the results.

    Main Parameters:
    ----------------
        
    chains: list of 'Chain' instances

    params: list of 2 parameter names

    labels: list of 2 str (optional)
        Alternative names of the parameters to show as axes labels.
        LaTeX typesetting is allowed.
        By default, the parameter names are used direcly.

    format: str (default: "-loglik")
        Quantity to plot in the color map:
            * '-loglik' for the log likelihood
            * 'delta-loglik' for the delta of the log likelihood with respect
                to the central value 'central_mloglik', that must be specified
            * 'chisq' for the chi squared, i.e., 2*(-loglik)
            * 'detachisq' same as 'delta-loglik' for the chi squared

    central_mloglik: float (optional)
        Central value of the -log(lik) in the color bar.
        If specified, values smaller than that one will be set to the inverse
        of the maximum loglik with respect to the central value.
        Example: value of the Planck baseline model loglik with the used
                 likelihoods (Planck+WMAPpol -> 4902.95).

    limits: list of 2 [min, max] floats (optional)
        To restrict the parameter ranges manually, give a list like
            limits = [[min_1, max_1], [min_2, max_2]]
        when 'min_i' ('max_i') is the maximum of the parameter 'i' = 1, 2.
        By default, the extrema of all the points in the chains are taken.

    n_grid: int (default: 100)
        Number of cells in the shortest side.

    aspect: float (defalut: 1 -- square plot)
        Aspect ratio of the plot: height/width.

    color_map: str (default: "jet_r")
        A valid 'matplotlib' colormap
        (see e.g. http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps ).

    black_and_white: bool (default: False)
        If True, select the "gray" color map. Overrides the "color_map" option.

    cb_orientation: "vertical" (default) of "horizontal"
        Orientation of the color bar. If horizontal, due to space constraints,
        the use of a custom ticks formatter may be needed: use the keyword

    bf_show: int (default: 1)
        If > 0, the best fit points are shown, either the overall (if 1) or that
        of each of the chains (if 2).

    regions_show: bool (default: True)
        If True, shows a rectangle (or the part of it within the plot limits)
        marking the border of the prior, if any.

    save_file: str (optional)
        If defined, instead of showing the plot, it is saved into the given file.


    Fine Tuninng Parameters:
    ------------------------

    fontsize_labels: float (default: 18)
        Font size of the axes and color bar labels.

    fontsize_ticks: float (default: 12)
        Font size of the ticks.

    padding: float (default 0.02)
        Fraction of padding around the plotted points.

    cb_ticks_formatter: str format (default: None)
        Formatter of the ticks of the color bars.
        Useful when horizontal orientation of the color bar is requested.

    cb_shrink: float (default: 1)
       Controls the size of the color bar.

    dpi: int (default: 150)
        Resolution used if the plot is saved to a file

    transparent: bool (default: True)
        Transparency of the background of the plot.

    not_yet: bool (default: False)
        If True, the output of the function is a string containing a
        "matplotlib.pyplot.savefig()" command.
        Useful if you want to add something else on top of the plot.

    bf_alpha=1, bf_radius=1, bf_thickness=1,
    bf_color_in="white", bf_color_out="black"
        Finely set the aspect of the best fit markers

    regions_color="0.5", regions_thickness=1, regions_style="--"
        Finely set the aspect of the prior ranges boxes.

    """
    # Make sense of input #####
    for chain in chains:
        assert isinstance(chain, Chain), (
            "The first argument must be a list of 'Chain' instances.")
        for i in [0, 1]:
            assert (params[i] in chain.varying_params() or 
                    params[i] in chain.derived_params()), (
                "The parameter %s is not on the chain %s."%(params[i], chain.name()))
    # Maxima and minima #####
    limits_new = [[None, None], [None, None]]
    if limits:
        for i in [0, 1]:
            if limits[i]:
                for j in [0, 1]:
                    limits_new[i][j] = limits[i][j]
    for i in [0, 1]:
        if not limits_new[i][0]:
            limits_new[i][0] = min([chain.get_min(params[i]) for chain in chains])
            limits_new[i][1] = max([chain.get_max(params[i]) for chain in chains])
    maxi = [limits_new[0][1], limits_new[1][1]]
    mini = [limits_new[0][0], limits_new[1][0]]
    assert maxi[0] > mini[0] and maxi[1] > mini[1],(
        "The given limits are not well formatted: min > max.")
    # Subdivisions #####
    n_grid = abs(int(n_grid))
    dims = [n_grid, n_grid]
    short_side = 0 if aspect <= 1 else 1
    dims[short_side] = int(dims[short_side]/float(aspect))
    # (notice the -1 in the "dims-1" in the next step)
    steps = [(maxi[i]-mini[i])/float(dims[i]) for i in [0, 1]]
    # Creating matrix with minima of minusloglik #####
    matrix = np.ones(shape=dims)
    indices = [(i,j) for i in range(dims[0])
                     for j in range(dims[1])]
    matrix = float("infinity") * matrix
    # Function to find the appropriate cell
    def to_matrix_element(p1value, p2value):
        i, j = (int(floor((p1value - mini[0]) / float(steps[0]))),
                int(floor((p2value - mini[1]) / float(steps[1]))))
        if i == dims[0]:
            i -= 1
        if j == dims[1]:
            j -= 1
        return i, j
    # Get the points into the matrix #####
    for chain in chains:
        for mloglik, p1value, p2value in zip(chain.points("mloglik"),
                                             chain.points(params[0]),
                                             chain.points(params[1])):
            i, j = to_matrix_element(p1value, p2value)
            matrix[i,j] = min(mloglik, matrix[i,j])
    # Prepare: center around central_mloglik, and NaN where there is no sample
    maxloglik = matrix.min()
    if central_mloglik:
        assert central_mloglik > maxloglik, (
            "The central -loglik value provided, %e, "%central_mloglik +
            "is smaller than the maximum -loglik, %e"%maxloglik)
        minloglik = central_mloglik - (maxloglik - central_mloglik)
        for i, j in indices:
            if matrix[i, j] < float("infinity"):
                matrix[i, j] = min(matrix[i, j], minloglik)
    else:
        minloglik = 0
        for i, j in indices:
            if matrix[i, j] != float("infinity"):
                minloglik = max(matrix[i, j], minloglik)
    cmap_range = 1+floor(minloglik), floor(maxloglik)
    for i, j in indices:
        if matrix[i, j] == float("infinity"):
            matrix[i, j] = float("nan")
    # Format #####
    factor = 2 if "chisq" in format else 1
    if "delta" in format:
        assert central_mloglik, ("A central log-likehood must be specified" +
                                 " if a 'delta'-like plot is requested.")
        delta = -1*central_mloglik
    else:
        delta = 0
    for i, j in indices:
        matrix[i, j] = factor*(matrix[i, j]+delta)
    # TODO: this variable is finally not used. Use it.
    cmap_range = [factor * (a+delta) for a in cmap_range]
    cmap_range[0] = int(floor(cmap_range[0]))
    cmap_range[1] = int(floor(cmap_range[1]))
    # Plot #####
    # The matrix must be transposed: the 0th component is the x axis
    matrix = matrix.transpose()
    fig = plt.figure()
    ax  = plt.axes()
    sq_aspect =  (maxi[0] - mini[0]) / (maxi[1] - mini[1])
    cmap = "gray" if black_and_white else color_map
    imsh = plt.imshow(matrix, cmap=cmap, interpolation="nearest", origin="lower",
                      aspect=aspect*sq_aspect, zorder=0,
                      extent = (mini[0], maxi[0], mini[1], maxi[1]))
    paddings = [padding*(limits_new[i][0]+limits_new[i][1]) for i in [0, 1]]
    paddings[1] *= aspect
    limits_plot = [[limits_new[0][0]-paddings[0], limits_new[0][1]+paddings[0]],
                   [limits_new[1][0]-paddings[1], limits_new[1][1]+paddings[1]]]
    ax.set_xlim(limits_plot[0][0], limits_plot[0][1])
    ax.set_ylim(limits_plot[1][0], limits_plot[1][1])
    # Color bar #####
    cb_options = {}
    cb_options["shrink"] = cb_shrink
    cb_options["orientation"] = cb_orientation
    if cb_orientation == "horizontal":
        cb_options["pad"] = 0.125
        cb_options["fraction"] = 0.05
    elif cb_orientation == "vertical":
        pass
    else:
        raise ValueError("The keyword 'cb_orientation' must be "+
                         "'horizontal' or 'vertical'")
    if cb_ticks_formatter:
        cb_options["format"]=cb_ticks_formatter
    cb = plt.colorbar(imsh, **cb_options)
    # Axes labels and ticks #####
    plt.tick_params(labelsize=fontsize_ticks)
    if labels:
        assert len(labels) == 2, "There must be 2 labels."
        plt.xlabel(labels[0], fontsize = fontsize_labels, fontweight = "bold")
        plt.ylabel(labels[1], fontsize = fontsize_labels, fontweight = "bold")
    cb.ax.tick_params(labelsize=fontsize_ticks)
    cb_label = "\ln\mathcal{L}$"
    if "delta" in format:
        cb_label = "\Delta"+cb_label
    if "chisq" in format:
        cb_label = "2\,"+cb_label
    cb_label = "$-"+cb_label
    cb.set_label(cb_label, fontsize = fontsize_labels, fontweight = "bold")
    # Show best fit markers #####
    if bf_show:
        bfs = {}
        for chain in chains:
            bfs[chain.name()] = [[b[chain.index_of_param(a)]
                for a in ["mloglik", params[0], params[1]]]
            for b in  chain.best_fit(more_than_one=True)]
        bf_plot = [item for sublist in bfs.values() for item in sublist]
        if bf_show == 1:
            bf_plot = [sorted(bf_plot, key=lambda x: x[0])[0]]
        # Plot them
        for bf in bf_plot:
            plt.scatter([bf[1]],[bf[2]], s=80*bf_radius, alpha=bf_alpha,
                        edgecolor=bf_color_out, facecolor=bf_color_in,
                        linewidths=1.5*bf_thickness)
            print "Best fit: % .6e ; %s = % .6e  ; %s = % .6e"%(
                factor*(bf[0]+delta), params[0], bf[1], params[1], bf[2])
    # Show chain priors limits #####
    if regions_show:
        for chain in chains:
            limits = [[chain.get_limits(params[i])[j] for j in [0, 1]] for i in [0, 1]]
            # Left
            if limits[0][0] >= limits_plot[0][0]:
                plt.plot([limits[0][0], limits[0][0]],
                         [max(limits[1][0], limits_plot[1][0]),
                          min(limits[1][1], limits_plot[1][1])],
                         color=regions_color, linewidth=2*regions_thickness,
                         linestyle=regions_style, zorder=1)
            # Right
            if limits[0][1] <= limits_plot[0][1]:
                plt.plot([limits[0][1], limits[0][1]],
                         [max(limits[1][0], limits_plot[1][0]),
                          min(limits[1][1], limits_plot[1][1])],
                         color=regions_color, linewidth=2*regions_thickness,
                         linestyle=regions_style, zorder=1)
            # Bottom
            if limits[1][0] >= limits_plot[1][0]:
                plt.plot([max(limits[0][0], limits_plot[0][0]),
                          min(limits[0][1], limits_plot[0][1])],
                         [limits[1][0], limits[1][0]],
                         color=regions_color, linewidth=2*regions_thickness,
                         linestyle=regions_style, zorder=1)
            # Top
            if limits[1][1] <= limits_plot[1][1]:
                plt.plot([max(limits[0][0], limits_plot[0][0]),
                          min(limits[0][1], limits_plot[0][1])],
                         [limits[1][1], limits[1][1]],
                         color=regions_color, linewidth=2*regions_thickness,
                         linestyle=regions_style, zorder=1)
    # Plotting #####
    if not save_file:
        plt.show()
        return
    else:
        plotting_command = ("plt.savefig('%s', transparent = %s, dpi = %s, "%(
                             save_file, transparent, int(dpi))+
                            "bbox_inches='tight', pad_inches=0.1)")
        if not_yet:
            return plotting_command
        else:
            eval(plotting_command)
            plt.close()    
            return
