#################################################################
# Tools to plot various numerical approximations to likelihoods #
#################################################################

import numpy as np
import matplotlib.pyplot as plt
from math import floor

# Local import
from Chain import Chain

### Plot of 2D likelihoods
def plot_lik_2D(mode, chains, params,
                # Main customisation parameters
                labels=None, format="-loglik", central_mloglik=None,
                limits=None, n_grid=100, aspect=1,
                color_map="jet_r", black_and_white=False,
                cb_orientation="vertical",
                bf_show=1, regions_show=True, save=True, axes=None,
                # Fine tuning
                fontsize_labels=18, fontsize_ticks=12,
                cb_ticks_formatter=None, cb_shrink=float(1),
                padding = 0.02, dpi=150,
                transparent=False, transparent_frame=False,
                bf_alpha=1, bf_radius=1, bf_thickness=1,
                bf_color_in="white", bf_color_out="black",
                regions_color="0.5", regions_thickness=1, regions_style="--",
                ):
    """
    Plots the [marginal|mean|profile] likelihood of the given chains
    with respect to the given parameters, on a grid, meaning

    * Mean:      (sum_i #_i * -loglik_i) / (sum_i #_i)
    * Marginal:  log(e * sum_i #_i)
    * Profile:   max(-loglik_i)

    being the sums over 'i' extended to all chain points falling within a given
    cell, and being '#_i' the number of stops of the chain point 'i'.

    The plotting region, unless stated otherwise (see keyword 'limits'), includes
    every point reached by the chain (in general a smaller region than that
    defined by the prior) plus a small padding (see keyword 'padding').

    A caveat: when mixing chains, the user is left to his/her own responsibility
    with respect to the correctness of the results.


    Mandatory arguments:
    --------------------
        
    mode: one of ["marginal", "mean", "profile"]

    chains: list of 'Chain' instances

    params: list of 2 parameter names


    Main customisation parameters:
    ------------------------------

    labels: list of 2 str (default: None)
        Alternative names of the parameters to show as axes labels.
        LaTeX typesetting is allowed.
        By default, the parameter names are used direcly.

    format: str (default: "-loglik")
        Quantity to plot in the color map ("mean" and "profile")
        and to print for the best fits:
            * '-loglik' for the log likelihood
            * 'delta-loglik' for the delta of the log likelihood with respect
                to the central value 'central_mloglik', that must be specified
            * 'chisq' for the chi squared, i.e., 2*(-loglik)
            * 'detachisq' same as 'delta-loglik' for the chi squared

    central_mloglik: float (default: None)
        Central value of the -log(lik) in the color bar.
        If specified, values smaller than that one will be set to the inverse
        of the maximum loglik with respect to the central value.
        Example: value of the Planck baseline model loglik with the used
                 likelihoods (Planck+WMAPpol -> 4902.95).
        ("mean" and "profile" only)

    limits: list of 2 [min, max] floats (default: None)
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

    save: bool or str (default: True, i.e. 'show()')
        This keyword defines what to do with the resulting plot. Three outcomes
        are possible:
        * "path/name.extension": saves the plot to the given file name.
        * True: shows the plot.
        * False: does nothing, and returns a tuple containing the axes instance
            and a keyword dictionary that can be passed to 'pyplot.savefig()'.
            This can be useful if you want, e.g. to plot something on top of it.

    axes: matplotlib.axes.Axes (default: None)
        Allows to specify the axes in which the figure must be plotted.
        Useful for including the plot as a subplot in a bigger figure.

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

    transparent: bool (default: False)
        Transparency of the background of the plot.

    transparent_frame: bool (default: False)
        Transparency of the frame of the plot.
        If 'transparent=True', this one is set to True too.

    bf_alpha=1, bf_radius=1, bf_thickness=1,
    bf_color_in="white", bf_color_out="black"
        Finely set the aspect of the best fit markers

    regions_color="0.5", regions_thickness=1, regions_style="--"
        Finely set the aspect of the prior ranges boxes.

    """
    # Make sense of input #####
    if isinstance(chains, Chain):
        chains = [chains]
    for chain in chains:
        assert isinstance(chain, Chain), (
            "The first argument must be a list of 'Chain' instances.")
        for i in [0, 1]:
            assert (params[i] in chain.varying_parameters() or 
                    params[i] in chain.derived_parameters()), (
                "The parameter %s is not on the chain %s."%(params[i], chain.name()))
    # Format of the color scale (profile and mean) and the best fit (all) #####
    factor = 2 if "chisq" in format else 1
    if "delta" in format:
        assert central_mloglik, ("A central log-likehood must be specified" +
                                 " if a 'delta'-like plot is requested.")
        delta = -1*central_mloglik
    else:
        delta = 0
    # Maxima and minima #####
    limits_new = [[None, None], [None, None]]
    if limits:
        for i in [0, 1]:
            if limits[i]:
                for j in [0, 1]:
                    limits_new[i][j] = limits[i][j]
    for i in [0, 1]:
        if limits_new[i][0] is None:
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
    # Creating matrix #####
    matrix = np.ones(shape=dims)
    indices = [(i,j) for i in range(dims[0])
                     for j in range(dims[1])]
    if mode == "profile":
        matrix = float("infinity") * matrix
    elif mode == "mean":
        matrix = 0*matrix
        matrix_num = matrix.copy()
    elif mode == "marginal":
        matrix = 0 * matrix
    # Function to find the appropriate cell
    def to_matrix_element(p1value, p2value):
        i, j = (int(floor((p1value - mini[0]) / float(steps[0]))),
                int(floor((p2value - mini[1]) / float(steps[1]))))
        # Handling extremes
        if i == dims[0]:
            i -= 1
        if j == dims[1]:
            j -= 1
        return i, j
    # Get the points into the matrix #####
    for chain in chains:
        for num, mloglik, p1value, p2value in zip(chain.points("#"),
                                                  chain.points("mloglik"),
                                                  chain.points(params[0]),
                                                  chain.points(params[1])):
            if (p1value < mini[0] or p1value > maxi[0] or
                p2value < mini[1] or p2value > maxi[1]):
                continue
            i, j = to_matrix_element(p1value, p2value)
            if mode == "profile":
                matrix[i,j] = min(mloglik, matrix[i,j])
            elif mode == "mean":
                matrix[i,j] += mloglik*num
                matrix_num[i,j] += num
            elif mode == "marginal":
                matrix[i,j] += num
    # Log and clipping of the marginal matrix
    if mode == "marginal":
        matrix = np.log(np.e*matrix)
        maxlogsteps = matrix.max()
        matrix = matrix.clip(-maxlogsteps, maxlogsteps)
    # Normalising the "mean" matrix:
    # (and filling with infinities to treat is as the profile one)
    if mode == "mean":
        for i, j in indices:
            if matrix_num[i, j] == 0:
                matrix[i, j] = float("infinity")
            else:
                matrix[i, j] /= matrix_num[i, j]
    # Centering and reducing the range  -- infinity to NaN
    if mode in ["profile", "mean"]:
        if central_mloglik:
            maxloglik = matrix.min()
            assert central_mloglik > maxloglik, (
                "The central -loglik value provided, %e, "%central_mloglik +
                "is smaller than the maximum -loglik, %e"%maxloglik)
            minloglik = central_mloglik - (maxloglik - central_mloglik)
        else:
            minloglik = float("infinity")
        for i, j in indices:
            if matrix[i, j] == float("infinity"):
                matrix[i, j] = float("nan")
            else:
                matrix[i, j] = min(matrix[i, j], minloglik)
                if "delta" in format:
                    matrix[i, j] -= central_mloglik
                if "chisq" in format:
                    matrix[i, j] *= 2
    # Plot #####
    # The matrix must be transposed: the 0th component is the x axis
    matrix = matrix.transpose()
    # Create axes, in none given
    if not(axes):
        fig = plt.figure()
        fig.frameon = not(transparent_frame or transparent)
        axes  = plt.axes()
    sq_aspect =  (maxi[0] - mini[0]) / (maxi[1] - mini[1])
    if mode == "marginal":
        if color_map[-2:] == "_r":
            color_map = color_map[:-2]
        else:
            color_map += "_r"
    cmap = "gray" if black_and_white else color_map
    imsh = axes.imshow(matrix, cmap=cmap, interpolation="nearest", origin="lower",
                      aspect=aspect*sq_aspect, zorder=0,
                      extent = (mini[0], maxi[0], mini[1], maxi[1]))
    paddings = [padding*abs(limits_new[i][1]-limits_new[i][0]) for i in [0, 1]]
    paddings[short_side] = paddings[short_side]*float(aspect)
    limits_plot = [[limits_new[0][0]-paddings[0], limits_new[0][1]+paddings[0]],
                   [limits_new[1][0]-paddings[1], limits_new[1][1]+paddings[1]]]
    axes.set_xlim(limits_plot[0][0], limits_plot[0][1])
    axes.set_ylim(limits_plot[1][0], limits_plot[1][1])
    # Axes labels and ticks #####
    axes.tick_params(labelsize=fontsize_ticks)
    if not labels:
        labels = params
    assert len(labels) == 2, "There must be 2 labels."
    axes.set_xlabel(labels[0], fontsize = fontsize_labels, fontweight = "bold")
    axes.set_ylabel(labels[1], fontsize = fontsize_labels, fontweight = "bold")
    # Color bar #####
    cb_options = {}
    assert cb_orientation in ["horizontal", "vertical"], (
        "The keyword 'cb_orientation' must be 'horizontal' or 'vertical'")
    cb_options["orientation"] = cb_orientation
    if cb_orientation == "horizontal":
        cb_options["pad"] = 0.125
        cb_options["fraction"] = 0.05
        cb_options["shrink"] = 1
    elif cb_orientation == "vertical":
        cb_options["shrink"] = aspect+2*padding
    cb_options["shrink"] *= cb_shrink
    if cb_ticks_formatter:
        cb_options["format"]=cb_ticks_formatter
    cb = plt.colorbar(imsh, ax=axes, **cb_options)
    cb.ax.tick_params(labelsize=fontsize_ticks)
    if mode in ["profile", "mean"]:
        cb_label = "\ln\mathcal{L}"
        if "delta" in format:
            cb_label = "\Delta"+cb_label
        if "chisq" in format:
            cb_label = "2\,"+cb_label
        cb_label = r"$-%s$"%cb_label
    if mode == "marginal":
        cb_label = r"$\propto\log\left(\#\mathrm{steps}\right)$"
    cb.set_label(cb_label, fontsize = fontsize_labels, fontweight = "bold")
    # Show best fit markers #####
# TODO: Fix this!
#    if bf_show:
#        bfs = {}
#        for chain in chains:
#            bfs[chain.name()] = [[b[chain.index_of_param(a)]
#                for a in ["mloglik", params[0], params[1]]]
#            for b in  chain.best_fit(more_than_one=True)]
#        bf_plot = [item for sublist in bfs.values() for item in sublist]
#        if bf_show == 1:
#            bf_plot = [sorted(bf_plot, key=lambda x: x[0])[0]]
#        # Plot them
#        for bf in bf_plot:
#            axes.scatter([bf[1]],[bf[2]], s=80*bf_radius, alpha=bf_alpha,
#                        edgecolor=bf_color_out, facecolor=bf_color_in,
#                        linewidths=1.5*bf_thickness)
#            print "Best fit: % .6e ; %s = % .6e  ; %s = % .6e"%(
#                factor*(bf[0]+delta), params[0], bf[1], params[1], bf[2])
    # Show chain priors limits #####
    if regions_show:
        for chain in chains:
            limits = [[chain.get_limits(params[i])[j] for j in [0, 1]] for i in [0, 1]]
            # Left
            if limits[0][0] >= limits_plot[0][0]:
                axes.plot([limits[0][0], limits[0][0]],
                         [max(limits[1][0], limits_plot[1][0]),
                          min(limits[1][1], limits_plot[1][1])],
                         color=regions_color, linewidth=2*regions_thickness,
                         linestyle=regions_style, zorder=1)
            # Right
            if limits[0][1] <= limits_plot[0][1]:
                axes.plot([limits[0][1], limits[0][1]],
                         [max(limits[1][0], limits_plot[1][0]),
                          min(limits[1][1], limits_plot[1][1])],
                         color=regions_color, linewidth=2*regions_thickness,
                         linestyle=regions_style, zorder=1)
            # Bottom
            if limits[1][0] >= limits_plot[1][0]:
                axes.plot([max(limits[0][0], limits_plot[0][0]),
                          min(limits[0][1], limits_plot[0][1])],
                         [limits[1][0], limits[1][0]],
                         color=regions_color, linewidth=2*regions_thickness,
                         linestyle=regions_style, zorder=1)
            # Top
            if limits[1][1] <= limits_plot[1][1]:
                axes.plot([max(limits[0][0], limits_plot[0][0]),
                          min(limits[0][1], limits_plot[0][1])],
                         [limits[1][1], limits[1][1]],
                         color=regions_color, linewidth=2*regions_thickness,
                         linestyle=regions_style, zorder=1)
    # Ticks #####
    from matplotlib.ticker import AutoMinorLocator
    axes.xaxis.set_minor_locator(AutoMinorLocator(10))
    axes.yaxis.set_minor_locator(AutoMinorLocator(10))
    # Plotting #####   
    # Options
    options = {"dpi": int(dpi), "transparent": transparent,
               "bbox_inches": "tight", "pad_inches": 0.1}
    # Save
    if isinstance(save, basestring):
        plt.savefig(save, **options)
        plt.close()    
    # Show
    elif save:
        plt.show()
        plt.close()
    # Return
    else:
        return axes, options
