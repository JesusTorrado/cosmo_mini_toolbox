##############################################
# Tool to plot and compare CMB power spectra #
##############################################

import os
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle

# Local
from CMBspectrum import CMBspectrum
import PlanckLogLinearScale

nonpos = "mask"

def plot_Cl_CMB(CMB_spectra,
                # Main customisation parameters
                Deltas=1, pol="TT", l_max=None, data_points=True, Deltas_max=175,
                title=None, scale = "planck", aspect="auto", save_file=None,
                black_and_white=False,
                # Fine tuning parameters
                lensed=True, l_prefactor=True,
                ticks_fontsize=10, labels_fontsize=14, title_fontsize=14,
                transparent=False, transparent_frame=False,
                dpi=150, not_yet=False,
               ):
    """
    A tool to plot (absolute) comparisons between different CMB spectra.
    The first one on the list is taken as a reference.


    Mandatory arguments:
    --------------------

    CMB_spectra: list of 'CMBspectrum' instances.
        The first one is taken as the reference one, if more than one is given.
        If a name has been defined for them (see doc. of the class) it is used
        as a label in the legend.

    
    Main customisation parameters:
    ------------------------------

    Deltas: one of [-1, 0, 1 (default)]
        *  1: plot both the spectra and the differences
        *  0: plot only the full spectra
        * -1: plot only the differences

    pol: one of ["TT" (default), "TE", "EE", "BB", "phiphi", "Tphi", "Ephi"]
        Poarisation (or lensing) spectra to plot.

    l_max: int (default: None)
        Maximum multipole to plot.
        If not defined, the maximum common l maximum is used.

    data_points: bool (default: True)
        Whether to plot the TT power spectrum data points of Planck.
        N.B.: They are just for illustration purposes: they are strongly correlated.

    Deltas_max: float (default: 175)
        Range of the y axis in the differences plot.

    title: str (default: None)
        Title of the plot

    scale: one of ["linear", "log", "planck" (default)]
        Scale for the multipoles.
        "planck" is the scale used in many public plots of the ESA Planck CMB
        spectrum: logarithmic in base 10 up to l=50, and linear from there on.

    aspect: matplotlib aspect value (default: "auto")
        Aspect ratio of the plot.
        Used only when 'Deltas'=[0, -1], ignored otherwise.

    black_and_white: bool (default: False)
        If True, select the "gray" color map. Overrides the "color_map" option.

    save_file: str (default: None)
        If defined, instead of showing the plot, it is saved into the given
        file name.


    Fine tuning parameters:
    -----------------------

    lensed: bool (default: True)
        Whether to used lensed or unlensed (if False) power spectra.

    l_prefactor: bool (default: True)
        Whether to plot the power spectrum multiplied by the l(l+1)/(2pi)
        prefactor.

    ticks_fontsize=10, labels_fontsize=14, title_fontsize=14,
        Fontsize of the axes ticks, labels and the plot title

    transparent: bool (default: False)
        Transparency of the background of the plot.

    transparent_frame: bool (default: False)
        Transparency of the frame of the plot.
        If 'transparent=True', this one is set to True too.

    dpi: int (default: 150)
        Resolution used if the plot is saved to a file

    not_yet: bool (default: False)
        If True, the output of the function is a string containing a
        "matplotlib.pyplot.savefig()" command.
        Useful if you want to add something else on top of the plot.

    """
    # Tests on the input ####
    assert all(lambda s: isinstance(s, CMBspectrum) for s in CMB_spectra), \
        "Some of the spectra provided are not instances of 'CMBspectrum'."
    # Prepare data for the plot #####
    spectrum_types = [a.lower() for a in ["TT", "TE", "EE", "BB",
                                         "phiphi", "Tphi", "Ephi"]]
    assert pol.lower() in spectrum_types, \
        ("Spectrum type (or polarisation) not recognised. Choose one of " + 
         str(spectra_types))
    if pol.lower() != "tt":
        data_points = False
    l, Cl = {}, {}
    for spectrum in CMB_spectra:
        if lensed:
            l [spectrum.name()] = spectrum.ll()
            Cl[spectrum.name()] = spectrum.lCl(pol, units="muK",
                                               l_prefactor=l_prefactor)
        else:
            l [spectrum.name()] = spectrum.ul()
            Cl[spectrum.name()] = spectrum.uCl(pol, units="muK",
                                               l_prefactor=l_prefactor)
    l_maxes = [l[spectrum.name()][-1] for spectrum in CMB_spectra]
    if l_max:
        assert all([m>=l_max for m in l_maxes]), \
            "Some of the spectra provided do not reach the l_max provided."
    else:
        for spectrum in CMB_spectra:
            l_max = min(l_maxes)
            i= np.where(l[spectrum.name()]==l_max)[0][0]
            l [spectrum.name()] = l [spectrum.name()][:i+1]
            Cl[spectrum.name()] = Cl[spectrum.name()][:i+1]
    # Prepare data points
    if data_points:
        data_folder = os.path.join(os.path.split(__file__)[0], "../data")
        data_lowl  = np.loadtxt(os.path.join(data_folder,
                                             "planck_spectrum_lowl.txt"))
        data_lowl = data_lowl.transpose()
        data_highl = np.loadtxt(os.path.join(data_folder,
                                             "planck_spectrum_highl.txt"))
        data_highl = data_highl.transpose()
    # Functions to prepare the spectra #####
    # Upper
    def plot_Cl(axes):
        for i, spectrum in enumerate(CMB_spectra):
            axes.plot(l [spectrum.name()], Cl[spectrum.name()],
                      color=next(colour_cycler), linestyle=next(style_cycler),
                      label = spectrum.name(), zorder = i)
        if data_points :
            axes.errorbar(data_lowl[0], data_lowl[1],
                          yerr = [data_lowl[3], data_lowl[2]],
                          fmt = ".", color = colour_data_lowl,  zorder = -2)
            axes.errorbar(data_highl[0], data_highl[1], yerr = data_highl[2],
                          fmt = ".", color = colour_data_highl, zorder = -1)
    # Lower
    def compare(name1, name2, invert = False) :
        diffs = []
        for l_i, cl_i in zip(l[name2], Cl[name2]) :
            i= np.where(l[name1]==l_i)[0][0]
            diffs.append((cl_i-Cl[name1][i])*(-1 if invert else 1))
        return l[name2], diffs
    def compare_data(name1, l2, Cl2, invert = False) :
        diffs = []
        for l_i, cl_i in zip(l2, Cl2) :
            i= np.where(l[name1]==l_i)[0][0]
            diffs.append((cl_i-Cl[name1][i])*(-1 if invert else 1))
        return l2, diffs
    def plot_Deltas(axes):
        # First one
        axes.plot([0 for a in l[CMB_spectra[0].name()]],
                   color=next(colour_cycler), linestyle=next(style_cycler),
                  label = CMB_spectra[0].name(), zorder = 0)
        # Rest
        for i, spectrum in enumerate(CMB_spectra[1:]):
            l_cmp, cl_cmp = compare(CMB_spectra[0].name(), spectrum.name())
            axes.plot(l_cmp, cl_cmp,
                      color=next(colour_cycler), linestyle=next(style_cycler),
                      label=spectrum.name(), zorder = i+1)

        # Data points
        if data_points :
            l_cmp_data, cmp_data = compare_data(CMB_spectra[0].name(),
                                                data_lowl[0], data_lowl[1])
            axes.errorbar(l_cmp_data, cmp_data,
                          yerr = [data_lowl[3], data_lowl[2]],
                          fmt = ".", color = colour_data_lowl, zorder = -2)
            l_cmp_data, cmp_data = compare_data(CMB_spectra[0].name(),
                                                data_highl[0], data_highl[1])
            axes.errorbar(l_cmp_data, cmp_data, yerr = data_highl[2],
                          fmt = ".", color = colour_data_highl, zorder = -1)    
    # Prepare the plot #####
    # Arrange
    fig, axarr = plt.subplots(2, sharex=(Deltas==1))
    ax_Cl     = axarr[0]
    ax_Deltas = axarr[1]
    if Deltas == 0:
        fig.delaxes(ax_Deltas)
        ax_Deltas = None
    elif Deltas == -1:
        fig.delaxes(ax_Cl)
        ax_Cl = None
    axes_list = [ax for ax in [ax_Cl, ax_Deltas] if ax]
    # Colours and styles
    list_of_colours = ["blue", "red", "green", "orange",
                       "purple", "darkcyan", "brown"]
    list_of_styles  = ["-"]
    if black_and_white:
        list_of_colours = ["black"]
        list_of_styles  = ["-", ":", "--", "-."]
        if len(CMB_spectra) > len(list_of_styles):
            print ("WARNING: if more that %d spectra "%len(list_of_styles) +
                  "are plot in balck and white, the line styles will be repeated.")
    colour_data_lowl  = "0.50"
    colour_data_highl = "0.30"
    # Do plot
    if ax_Cl:
        # Generate the cycles
        colour_cycler = cycle(list_of_colours)
        style_cycler  = cycle(list_of_styles)
        plot_Cl(ax_Cl)
    if ax_Deltas:
        # (Re)generate the cycles
        colour_cycler = cycle(list_of_colours)
        style_cycler  = cycle(list_of_styles)
        plot_Deltas(ax_Deltas)
    # Formatting #####
    # Axes labels
    if ax_Cl:
        ax_Cl.set_ylabel((r"$\frac{\ell(\ell+1)}{2\pi}\,C_\ell^{"+
                          r"\mathrm{"+pol+r"}}$ $(\mu\mathrm{K}^2)$"),
                         fontsize=labels_fontsize)
        if not(ax_Deltas):
            ax_Cl.set_xlabel(r"Multipole, $\ell$", fontsize=labels_fontsize)
    if ax_Deltas:
        ax_Deltas.set_ylabel((r"$\frac{\ell(\ell+1)}{2\pi}\,\Delta C_\ell^{"+
                              r"\mathrm{"+pol+r"}}$ $(\mu\mathrm{K}^2)$"),
                         fontsize=labels_fontsize)
        ax_Deltas.set_xlabel(r"$\mathrm{Multipole,}\,\ell$",
                             fontsize=labels_fontsize)
    # Legend
    if ax_Cl:
        leg = ax_Cl.legend(loc = (2 if scale == "log" else 1),
                           fancybox=True, prop={'size':10})
    else:
        leg = ax_Deltas.legend(loc = 1, fancybox=True, prop={'size':10})
    # Title
    if title :
        ax = ax_Cl if ax_Cl else ax_Deltas
        ax.set_title(title, horizontalalignment='center',
                     y=1.1, fontsize=title_fontsize)
    # Limits and ticks
    # If 2 plots
    if ax_Cl and ax_Deltas:
        ax_Cl.xaxis.set_ticks_position('top')
        fig.subplots_adjust(hspace = 0.001)
    for ax in axes_list:
        ax.set_xlim(2, l_max)
        plt.setp(ax.get_yticklabels(), fontsize=ticks_fontsize)
        plt.setp(ax.get_xticklabels(), fontsize=ticks_fontsize)
    if ax_Deltas:
        ax_Deltas.set_ylim(-1*Deltas_max, Deltas_max)
    # Axes' scale
    scales = ["linear", "log", "planck"]
    assert scale.lower() in scales, \
        "'scale' must be in "+str(scales)
    for ax in axes_list:
        ax.set_xscale(scale)
    # Aspect
    if Deltas != 1:
        for ax in axes_list:
            try:
                ymin, ymax = ax.get_ylim()
                xmin, xmax = ax.get_xlim()
# Attemp to modify aspect independently
#                print ax.get_aspect()
#                ax.set_aspect(0.25, adjustable = "box", anchor = "C")
#                aspect = float(aspect)
#                ax.set_aspect(aspect*float(xmax-xmin)/float(ymax-ymin))
            except ValueError:
                ax.set_aspect(aspect)
    # Transparency
    fig.frameon = not(transparent_frame or transparent)
    # Legend
    leg.get_frame().set_alpha(float(0.8))
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
