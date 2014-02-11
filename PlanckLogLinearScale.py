###############################################################################
# Log+Linear scale used to plot the CMB power spectrum by the ESA Planck team #
###############################################################################

# Note: import the whole module, not just the class,
#       since the scale must be registered (see last line)

import numpy as np
from matplotlib import scale as mscale
from matplotlib.transforms import Transform
from matplotlib.ticker import Formatter, FixedLocator

class PlanckLogLinearScale(mscale.ScaleBase):
    # The scale class must have a member ``name`` that defines the
    # string used to select the scale.  For example,
    # ``gca().set_yscale("mercator")`` would be used to select this
    # scale.
    name = 'planck'
    def __init__(self, axis, **kwargs):
        """
        Any keyword arguments passed to ``set_xscale`` and
        ``set_yscale`` will be passed along to the scale's
        constructor.
        """
        mscale.ScaleBase.__init__(self)
    def get_transform(self):
        """
        Override this method to return a new instance that does the
        actual transformation of the data.
        """
        return self.PlanckLogLinearTransform()
    def set_default_locators_and_formatters(self, axis):
         """
         Override to set up the locators and formatters to use with the
         scale.  This is only required if the scale requires custom
         locators and formatters.  Writing custom locators and
         formatters is rather outside the scope of this example, but
         there are many helpful examples in ``ticker.py``.

#         In our case, the Mercator example uses a fixed locator from
#         -90 to 90 degrees and a custom formatter class to put convert
#         the radians to degrees and put a degree symbol after the
#         value::
         """
# EXAMPLE CODE
#         class DegreeFormatter(Formatter):
#             def __call__(self, x, pos=None):
#                 # \u00b0 : degree symbol
#                 return "%d\u00b0" % ((x / np.pi) * 180.0)
#         deg2rad = np.pi / 180.0
#         axis.set_major_locator(FixedLocator(np.arange(1, 2500, 10)))
#         axis.set_major_formatter(DegreeFormatter())
#         axis.set_minor_formatter(DegreeFormatter())
         change = 50
         axis.set_major_locator(FixedLocator(np.concatenate((np.array([2,5,10,20]),
                                                             np.array([change]),
                                                             np.arange(500, 2500, 500)))))
         axis.set_minor_locator(FixedLocator(np.concatenate((np.arange(2,10),
                                                             np.arange(10,50,10),
                                                             np.arange(floor(change/100), 2500, 100)))))
    class PlanckLogLinearTransform(Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True
        has_inverse = True
        def __init__(self):
            Transform.__init__(self)
        def transform_non_affine(self, a):
            """
            This transform takes an Nx1 ``numpy`` array and returns a
            transformed copy.  Since the range of the Mercator scale
            is limited by the user-specified threshold, the input
            array must be masked to contain only valid values.
            ``matplotlib`` will handle masked arrays and remove the
            out-of-range data from the plot.  Importantly, the
            ``transform`` method *must* return an array that is the
            same shape as the input array, since these values need to
            remain synchronized with values in the other dimension.
            """
            change = 50.
            factor = 750.
            lower   = a[np.where(a<=change)]
            greater = a[np.where(a> change)]
            if lower.size:
                lower_after   = factor*np.log10(lower)
            if greater.size:
                greater_after = factor*np.log10(change) + (greater-change)
            # Only low
            if not(greater.size):
                return lower_after
            # Only high
            if not(lower.size):
                return greater_after
            return np.concatenate((lower_after, greater_after))
        def inverted(self):
            """
            Override this method so matplotlib knows how to get the
            inverse transform for this transform.
            """
            return InvertedPlanckLogLinearTransform()
    class InvertedPlanckLogLinearTransform(Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True
        has_inverse = True
        def transform_non_affine(self, a):
            change = 50.
            factor = 750.
            lower   = a[np.where(a<=factor*np.log10(change))]
            greater = a[np.where(a> factor*np.log10(change))]
            if lower.size:
                lower_after   = np.power(10.0, lower/factor)
            if greater.size:
                greater_after = greater + change - factor*np.log10(change)
            # Only low
            if not(greater.size):
                return lower_after
            # Only high
            if not(lower.size):
                return greater_after
            return np.concatenate((lower_after, greater_after))
        def inverted(self):
            return PlanckLogLinearTransform()

# Finished. Register the scale!
mscale.register_scale(PlanckLogLinearScale)
