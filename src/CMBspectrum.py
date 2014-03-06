import os as os
import numpy as np

T_CMB = 2.726

class CMBspectrum():
    """
    Stores a CMB power spectrum from a CLASS or CAMB file.
    The spectrum is internally stored dimensionless and with the l(l+1)/(2pi) prefactor.

    Parameters
    ----------
    folder   : str
               Name of the folder containing the power spectrum.
    prefix   : str
               If the spectrum names do have a prefix, i.e. "test_", give it here.
               Otherwise, CLASS and CAMB standard names are assumed, i.e. "folder/cl.dat".
               In the case of CAMB, the prefix is to be given *without* the final "_".
    name     : str, optional
               Name given to the power spectrum for plotting purposes.
               If not specified, file name without extension is used.
    code     : str
               'CLASS' or 'CAMB' (any capitalisation).
    """
    def __init__(self, folder, prefix=None, name=None, code="CLASS"):
        if not prefix:
            prefix = ""
        if code.lower() == "class":
            self._code = "class"
        if code.lower() == "camb":
            self._code = "camb"
        self._load_data(folder, prefix)
        if name:
            self._name = name
        else:
            if prefix:
                for i in ["_", "-"]:
                    if prefix[-1] == i:
                        prefix.rstrip(i)
                self._name = prefix
            else:
                self._name = folder.rstrip("/").split("/")[-1]

    ### Load spectra: CLASS & CAMB
    def _load_data(self, folder, prefix):
        if self._code == "class":
            pname = os.path.join(folder, prefix + "parameters.ini")
            #uname = os.path.join(folder, prefix + "unused_parameters")
            cname = os.path.join(folder, prefix + "cl.dat")
            lname = os.path.join(folder, prefix + "cl_lensed.dat")
        if self._code == "camb":
            pname = os.path.join(folder, prefix + "_params.ini")
            cname = os.path.join(folder, prefix + "_scalCls.dat")
            lname = os.path.join(folder, prefix + "_lensedCls.dat")
        # 1. Parameters
        try:
            pfile = open(pname, "r")
        except IOError:
            raise IOError("The parameters file does not exist: '%s'"%pname)
        self._parameters = {}
        for line in pfile:
            if "=" in line and line.lstrip()[0] != "#":
                left, right = [a.strip() for a in line.split("=")[0:2]]
                if right:
                    self._parameters[left] = right
        pfile.close()
        # 2. Spectrum
        self._lensed = False
        if self._code == "class":
            columns_all = ["l", "TT", "EE", "TE", "BB", "phiphi", "Tphi", "Ephi"]
            if not "tCl" in self.parameter("output"):
                raise NotImplementedError("Not implemented when not using 'tCl' in 'output'.")
            if not "pCl" in self.parameter("output"):
                columns_all = [a for a in columns_all if not "E" in a and not "B" in a]
            if not "lCl" in self.parameter("output"):
                columns_all = [a for a in columns_all if not "phi" in a]
            else:
                if "y" in self.parameter("lensing"):
                    self._lensed = True
        if self._code == "camb":
            columns_all = ["l", "TT", "EE", "BB", "TE", "phiphi", "Tphi", "Ephi"]
            if self.parameter("do_lensing") != "T":
                columns_all = [a for a in columns_all if not "phi" in a]
            else:
                self._lensed = True                
        self._columns = columns_all
        self._columns_indices = dict([a,i] for i,a in enumerate(self._columns))
        # 2.1. unlensed
        try:
            data = np.loadtxt(cname)
        except IOError:
            raise IOError("The spectrum file does not exist: '%s'"%cname)
        self._uCl = np.transpose(data)
        # 2.2. lensed
        if self._lensed:
            try:
                data = np.loadtxt(lname)
            except IOError:
                raise IOError("The lensed spectrum file does not exist: '%s'"%cname)
            self._lCl = np.transpose(data)
        # 3. Units
        if self._code == "camb":
            if self.parameter("CMB_outputscale"):
                scale = float(self.parameter("CMB_outputscale"))
            else:
                scale = 7.4311e12 # CAMB default
            self._uCl[1:] /= scale
            self._lCl[1:] /= scale
        self._l_prefactor = True
        self._units = "1"

    ### Set and Retrieve name
    def set_name(self, name):
        self._name = name
    def name(self):
        return self._name

    ### Retrieve Cosmological Code parameters
    def parameter(self, name):
        return self._parameters.get(name)
    def param(self, name):
        return self.parameter(name)
    def parameters(self):
        return self._parameters
    def params(self):
        return self.parameters()

    ### Retrieve spectrum, unlensed
    def uCl(self, column, units="1", l_prefactor=True, T_CMB = T_CMB):
        # l-prefactor
        if ((self._l_prefactor and l_prefactor) or
            (not self._l_prefactor and not l_prefactor)):
            factor = 1
        elif self._l_prefactor and not l_prefactor:
            factor = 1/(self.ul()*(self.ul()+1)/(2.*np.pi))
        elif not self._l_prefactor and l_prefactor:
            factor = self.ul()*(self.ul()+1)/(2.*np.pi)
        try:
            return (factor * self._units_from_to(self._units, units, T_CMB) * 
                    self._uCl[self._columns_indices[column]])
        except KeyError:
            raise KeyError("It seems the '%s' spectrum "%column+
                           "has not been calculated. "+
                           "The calcluated spectra are "+str(self._columns[1:]))
    def ul(self):
        return self.uCl("l")
    def ulmax(self):
        return self.uCl("l")[-1]

    ### Retrieve spectrum, lensed
    def lCl(self, column, units="1", l_prefactor=True, T_CMB = T_CMB):
        assert  self._lensed, "No lensed spectrum has been calculated."
        if column=="l":
            return self._lCl[self._columns_indices[column]]
        # l-prefactor
        if ((self._l_prefactor and l_prefactor) or
            (not self._l_prefactor and not l_prefactor)):
            factor = 1
        elif self._l_prefactor and not l_prefactor:
            factor = 1/(self.ll()*(self.ll()+1)/(2.*np.pi))
        elif not self._l_prefactor and l_prefactor:
            factor = self.ll()*(self.ll()+1)/(2.*np.pi)
        try:
            return (factor * self._units_from_to(self._units, units, T_CMB) * 
                    self._lCl[self._columns_indices[column]])
        except KeyError:
            raise KeyError("It seems the '%s' lensed spectrum "%column+
                           "has not been calculated. "+
                           "The calcluated spectra are "+str(self._columns[1:]))
    def ll(self):
        return self.lCl("l")
    def llmax(self):
        return self.lCl("l")[-1]

    ### Internal
    # Change of units
    def _units_from_to(self, unit1, unit2, T_CMB = T_CMB):
        if unit1 == unit2:
            return 1
        factor_from_adim_to = {"1":   1,
                               "K":   (float(T_CMB))**2,
                               "muK": (float(T_CMB)*10**6)**2}
        return (factor_from_adim_to[unit1]**-1)*factor_from_adim_to[unit2]
