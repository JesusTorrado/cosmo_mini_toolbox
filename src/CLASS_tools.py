######################################
# Small set of tools for interacting #
# with the CLASS Boltzmann code      #
######################################

import os
from tempfile import mkdtemp
from shutil import rmtree

from CMBspectrum import CMBspectrum

def read_CLASS_param_file(param_file):
    params = {}
    with open(param_file, "r") as pfile:
        for line in param_file :
            # Take everything at the left of a comment symbol
            aux = line.split("#")[0]
            if not(aux):
                continue
            if not "=" in aux :
                continue
            aux = [a.strip() for a in aux.split("=")]
            if aux[1]:
                params[aux[0]] = aux[1]
    return params

def write_CLASS_param_file(params, name, folder, verbose=1,
                             output="tCl,pCl,lCl"):
    """
    Generates a CLASS parameter file from a given dictionary of parameter
    values "params" with a given name "[name].ini" in the folder "folder".

    Returns the full name of the written parameters file, including path.

    Optional arguments:
    -------------------

    verbose: int (default: 1)
        If True, sets the *_verbose parameters of CLASS to the value of the
        given value (0 for silent, 1 for verbose, >1 for extra verbose)

    output: str (default: "tCl,pCl,lCl")
        Value of the "output" parameter: which spectra is to be requested.

    """
    # products
    params["output"] = output
    if "lCl" in output:
        params["lensing"] = "yes"
    # verbose
    verbose_params = ["background", "thermodynamics", "perturbations",
                      "bessels", "transfer", "primordial", "spectra",
                      "nonlinear", "lensing", "output"]
    for p in verbose_params:
        params[p+"_verbose"] = verbose
    # file
    assert os.path.exists(folder) ,\
        "The given folder does not exist: " + folder
    full_name = os.path.join(folder, name+".ini")
    with open(full_name, "w") as param_file:
        for param, value in params.items():
            param_file.write(param + " = " + str(value) + "\n")
    return full_name

# Run Class and output lines
def run_CLASS(class_folder, param_file, precision_file=None, verbose=False):
    """
    Runs CLASS (the 'class' binary at the given 'class_folder') as a subprocess,
    using the parameter file 'param_file' and, if specified, the precision file
    'precision_file'.

    If 'verbose' is set to True, the output of CLASS is printed.

    Returns the output of CLASS as a list of lines.
    """
    assert os.path.isfile(param_file), "The given parameter file does not exist!"
    class_command = " ".join([os.path.join(class_folder, "class"), param_file,
                              precision_file if precision_file else ""])
    if verbose:
        print "Running Class as: '%s'"%class_command
    lines = os.popen(class_command).readlines()
    if verbose:
        print "".join(lines)
    return lines

# Create a CMBspectrum instance from CLASS
def CMBspectrum_from_param_file_CLASS(class_folder, param_file,
                                      precision_file=None, verbose=False):
    """
    Generates a CMBspectrum instance from a CLASS param file 'param_file',
    running the CLASS instance in 'class_folder'.

    If 'param_file' is a dictionary, a temporary param file is generated
    automatically.

    For help on the arguments, see the documentation of 'run_CLASS'.
    """
    folder = mkdtemp()
    # If param_file is a file -> dictionary:
    if not isinstance(param_file, dict):
        param_file = read_CLASS_param_file(param_file)
    # Else, use it directly as a dictionary
    param_file["root"] = folder+"/" # BAD behaviour of CLASS
    param_file["write parameters"] = "yes"
    param_file_name = "tmp_params"
    param_file = write_CLASS_param_file(param_file, param_file_name,
                                        folder, verbose=1,
                                        output="tCl,pCl,lCl")
    lines = run_CLASS(class_folder, param_file, precision_file=precision_file,
                      verbose=verbose)
    spectrum = CMBspectrum(folder, prefix=None, code="CLASS")
    # Erase the tmp folder
    rmtree(folder)
    return spectrum
