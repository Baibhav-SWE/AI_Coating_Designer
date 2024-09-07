import numpy as np
from nload import nload

def set_stack(mat, thick, lam, angle, *args):
    """
    Creates dictionary for the multilayered stack

    Parameters:
    mat (list): List of strings with names for files containing n&k for materials
    thick (list): Thicknesses of layers, notice first and last layer are infinite
    lam (np.array): Wavelengths
    angle (float): Angle of incidence
    *args: Optional arguments:
           1. incoh - threshold thickness for incoherence in layers

    Returns:
    dict: Dictionary of stack ready for further analysis
    """
    # Modifying input if needed
    lam = np.atleast_2d(lam).T

    # Interpolating N&K values
    mat_unique, ic = np.unique(mat, return_inverse=True)
    n_unique = np.zeros((len(lam), len(mat_unique)), dtype=complex)

    for imat, material in enumerate(mat_unique):
        n_unique[:, imat] = nload(material, lam)
        print(f"Material: {material}, first few n values: {n_unique[:5, imat]}")

    n = n_unique[:, ic]
    print(f"n shape: {n.shape}, first few values of first layer: {n[:5, 0]}")

    # Creating stack dictionary
    if len(args) == 1:
        nincoh = np.array(thick) < args[0]
    else:
        nincoh = np.ones_like(thick, dtype=bool)

    stack = {
        "nk": n,
        "thick": thick,
        "nincoh": nincoh,
        "wavelength": lam,
        "angle": angle
    }

    return stack
