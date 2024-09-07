import numpy as np
from scipy import interpolate
import os

def nload(file, lam):
    """
    Interpolates n&k data from file to a specific wavelength range

    Parameters:
    file (str): 1. Name of file with n&k data, two columns: wavelength (nm) | n
                2. Constant real number if starts with "const="
                3. n=1 if set to "air"
    lam (np.array): Wavelength (nm) for interpolation

    Returns:
    np.array: Interpolated complex refractive index: ref = n + 0j
    """
    print(f"Loading data for {file}")
    print(f"file: {file}, type: {type(file)}")
    print(f"lam: type: {type(lam)}, shape: {lam.shape}")

    # Ensure lam is a 1D array
    lam = lam.ravel()

    if file == "air":
        ref = np.ones_like(lam, dtype=complex)
        print(f"ref: type: {type(ref)}, shape: {ref.shape}")
        return ref
    elif file.startswith("const="):
        try:
            const_value = float(file.split("=")[1])
            print(f"const_value: {const_value}, type: {type(const_value)}")
            ref = np.full_like(lam, const_value, dtype=complex)
            print(f"ref: type: {type(ref)}, shape: {ref.shape}")
            return ref
        except ValueError:
            raise ValueError(f"Invalid constant value: {file}")
    else:
        if not os.path.exists(file):
            raise FileNotFoundError(f"File not found: {file}")
        
        data = np.loadtxt(file, delimiter=',')
        print(f"data: type: {type(data)}, shape: {data.shape}")

        wavelength = data[:, 0]
        n = data[:, 1]
        print(f"wavelength: type: {type(wavelength)}, shape: {wavelength.shape}")
        print(f"n: type: {type(n)}, shape: {n.shape}")

        # Create interpolation function
        f = interpolate.interp1d(wavelength, n, kind='linear', fill_value='extrapolate')
        print(f"f: type: {type(f)}")
        
        # Interpolate n values for given wavelengths
        n_interp = f(lam)
        print(f"Interpolated n values for {file}: first few: {n_interp[:5]}, last few: {n_interp[-5:]}")

        # Create complex refractive index (assuming k=0 as it's not provided in the CSV)
        ref = n_interp + 0j
        print(f"ref: type: {type(ref)}, shape: {ref.shape}")

        return ref
