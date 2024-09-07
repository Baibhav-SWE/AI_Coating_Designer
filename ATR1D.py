import numpy as np

def ATR1D(stack):
    """
    Calculates absorptance, transmittance and reflectance for 1D multilayered stack

    Parameters:
    stack (dict): Stack dictionary, see 'set_stack.py' for details

    Returns:
    tuple: (A, T, R) where each is a dict containing 's', 'p', and 'sp' components
    """
    # Checking input and allocating workspace
    lam = stack["wavelength"]
    d = stack["thick"]
    theta0 = stack["angle"]
    n = stack["nk"]

    rs = np.zeros((len(lam), len(d)+1), dtype=complex)
    rp = np.zeros((len(lam), len(d)+1), dtype=complex)
    ts = np.zeros((len(lam), len(d)+1), dtype=complex)
    tp = np.zeros((len(lam), len(d)+1), dtype=complex)
    t = {"s": [[None for _ in range(len(d)+1)] for _ in range(len(lam))],
         "p": [[None for _ in range(len(d)+1)] for _ in range(len(lam))]}
    ph = [[None for _ in range(len(d))] for _ in range(len(lam))]

    R = {"s": np.zeros(len(lam)), "p": np.zeros(len(lam))}
    T = {"s": np.zeros(len(lam)), "p": np.zeros(len(lam))}
    I = {"s": [None] * len(lam), "p": [None] * len(lam)}

    # Calculating phase shifts for each layer
    theta = np.zeros(n.shape, dtype=complex)
    theta[:, 0] = theta0

    for i in range(1, n.shape[1]):
        theta[:, i] = np.arcsin(n[:, i-1] / n[:, i] * np.sin(theta[:, i-1]))  # Snell's law

    delt = 2 * np.pi * n[:, 1:-1] * np.outer(1/lam, d) * np.cos(theta[:, 1:-1])
    edp = np.exp(1j * delt)
    edm = 1 / edp

    # Calculating Fresnel coefficients for each boundary
    for i in range(len(d) + 1):
        rs[:, i] = (n[:, i] * np.cos(theta[:, i]) - n[:, i+1] * np.cos(theta[:, i+1])) / \
                   (n[:, i] * np.cos(theta[:, i]) + n[:, i+1] * np.cos(theta[:, i+1]))
        rp[:, i] = (n[:, i] * np.cos(theta[:, i+1]) - n[:, i+1] * np.cos(theta[:, i])) / \
                   (n[:, i] * np.cos(theta[:, i+1]) + n[:, i+1] * np.cos(theta[:, i]))
        ts[:, i] = 2 * n[:, i] * np.cos(theta[:, i]) / \
                   (n[:, i] * np.cos(theta[:, i]) + n[:, i+1] * np.cos(theta[:, i+1]))
        tp[:, i] = 2 * n[:, i] * np.cos(theta[:, i]) / \
                   (n[:, i] * np.cos(theta[:, i+1]) + n[:, i+1] * np.cos(theta[:, i]))

    # Calculating transfer matrix elements for each layer
    for i in range(len(lam)):
        t["s"][i][0] = np.array([[1, rs[i, 0]], [rs[i, 0], 1]]) / ts[i, 0]
        t["p"][i][0] = np.array([[1, rp[i, 0]], [rp[i, 0], 1]]) / tp[i, 0]
        for j in range(1, len(d)+1):
            ph[i][j-1] = np.array([[edm[i, j-1], 0], [0, edp[i, j-1]]])
            t["s"][i][j] = np.array([[1, rs[i, j]], [rs[i, j], 1]]) / ts[i, j]
            t["p"][i][j] = np.array([[1, rp[i, j]], [rp[i, j], 1]]) / tp[i, j]

    # Getting S & P transmittance & reflectance
    cfTs = np.real(n[:, -1] * np.cos(theta[:, -1]) / n[:, 0] / np.cos(theta0))
    cfTp = np.real(np.conj(n[:, -1]) * np.cos(theta[:, -1]) / np.conj(n[:, 0]) / np.cos(theta0))

    nincoh = stack["nincoh"]
    for i in range(len(lam)):
        I["s"][i] = np.eye(2)
        I["p"][i] = np.eye(2)
        j = 0
        while j < len(d):
            tmp = {"s": t["s"][i][j], "p": t["p"][i][j]}
            while nincoh[j]:
                tmp["s"] = tmp["s"] @ ph[i][j] @ t["s"][i][j+1]
                tmp["p"] = tmp["p"] @ ph[i][j] @ t["p"][i][j+1]
                j += 1
                if j >= len(d):
                    I["s"][i] = I["s"][i] @ np.abs(tmp["s"])**2
                    I["p"][i] = I["p"][i] @ np.abs(tmp["p"])**2
                    break
            if j < len(d) and not nincoh[j]:
                I["s"][i] = I["s"][i] @ np.abs(tmp["s"])**2 @ np.abs(ph[i][j])**2
                I["p"][i] = I["p"][i] @ np.abs(tmp["p"])**2 @ np.abs(ph[i][j])**2
                j += 1
                if j >= len(d):
                    I["s"][i] = I["s"][i] @ np.abs(t["s"][i][j])**2
                    I["p"][i] = I["p"][i] @ np.abs(t["p"][i][j])**2
                    break
                elif not nincoh[j]:
                    I["s"][i] = I["s"][i] @ np.abs(t["s"][i][j])**2 @ np.abs(ph[i][j])**2
                    I["p"][i] = I["p"][i] @ np.abs(t["p"][i][j])**2 @ np.abs(ph[i][j])**2
                    j += 1
                    if j >= len(d):
                        I["s"][i] = I["s"][i] @ np.abs(t["s"][i][j])**2
                        I["p"][i] = I["p"][i] @ np.abs(t["p"][i][j])**2
                        break
        R["s"][i] = I["s"][i][1, 0] / I["s"][i][0, 0]
        R["p"][i] = I["p"][i][1, 0] / I["p"][i][0, 0]
        T["s"][i] = cfTs[i] / I["s"][i][0, 0]
        T["p"][i] = cfTp[i] / I["p"][i][0, 0]

    # Constructing output
    T["sp"] = (T["s"] + T["p"]) / 2
    R["sp"] = (R["s"] + R["p"]) / 4 * 2
    A = {
        "s": 1 - T["s"] - R["s"],
        "p": 1 - T["p"] - R["p"],
        "sp": 1 - T["sp"] - R["sp"]
    }

    return A, T, R