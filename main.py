import numpy as np
import matplotlib.pyplot as plt
from set_stack import set_stack
from ATR1D import ATR1D

# General input data
lam = np.linspace(380, 1080, 851)
print("new statement")
Bragg1 = ["MgF2.csv","MgF2.csv","MgF2.csv","MgF2.csv"]  # First Bragg stack
Bragg2 = ["MgF2.csv"]  # Second Bragg stack
matPSC = Bragg1 + ["GLS_NEW.csv"] + Bragg2  # Combining Bragg1 and Bragg2 with glass in between
dgls = 1000
dBragg1 = [110.18,100,300,400]  # Thicknesses for Bragg1 (SiO2/TiO2 x2)
dBragg2 = [100]  # Thicknesses for Bragg2 (TiO2/SiO2 x2)
dPSC = dBragg1 + [dgls] + dBragg2
incoh = 1e3  # incoherent layer is 1000nm thickness
theta = 0
dAZO=0
dEVA=0
dARC=0      
# PSC with AZO and EVA on PV
materials = ["air"] + matPSC + ["air", "air", "air", "air"]
d = dPSC + [dAZO, dEVA, dARC]  # Note: dAZO, dEVA, dARC are not defined in the original script
stack = set_stack(materials, d, lam, theta, incoh)
_, T_PSC_AZO_EVA_PV, _ = ATR1D(stack)

# Add these print statements
print(f"T['sp'] shape: {T_PSC_AZO_EVA_PV['sp'].shape}")
print(f"First 5 T values: {T_PSC_AZO_EVA_PV['sp'][:5]}")
print(f"Last 5 T values: {T_PSC_AZO_EVA_PV['sp'][-5:]}")


# Plotting results
plt.figure()
plt.plot(lam, T_PSC_AZO_EVA_PV['sp'], 'b', linewidth=1.5, label='Transmittance (T)')

print(T_PSC_AZO_EVA_PV['sp'])
plt.xlabel('Wavelength [nm]')
plt.ylabel('Fraction')
plt.xlim(lam[0], lam[-1])
plt.ylim(0.93, 0.98)
plt.legend()
plt.title('TRA Plot for PSC Structure with Separate Bragg Stacks')

plt.show()
