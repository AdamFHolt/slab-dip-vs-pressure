#!/bin/python
import numpy as np

### conditions for ref. visc ###
visc_ref   = 2.5e20
slab_visc_contrast = 250
lm_visc_contrast   = 40

### other params
temp_ref   = 1573.
Ediff = 0 
ndiff = 1.       
R = 8.314

# compute background upper mantle diffusion creep prefactor
Adiff = 0.5 * (1./(visc_ref)) * np.exp((Ediff)/(R*temp_ref))
visc_check_diff = 0.5 * (1./(Adiff)) * np.exp((Ediff)/(R*temp_ref))
print("Bacground diffusion prefactor     = %e. (A check: %e = %e)" % (Adiff,visc_ref,visc_check_diff))

Adiff_lith    = 0.5 * (1./(visc_ref*slab_visc_contrast)) * np.exp((Ediff)/(R*temp_ref))
Adiff_lm      = 0.5 * (1./(visc_ref*lm_visc_contrast))   * np.exp((Ediff)/(R*temp_ref))
Adiff_lith_lm = 0.5 * (1./(visc_ref*slab_visc_contrast*lm_visc_contrast)) * np.exp((Ediff)/(R*temp_ref))

visc_check_lith    = 0.5 * (1./(Adiff_lith)) * np.exp((Ediff)/(R*temp_ref))
visc_check_lm      = 0.5 * (1./(Adiff_lm))   * np.exp((Ediff)/(R*temp_ref))
visc_check_lith_lm = 0.5 * (1./(Adiff_lith_lm))   * np.exp((Ediff)/(R*temp_ref))

print("Lithosphere diffusion prefactor   = %e. (A check: %e = %e)" % (Adiff_lith,visc_ref*slab_visc_contrast,visc_check_lith))
print("Lower mantle diffusion prefactor  = %e. (A check: %e = %e)" % (Adiff_lm,  visc_ref*lm_visc_contrast,  visc_check_lm))
print("Lith in the Lower mantle diffusion prefactor  = %e. (A check: %e = %e)" % (Adiff_lith_lm,  visc_ref*slab_visc_contrast*lm_visc_contrast,  visc_check_lith_lm))
