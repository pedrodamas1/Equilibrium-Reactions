
'''
SYSTEMATIC ACID-BASE CALCULATIONS

To run use: python -m examples.HCl_NaOH_titration

Reactions in place:
    HCl <=> H+ + Cl-
    NaOH <=> Na+ + OH-
    H2O <=> H+ + OH-
'''

from molecule import Molecule
from reaction import Reaction
from system import System
import numpy as np

# Create some substances
hydrochloric_acid = Molecule(name = 'HCl', charge = 0.)
hydron = Molecule(name = 'H+', charge = +1)
chloride = Molecule(name = 'Cl-', charge = -1)
sodium_hydroxide = Molecule(name = 'NaOH')
sodium_ion = Molecule(name = 'Na+', charge = +1)
hydroxide = Molecule(name = 'OH-', charge = -1)

# Create the reactions present
hydrochloric_dissociation = Reaction(species = {hydrochloric_acid: -1, hydron: 1, chloride: 1}, equilibrium_constant = 1.**6.1)
sodium_hydroxide_dissociation = Reaction(species = {sodium_hydroxide: -1, sodium_ion: 1, hydron: 1}, equilibrium_constant = 1.**-0.2)
water_dissociation = Reaction(species = {hydroxide:1, hydron:1}, equilibrium_constant = 1.0e-14)

# # Create and solve the reaction system
# system = System(reactions = {hydrochloric_dissociation, sodium_hydroxide_dissociation, water_dissociation}, conservation = {'Cl':0.02, 'Na':0.015})
# system.solve(initial_guess=np.ones(len(system.molecules))*1e-2)

# # Use this website to validate results: https://www.aqion.onl/
# print(f'The equilibrium concentration of hydron is [H+]={hydron.concentration:.5f} and the pH is {-np.log10(hydron.concentration):.3f}')

# Create and solve the reaction system
system = System(reactions = {hydrochloric_dissociation, sodium_hydroxide_dissociation, water_dissociation}, conservation = {'Cl':1.0, 'Na':0.9})

success = system.brute_solve(0.9, 1,1)
cH = hydron.concentration
pH = -np.log10(hydron.concentration)
pH_arr = [pH]
print(f'Success: {success} [H+]={cH:.10f} pH: is {pH:.5f}')

# Now use the initial result to move along small steps
guess = system.concentrations
ci_arr = np.linspace(0.9, 1.1, 100)
for ci in ci_arr[1:]:
	system.conservation['Na'] = ci
	sucess = system.solve(initial_guess=guess)
	guess = system.concentrations
	cH = hydron.concentration
	pH = -np.log10(cH)
	pH_arr.append(pH)
	print(f'Success: {success} [H+]={cH:.10f} pH: is {pH:.5f}')

import matplotlib.pyplot as plt

plt.plot(ci_arr, pH_arr)
plt.xlabel("Concentration of NaOH [M]")
plt.ylabel("pH")
plt.grid()
plt.show()


exit()



from scipy.optimize import fsolve

def my_func(zz, *args):

	C_HCl, C_NaOH = args
	H, Cl, HCl, Na, OH, NaOH = zz
	pka, pkb, pkw = -6.1, 0.2, 14

	f = np.zeros(6)

	f[0] = pka + H + Cl - HCl
	f[1] = pkb + Na + OH - NaOH 
	f[2] = pkw + H + OH
	f[3] = 10**Cl + 10**HCl - C_HCl
	f[4] = 10**Na + 10**NaOH - C_NaOH
	f[5] = 10**Cl + 10**OH - 10**Na - 10**H

	return f


Ca= 1.0
pH_values, Cb_values = [], []

for Cb in np.linspace(0.9, 1.1, 100):

	init_guess = np.log10( np.array([Ca, Ca, Ca, Cb, Cb, Cb]) * 1.E-2)
	z = fsolve(my_func, init_guess, args=(Ca, Cb), full_output=True)
	
	if z[3] !=  "The solution converged.":

		print("The solution did not converge")

		init_guess = np.log10( np.array([Ca, Ca, Ca, Ca, Cb, Cb]) * 1.E-2)
		z = fsolve(my_func, init_guess, args=(Ca, Cb), full_output=True)

		if z[3] !=  "The solution converged.":
			print("Still not converging.")
		else:
			print("But now it did!")

	concentrations = 10**z[0]
	H, Cl, HCl, Na, OH, NaOH = concentrations
	pH = -np.log10(H)
	print(pH, Cb)

	if z[3] !=  "The solution converged.":
		pH = pH_values[-1]

	pH_values.append(pH)
	Cb_values.append(Cb)


import matplotlib.pyplot as plt

plt.plot(Cb_values, pH_values)
plt.xlabel("Concentration of NaOH [M]")
plt.ylabel("pH")
plt.grid()
plt.show()

