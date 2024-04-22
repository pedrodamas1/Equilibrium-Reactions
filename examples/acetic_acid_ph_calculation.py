
'''
To run use: python -m examples.acetic_acid_ph_calculation

Reactions in place:
    HC2H3O2 <=> C2H3O2- + H+
    H2O <=> OH- + H+
'''

from molecule import Molecule
from reaction import Reaction
from system import System
import numpy as np

# Create some substances
acetic_acid = Molecule(name = 'HC2H3O2')
hydron = Molecule(name = 'H+', charge = +1)
acetate = Molecule(name = 'C2H3O2-', charge = -1)
hydroxide = Molecule(name = 'OH-', charge = -1)

# Create the reactions present
acetic_dissociation = Reaction(species = {acetic_acid: -1, hydron: 1, acetate: 1}, equilibrium_constant = 1.8e-5)
water_dissociation = Reaction(species = {hydroxide:1, hydron:1}, equilibrium_constant = 1.0e-14)

# Create and solve the reaction system
system = System(reactions = {acetic_dissociation, water_dissociation}, conservation = {'C2H3O2':0.5})
system.solve()

# Use this website to validate results: https://www.aqion.onl/
print(f'The equilibrium concentration of hydron is [H+]={hydron.concentration:.5f} and the pH is {-np.log10(hydron.concentration):.3f}')

