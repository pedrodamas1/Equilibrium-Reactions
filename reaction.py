"""
Acetic acid dissociation example
HC2H3O2 <=> H+ + C2H3O2-
Ka = 1.8 * 10^-5
Initial condition: [HC2H3O2] = 0.50 M
"""

from atom import Atom
from molecule import Molecule
from typing import Dict, Union

class Reaction:
	def __init__(self, reactants: Dict[Union[Atom, Molecule], int], products: 
			  Dict[Union[Atom, Molecule], int], constant: float) -> None:
		self.reactants = reactants
		self.products = products
		self.constant = constant

	def get_coefficients(self):
		# must return [-moles1, -moles2, moles3, moles4]
		pass
	def get_charges(self):
		# must return the sign-sensitive charges
		pass

	def get_concentrations(self): 
		# must return the concentrations
		pass

		

if __name__ == '__main__':

	H = Atom('H')
	C = Atom('C')
	O = Atom('O')

	acetic_acid = Molecule(particles=((H,1), (C,2), (H,3), (O,2)))
	hydron = Molecule(particles=((H,1)), charge=+1)
	acetate = Molecule(particles=((C,2), (H,3), (O,2)), charge=-1)

	reaction1 = Reaction(
		reactants = {acetic_acid: 1},
		products={hydron:1, acetate:1},
		constant = 1.8e-5,
	)


# system = System(
#     reactions = Reactions(reaction1, reaction2, reaction3, etc),
#     initial_conditions = {'P': 0.4} # conserve phosphorus in phosphoric acid dissociation. Otherwise, how would it know?
# )
# system.solve()



# if __name__ == '__main__':
#     # HCl <=> H+ + Cl-
#     hydrochloric_acid = Molecule('HCl', 0)
#     hydron = Molecule('H', +1)
#     chloride = Molecule('Cl', -1)

#     hydrochloric_dissociation = Reaction(
#         reactants={
#             hydrochloric_acid:{'moles':1, 'concentration':0.1},
#         },
#         products={
#             hydron:{'moles':1, 'concentration':0.},
#             chloride:{'moles':1, 'concentration':0.},
#         },
#         constant=1e8,
#     )

#     system = System([hydrochloric_dissociation])
