from __future__ import annotations
from typing import Dict, Set
import numpy as np
from collections import Counter
from atom import Atom


class Molecule:

	def __init__(self, particles: Dict[Molecule, int], charge: int) -> None:
		self.particles = particles
		self.charge = charge

	def __str__(self) -> str:
		# Create a placeholder for the string representation
		text = ''

		# Loop over the level 1 particles and moles in this molecule
		for particle, moles in self.particles.items():

			# Defined the template depending on the particle type and the number of moles 
			if isinstance(particle, Molecule) and moles != 1:
				template = '({}){}'
			elif isinstance(particle, Atom) and moles != 1:
				template = '{}{}'
			else:
				template = '{}'

			# Add the representation of this particle to the output string
			text += template.format(particle, moles)
		
		# Return the string
		return text
	
	def print_ion(self):
		return f'{self}{self.charge:+d}'
	
	def count_species(self) -> Counter:
		# Create an empty counter instance
		counter = Counter()

		# Loop over the particles and moles in this molecule
		for particle, moles in self.particles.items():

			# If this particle is an atom, we just add its moles to the counter
			if isinstance(particle, Atom):
				counter[particle] += moles
			
			# If instead it is a molecule, we fetch its species recursively
			else:
				for subparticle, submoles in particle.count_species().items():
					counter[subparticle] += submoles * moles
		
		# Return the counter object
		return counter

	def get_species(self) -> Set[Atom]:
		# Create an empty set instance
		species = set()

		# Loop over the particles in this molecule
		for particle in self.particles:

			# If this particle is an atom, we just add it to the set
			if isinstance(particle, Atom):
				species.add(particle)
			
			# If instead it is a molecule, we fetch its species recursively
			else:
				for subparticle in particle.get_species():
					species.add(subparticle)
		
		# Return the set object
		return species


if __name__ == '__main__':
	# Create a few elements
	iron = Atom('Fe')
	oxygen = Atom('O')
	nitrogen = Atom('N')
	hydrogen = Atom('H')

	# Create a water molecule
	water = Molecule(particles = {hydrogen:2, oxygen:1}, charge = 0)
	print(water)
	print(water.count_species())
	print(water.get_species())

	# Create a level one molecule
	nitrate = Molecule(particles = {nitrogen:1, oxygen:3}, charge = -1)
	print(nitrate)
	print(nitrate.print_ion())
	print(nitrate.count_species())
	print(nitrate.get_species())

	# Create a level two molecule (or composite molecule)
	iron_nitrate = Molecule(particles = {iron:2, nitrate:3}, charge = +3)
	print(iron_nitrate)
	print(iron_nitrate.print_ion())
	print(iron_nitrate.count_species())
	print(iron_nitrate.get_species())





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
