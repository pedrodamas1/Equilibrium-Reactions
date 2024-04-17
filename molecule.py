from __future__ import annotations
from typing import Tuple, Set, Union
from collections import Counter
from atom import Atom


class Molecule:

	def __init__(self, particles: Tuple[Tuple[Union[Atom, Molecule], int]], 
			  charge: int = 0, concentration: float = 0.) -> None:
		self.particles = particles
		self.charge = charge
		self.concentration = concentration

	def __str__(self) -> str:
		# Create a placeholder for the string representation
		text = ''

		# Loop over the level 1 particles and moles in this molecule
		for particle, moles in self.particles:

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
		for particle, moles in self.particles:

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
		for particle, _ in self.particles:

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
	Fe = Atom('Fe')
	O = Atom('O')
	N = Atom('N')
	H = Atom('H')
	C = Atom('C')


	# Create a water molecule
	water = Molecule(particles = ((H,2), (O,1)))
	print(water)
	print(water.count_species())
	print(water.get_species())

	# Create a level one molecule
	nitrate = Molecule(particles = ((N,1), (O,3)), charge = -1)
	print(nitrate)
	print(nitrate.print_ion())
	print(nitrate.count_species())
	print(nitrate.get_species())

	# Create a level two molecule (or composite molecule)
	iron_nitrate = Molecule(particles = ((Fe,2), (nitrate,3)), charge = +3)
	print(iron_nitrate)
	print(iron_nitrate.print_ion())
	print(iron_nitrate.count_species())
	print(iron_nitrate.get_species())

	# Acetic acid: HC2H3O2 - notice how we have separate H representations
	acetic_acid = Molecule(particles=((H,1), (C,2), (H,3), (O,2)))
	print(acetic_acid)
	print(acetic_acid.get_species())
	print(acetic_acid.count_species())
