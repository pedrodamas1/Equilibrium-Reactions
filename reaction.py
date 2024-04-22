from typing import Dict
from molecule import Molecule

class Reaction:
	"""
	Represents a chemical reaction with participating species and equilibrium constant.
	"""

	def __init__(self, species: Dict[Molecule, int], equilibrium_constant: float) -> None:
		"""
		Initialize the Reaction object.

		Args:
			species (Dict[Molecule, int]): Dictionary containing species and their stoichiometric coefficients.
			equilibrium_constant (float): The equilibrium constant of the reaction.
		"""
		self.species = species
		self.equilibrium_constant = equilibrium_constant

	def __str__(self) -> str:
		"""
		Return a string representation of the chemical reaction.
		"""
		reactants = [f"{'' if abs(moles) == 1 else abs(moles)}{molecule}" for molecule, moles in self.species.items() if moles < 0]
		products = [f"{'' if abs(moles) == 1 else abs(moles)}{molecule}" for molecule, moles in self.species.items() if moles > 0]
		return ' + '.join(reactants) + ' <=> ' + ' + '.join(products)

	def __repr__(self) -> str:
		"""
		Return a string representation.
		"""
		return self.__str__()
