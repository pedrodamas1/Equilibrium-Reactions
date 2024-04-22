
class Molecule:
	"""
	Represents a molecule with its name, charge, and concentration.
	"""

	def __init__(self, name: str, charge: int = 0, concentration: float = 0.0) -> None:
		"""
		Initialize the Molecule object.

		Args:
			name (str): The name of the molecule.
			charge (int): The charge of the molecule (default is 0).
			concentration (float): The concentration of the molecule (default is 0.0).
		"""
		self.name = name
		self.charge = charge
		self.concentration = concentration

	def __str__(self) -> str:
		"""
		Return a string representation of the molecule.
		"""
		return self.name

	def __repr__(self) -> str:
		"""
		Return a string representation that can be used to recreate the object.
		"""
		return self.name
