from typing import Dict, Set
import numpy as np
from scipy.optimize import root


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


class System:
	"""
	Represents a chemical reaction system.
	"""

	def __init__(self, reactions: Set[Reaction], conservation: Dict[str, float]) -> None:
		"""
		Initialize the System object.

		Args:
			reactions (Set[Reaction]): Set of reactions in the system.
			conservation (Dict[str, float]): Dictionary containing conservation conditions.
		"""
		self.reactions = reactions
		self.conservation = conservation
	
	@property
	def molecules(self) -> Set[Molecule]:
		"""
		Get the set of all molecules involved in the system.
		"""
		return {molecule for reaction in self.reactions for molecule in reaction.species}

	@property
	def molar_matrix(self) -> np.ndarray:
		"""
		Calculate the molar matrix of the system.
		"""
		
		# Get the number of molecules
		N = len(self.molecules)

		# Create a dictionary of molecules to ordinal values
		species_dict = {molecule: i for i, molecule in enumerate(self.molecules)}

		# Create an empty matrix of zeros
		mat = np.zeros((len(self.reactions), N))

		# Loop over the reactions to fill into the matrix
		for i, reaction in enumerate(self.reactions):

			# Loop over the molecule and moles in each reaction
			for molecule, moles in reaction.species.items():

				# Get the column number corresponding to the current molecule
				j = species_dict[molecule]

				# Set the matrix cell to the number of moles (sign-sensitive)
				mat[i,j] = moles
		
		# Return the matrix
		return mat
	
	@property
	def conservation_matrix(self) -> np.ndarray:
		"""
		Calculate the conservation matrix of the system.
		"""

		# Get the number of molecules
		N = len(self.molecules)

		# Create a dictionary of molecules to ordinal values
		species_dict = {molecule: i for i, molecule in enumerate(self.molecules)}

		# Create an empty matrix of zeros
		mat = np.zeros((len(self.conservation), N))

		# Loop over each conservation condition
		for i, specie in enumerate(self.conservation):

			# Loop over each molecule in the system
			for molecule in self.molecules:

				# Check if this molecule contains our conserved specie
				if specie in molecule.name:

					# Get the column number corresponding to the current molecule
					j = species_dict[molecule]

					# Set the matrix cell to 1
					mat[i,j] = 1
		
		# Return the matrix
		return mat

	def solve(self, initial_guess: np.ndarray = None) -> bool:
		"""
		Solve the system of chemical reactions.

		Args:
			initial_guess (np.ndarray, optional): Initial guess for concentrations. Defaults to None.

		Returns:
			bool: True if the solver succeeds, False otherwise.
		"""
		
		def func(log_c_arr: np.ndarray) -> np.ndarray:
			"""
			Calculate the function for the solver.

			Args:
				log_c_arr (np.ndarray): Array of the logarithm of species concentrations.

			Returns:
				np.ndarray: Array of equations to solve.
			"""

			# Convert log_c_arr to C_arr = [C1, C2, C3, etc]
			c_arr = np.power(10, log_c_arr)

			# Get the reaction constants and the pK
			k_arr = np.array([reaction.equilibrium_constant for reaction in self.reactions])
			pk_arr = -np.log10(k_arr)

			# Get the molar matrix with sign-sensitive stoichiometric coefficients
			M_mat = self.molar_matrix

			# Equilibrium equations
			f1 = pk_arr + M_mat @ log_c_arr

			# Get the initial condition concentrations
			ci_arr = np.fromiter(self.conservation.values(), dtype=float)

			# Get conservation matrix which identifies conserved species
			C_mat = self.conservation_matrix

			# Mass conservation equation
			f2 = ci_arr - C_mat @ c_arr

			# Get the charge array of species
			q_arr = np.array([specie.charge for specie in self.molecules])

			# Charge conservation equations
			f3 = np.array([np.dot(q_arr, c_arr)])

			# Concatenate all and return
			return np.concatenate((f1,f2,f3))
		
		# If an initial guess was not provided, set one
		if initial_guess is None:
			initial_guess = np.ones(len(self.molecules))

		# Solve the system
		sol = root(func, np.log10(initial_guess))

		# Set the concetration value of each molecule
		if sol['success']:
			for molecule, c in zip(self.molecules, np.power(10, sol['x'])):
				molecule.concentration = c

		# Return the sucess of the solver
		return sol['success']


if __name__ == '__main__':
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
	print(f'The equilibrium concentration of hydron is [H+]={hydron.concentration:.5f} and \
	   the solution pH is {-np.log10(hydron.concentration):.3f}')



