from typing import Set, Dict
from reaction import Reaction
from molecule import Molecule
import numpy as np
from scipy.optimize import root, fsolve

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

		# Check if equations and unknowns are balanced
		N_unknowns = len(self.molecules)
		N_equations = len(self.reactions) + len(self.conservation) + 1
		is_balanced = N_unknowns == N_equations
		if not is_balanced:
			raise Exception(f'The number of unknows ({N_unknowns}) does not match the number of equations ({N_equations})')
		
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

		# Solve the system using the 'lm' algorithm (works better than 'hybr')
		sol = root(func, np.log10(initial_guess), method='lm')

		# Set the concetration value of each molecule
		if sol['success']:
			for molecule, c in zip(self.molecules, np.power(10, sol['x'])):
				molecule.concentration = c

		# Return the sucess of the solver
		return sol['success']

