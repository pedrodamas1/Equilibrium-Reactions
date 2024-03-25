from __future__ import annotations
from typing import Dict
import numpy as np
from collections import Counter

class Atom:
    def __init__(self, name: str) -> None:
        self.name = name
    def __str__(self) -> str:
        return self.name
    def __repr__(self) -> str:
        return self.__str__()

iron = Atom('Fe')
nitrogen = Atom('N')
oxygen = Atom('O')

class Molecule:
    def __init__(self, particles: Dict[Molecule, int], charge: int) -> None:
        self.particles = particles
        self.charge = charge
    def __str__(self) -> str:
        text = ''
        for particle, moles in self.particles.items():
            if isinstance(particle, Molecule) and moles != 1:
                template = '({}){}'
            elif isinstance(particle, Atom) and moles != 1:
                template = '{}{}'
            else:
                template = '{}'
            text += template.format(particle, moles)
        return text
    def get_ion(self):
        return f'{self}{self.charge:+d}'
    def count_species(self):
        counter = Counter()
        for particle, moles in self.particles.items():
            if isinstance(particle, Atom):
                counter[particle] += moles
            else:
                for key,val in particle.count_species().items():
                    counter[key] += val*moles
                # counter += Counter({key:val*moles for key,val in particle.count_species().items()})
        print([[particle]*moles for particle,moles in self.particles.items()])
        return counter




# K4(Fe(CN)6) — Potassium ferrocyanide

# (Co(NH3)5Cl)2+ — Cobalt pentaamine chloride

nitrate = Molecule({nitrogen:1, oxygen:3}, 0)
print(nitrate.count_species())
iron_nitrate = Molecule({iron:2, nitrate:3}, +3)
print(iron_nitrate.count_species())




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
