class NobleGasModl():
    def __init__(self):
        self.model_parameters = {
            'r_hop' : 3.1810226927827516, # hopping length scale
            't_ss' : 0.03365982238611262, # s-s hopping energy scale
            't_sp' : -0.029154833035109226, # s-p hopping energy scale
            't_pp1' : -0.0804163845390335, # 1st p-p hopping energy scale
            't_pp2' : -0.01393611496959445, # 2nd p-p hopping energy scale
            'r_pseudo' : 2.60342991362958, # pseudopotential length scale
            'v_pseudo' : 0.022972992186364977, # pseudopotential energy scale
            'dipole' : 2.781629275106456, # dipole strength of s-p transition
            'energy_s' : 3.1659446174413004, # onsite energy of s orbital
            'energy_p' : -2.3926873325346554, # onsite energy of p orbital
            'coulomb_s' : 0.3603533286088998, # Coulomb self-energy of monopole
            'coulomb_p' : -0.003267991835806299 # Coulomb self-energy of dipole
        }
        
        self.ionic_charge = 6
        self.orbital_types = ['s', 'px', 'py', 'pz']
        self.p_orbitals = self.orbital_types[1:]
        self.orbitals_per_atom = len(self.orbital_types)
        self.vec = {'px': [1, 0, 0], 'py': [0, 1, 0], 'pz': [0, 0, 1]}
        self.orbital_occupations = { 's':0, 'px':1, 'py':1, 'pz':1 }

    def __str__(self):
        return 'Noble Gas Model created with '+str(self.orbital_types)+' orbitals and an ionic charge of '+ str(self.ionic_charge) + '.'

    def atom(self, ao_index):
        '''Returns the atom index part of an atomic orbital index.


        Parameters
        ----------
        ao_index : int
            index of atomic orbital.

        Returns
        -------
        ao_index // orbitals_per_atom
            An integer which is the index of the atom that the atomic orbital is centered on.
        '''
        return ao_index // self.orbitals_per_atom

    def orb(self, ao_index):
        '''Returns the orbital type of an atomic orbital index.
        
        Parameters
        ----------
        ao_index : int
            index of atomic orbital.

        Returns
        -------
        ao_index % orbitals_per_atom
            A string which is the atomic orbital.
        '''
        orb_index = ao_index % self.orbitals_per_atom
        return self.orbital_types[orb_index]

    def ao_index(self, atom_p, orb_p):
        '''
        Returns the atomic orbital index for a given atom index and orbital type.
        
        Parameters
        ----------
        atom_p : int
            index of the atom.
        orb_p :  str
            type of the orbital
        Returns
        -------
        (atom_p * self.orbitals_per_atom ) + self.orbital_types.index(orb_p)
            The complete index of the atomic orbital.
        '''
        p = atom_p * self.orbitals_per_atom
        p += self.orbital_types.index(orb_p)
        return p


"""
Test to be put later of=n...
import numpy as np
atomic_coordinates = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 5.0]])
number_of_atoms = len(atomic_coordinates)
nbg1 = Noble_Gas_Model()
for index in range(number_of_atoms * nbg1.orbitals_per_atom):
    assert (index == nbg1.ao_index(nbg1.atom(index), nbg1.orb(index)))
"""
