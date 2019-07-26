class NobleGasModl():
    def __init__(self):
        self.model_parameters = {
         'coulomb_p': -0.010255409806855187,
         'coulomb_s': 0.4536486561938202,
         'dipole': 1.6692376991516769,
         'energy_p': -3.1186533988406335,
         'energy_s': 11.334912902362603,
         'r_hop': 2.739689713337267,
         'r_pseudo': 1.1800779720963734,
         't_pp1': -0.029546671673199854,
         't_pp2': -0.0041958662271044875,
         't_sp': 0.000450562836426027,
         't_ss': 0.0289251941290921,
         'v_pseudo': -0.015945813280635074}
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
