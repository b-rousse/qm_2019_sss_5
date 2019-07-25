class Noble_Gas_Model():
    def __init__(self, model_parameters, ionic_charge,orbital_types,p_orbitals,orbitals_per_atom,vec, orbital_occupations):
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
         'v_pseudo': -0.015945813280635074
                                }
        self.ionic_charge = 6
        self.orbital_types = ['s', 'px', 'py', 'pz']
        self.p_orbitals = orbital_types[1:]
        self.orbitals_per_atom = len(self.orbital_types)
        self.vec = {'px': [1, 0, 0], 'py': [0, 1, 0], 'pz': [0, 0, 1]}
        self.orbital_occupations = { 's':0, 'px':1, 'py':1, 'pz':1 }
        

    def __str__(self):
        return 'name: ' + self.surname + ', ' +self.name + '\ncourses: ' + str(self.courses)

    def atom(ao_index):
        '''Returns the atom index part of an atomic orbital index.

        This is a longer explanation of what the function does, NOT how the function does it. 

        Parameters
        ----------
        ao_index : int
            index of atomic orbital.

        Returns
        -------
        ao_index // orbitals_per_atom
            An integer which is the index of the atom that the atomic orbital is centered on.
        '''
        return ao_index // orbitals_per_atom

    def orb(ao_index):
        '''Returns the orbital type of an atomic orbital index.
        
        This is a more detailed description'''
        orb_index = ao_index % orbitals_per_atom
        return orbital_types[orb_index]

    def ao_index(atom_p, orb_p):
        '''Returns the atomic orbital index for a given atom index and orbital type.'''
        p = atom_p * orbitals_per_atom
        p += orbital_types.index(orb_p)
        return p
