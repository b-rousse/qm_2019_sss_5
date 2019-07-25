import numpy as np

class HartreeFock:
    def __init__(self, NobleGasModel, atomic_coordinates):
        self.atomic_coordinates = atomic_coordinates
        #self.gas_model = NobleGasModel...
        self.ndof = len(self.atomic_coordinates) * NobleGasModel.orbitals_per_atom
        self.interaction_matrix = self.calculate_interaction_matrix(self.atomic_coordinates, self.ndof, NobleGasModel.model_parameters)
        self.hamiltonian_matrix = self.calculate_hamiltonian_matrix(atomic_coordinates, self.ndof, NobleGasModel.model_parameters)
        self.density_matrix = None
        self.chi_tensor = None
        self.fock_matrix = self.calculate_fock_matrix(self.hamiltonian_matrix, self.interaction_matrix, self.density_matrix, self.chi_tensor)#check for error
#        self.density_matrix = self.calculate_density_matrix(self.fock_matrix)
#        self.chi_tensor = self.calculate_chi_tensor(self.atomic_coordinates, self.ndof ,NobleGasModel.model_parameters) 



    def calculate_interaction_matrix(self, atomic_coordinates, ndof, model_parameters):

        '''Returns the electron-electron interaction energy matrix for an input list of atomic coordinates.'''
        interaction_matrix = np.zeros( (ndof,ndof) )
        for p in range(ndof):
            for q in range(ndof):
                if atom(p) != atom(q):
                    r_pq = atomic_coordinates[atom(p)] - atomic_coordinates[atom(q)]
                    interaction_matrix[p,q] = coulomb_energy(orb(p), orb(q), r_pq)
                if p == q and orb(p) == 's':
                    interaction_matrix[p,q] = model_parameters['coulomb_s']
                if p == q and orb(p) in p_orbitals:
                    interaction_matrix[p,q] = model_parameters['coulomb_p']



        return interaction_matrix


    def coulomb_energy(self, o1, o2, r12):
        """Returns Coulomb Matrix element for a pair of multipoles of type o1 and o2 separated by vector r12

        Parameters
        ----------
        o1 : orbital 1, string 
        o2 : orbital 2, string
        r12 : hopping lenght scale, numpy array 

        Returns
        -------
        Coulomb Matrix Element

        """

        r12_length = np.linalg.norm(r12)
        if o1 == 's' and o2 == 's':
            coulomb_energy  = 1.0 / r12_length
        if o1 == 's' and o2 in p_orbitals:
            coulomb_energy = np.dot(vec[o2], r12) / r12_length**3
        if o2 == 's' and o1 in p_orbitals:
            coulomb_energy = -np.dot(vec[o1], r12) / r12_length**3
        if o1 in p_orbitals and o2 in p_orbitals:
            coulomb_energy = ( np.dot(vec[o1], vec[o2]) / r12_length**3
                   - 3.0 * np.dot(vec[o1], r12) * np.dot(vec[o2], r12) / r12_length**5 )
        return coulomb_energy



    def calculate_potential_vector(self, atomic_coordinates, ndof, model_parameters):

        '''Returns the electron-ion potential energy vector for an input list of atomic coordinates.'''

        potential_vector = np.zeros(ndof)
        for p in range(ndof):
            potential_vector[p] = 0.0
            for atom_i,r_i in enumerate(atomic_coordinates):
                r_pi = atomic_coordinates[atom(p)] - r_i
                if atom_i != atom(p):
                    potential_vector[p] += ( pseudopotential_energy(orb(p), r_pi, model_parameters) 
                                             - ionic_charge * coulomb_energy(orb(p), 's', r_pi) )
        return potential_vector




    def calculate_chi_tensor(self, atomic_coordinates, ndof  ,model_parameters):
        '''Returns the chi tensor for an input list of atomic coordinates'''

        chi_tensor = np.zeros( (ndof,ndof,ndof) )
        for p in range(ndof):
            for orb_q in orbital_types:
                q = ao_index(atom(p), orb_q) # p & q on same atom
                for orb_r in orbital_types:
                    r = ao_index(atom(p), orb_r) # p & r on same atom
                    chi_tensor[p,q,r] = chi_on_atom(orb(p), orb(q), orb(r), model_parameters)

        return chi_tensor



    def chi_on_atom(self, o1, o2, o3):
        '''Returns the value of the chi tensor for 3 orbital indices on the same atom.'''
        if o1 == o2 and o3 == 's':
            return 1.0
        if o1 == o3 and o3 in p_orbitals and o2 == 's':
            return model_parameters['dipole']
        if o2 == o3 and o3 in p_orbitals and o1 == 's':
            return model_parameters['dipole']
        return 0.0


    def hopping_energy(self, o1, o2, r12, model_parameters):
        '''Returns the hopping matrix element for a pair of orbitals of type o1 & o2 separated by a vector r12.'''
        r12_rescaled = r12 / model_parameters['r_hop']
        r12_length = np.linalg.norm(r12_rescaled)
        ans = np.exp( 1.0 - r12_length**2 )
        if o1 == 's' and o2 == 's':
            ans *= model_parameters['t_ss']
        if o1 == 's' and o2 in p_orbitals:
            ans *= np.dot(vec[o2], r12_rescaled) * model_parameters['t_sp']
        if o2 == 's' and o1 in p_orbitals:
            ans *= -np.dot(vec[o1], r12_rescaled)* model_parameters['t_sp']
        if o1 in p_orbitals and o2 in p_orbitals:
            ans *= ( (r12_length**2) * np.dot(vec[o1], vec[o2]) * model_parameters['t_pp2']
                     - np.dot(vec[o1], r12_rescaled) * np.dot(vec[o2], r12_rescaled)
                     * ( model_parameters['t_pp1'] + model_parameters['t_pp2'] ) )
        return ans


    def calculate_hamiltonian_matrix(self, atomic_coordinates, ndof ,model_parameters):
        '''Returns the 1-body Hamiltonian matrix for an input list of atomic coordinates.'''

        hamiltonian_matrix = np.zeros( (ndof,ndof) )
        potential_vector = calculate_potential_vector(atomic_coordinates, model_parameters)
        for p in range(ndof):
            for q in range(ndof):
                if atom(p) != atom(q):
                    r_pq = atomic_coordinates[atom(p)] - atomic_coordinates[atom(q)]
                    hamiltonian_matrix[p,q] = hopping_energy(orb(p), orb(q), r_pq, model_parameters)
                if atom(p) == atom(q):
                    if p == q and orb(p) == 's':
                        hamiltonian_matrix[p,q] += model_parameters['energy_s']
                    if p == q and orb(p) in p_orbitals:
                        hamiltonian_matrix[p,q] += model_parameters['energy_p']
                    for orb_r in orbital_types:
                        r = ao_index(atom(p), orb_r)
                        hamiltonian_matrix[p,q] += ( chi_on_atom(orb(p), orb(q), orb_r, model_parameters)
                                                     * potential_vector[r] )

        return hamiltonian_matrix



    def pseudopotential_energy(self, o, r, model_parameters):
        '''Returns the energy of a pseudopotential between a multipole of type o and an atom separated by a vector r.'''
        ans = model_parameters['v_pseudo']
        r_rescaled = r / model_parameters['r_pseudo']
        r_length = np.linalg.norm(r_rescaled)
        ans *= np.exp( 1.0 - r_length**2 )
        if o in p_orbitals:
            ans *= -2.0 * np.dot(vec[o], r_rescaled)
        return ans


    def calculate_atomic_density_matrix(self, atomic_coordinates, ndof):
        '''Returns a trial 1-electron density matrix for an input list of atomic coordinates.'''

        density_matrix = np.zeros( (ndof,ndof) )
        for p in range(ndof):
            density_matrix[p,p] = orbital_occupation[orb(p)]
        return density_matrix


    def calculate_fock_matrix(self, hamiltonian_matrix, interaction_matrix, density_matrix, chi_tensor):

        '''Returns the Fock matrix defined by the input Hamiltonian, interaction, & density matrices.'''
        fock_matrix = hamiltonian_matrix.copy()
        fock_matrix += 2.0*np.einsum('pqt,rsu,tu,rs',
                                     chi_tensor, chi_tensor, interaction_matrix, density_matrix, optimize=True)
        fock_matrix -= np.einsum('rqt,psu,tu,rs',
                                 chi_tensor, chi_tensor, interaction_matrix, density_matrix, optimize=True)
        return fock_matrix


    def calculate_density_matrix(self, fock_matrix):
        '''Returns the 1-electron density matrix defined by the input Fock matrix.'''
    
        num_occ = (ionic_charge//2)*np.size(fock_matrix,0)//orbitals_per_atom
        orbital_energy, orbital_matrix = np.linalg.eigh(fock_matrix)
        occupied_matrix = orbital_matrix[:,:num_occ]
        density_matrix = occupied_matrix @ occupied_matrix.T

        return density_matrix

    def scf_cycle(self, hamiltonian_matrix, interaction_matrix, density_matrix, chi_tensor, 
                max_scf_iterations = 100, mixing_fraction = 0.25, convergence_tolerance = 1e-4):
        '''Returns converged density & Fock matrices defined by the input Hamiltonian, interaction, & density matrices.'''
        old_density_matrix = density_matrix.copy()
        for iteration in range(max_scf_iterations):
            new_fock_matrix = calculate_fock_matrix(hamiltonian_matrix, interaction_matrix, old_density_matrix, chi_tensor)
            new_density_matrix = calculate_density_matrix(new_fock_matrix)

            error_norm = np.linalg.norm( old_density_matrix - new_density_matrix )
            if error_norm < convergence_tolerance:
                return new_density_matrix, new_fock_matrix

            old_density_matrix = (mixing_fraction * new_density_matrix
                                  + (1.0 - mixing_fraction) * old_density_matrix)
        print("WARNING: SCF cycle didn't converge")
        return new_density_matrix, new_fock_matrix


    def calculate_energy_ion(self, atomic_coordinates):
        '''Returns the ionic contribution to the total energy for an input list of atomic coordinates.'''
        energy_ion = 0.0
        for i, r_i in enumerate(atomic_coordinates):
            for j, r_j in enumerate(atomic_coordinates):
                if i < j:
                    energy_ion += (ionic_charge**2)*coulomb_energy('s', 's', r_i - r_j)
        return energy_ion


    def calculate_energy_scf(self, hamiltonian_matrix, fock_matrix, density_matrix):
        '''Returns the Hartree-Fock total energy defined by the input Hamiltonian, Fock, & density matrices.'''
        energy_scf = np.einsum('pq,pq',hamiltonian_matrix + fock_matrix,density_matrix)
        return energy_scf



