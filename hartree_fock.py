import numpy as np
#from noble_gas_model import NobleGasModl
class HartreeFock:
    def __init__(self, NobleGasModel, atomic_coordinates):
        self.atomic_coordinates = atomic_coordinates
        #self.gas_model = NobleGasModel...
        self.ndof = len(self.atomic_coordinates) * NobleGasModel.orbitals_per_atom
        self.interaction_matrix = self.calculate_interaction_matrix(NobleGasModel)
        self.hamiltonian_matrix = self.calculate_hamiltonian_matrix(NobleGasModel)
        self.chi_tensor = self.calculate_chi_tensor(NobleGasModel)
        self.density_matrix = self.calculate_atomic_density_matrix(NobleGasModel)
        self.fock_matrix = self.calculate_fock_matrix() #check for error
        #self.density_matrix = self.calculate_density_matrix(self.fock_matrix)
        #self.chi_tensor = self.calculate_chi_tensor(self.atomic_coordinates, self.ndof ,NobleGasModel.model_parameters) 

    def calculate_interaction_matrix(self, NobleGasModel):

        '''Returns the electron-electron interaction energy matrix for an input list of atomic coordinates.'''
        interaction_matrix = np.zeros( (self.ndof,self.ndof) )
        for p in range(self.ndof):
            for q in range(self.ndof):
                if NobleGasModel.atom(p) != NobleGasModel.atom(q):
                    r_pq = self.atomic_coordinates[NobleGasModel.atom(p)] - self.atomic_coordinates[NobleGasModel.atom(q)]
                    interaction_matrix[p,q] = self.calculate_coulomb_energy(NobleGasModel.orb(p), NobleGasModel.orb(q), r_pq, NobleGasModel)
                if p == q and NobleGasModel.orb(p) == 's':
                    interaction_matrix[p,q] = NobleGasModel.model_parameters['coulomb_s']
                if p == q and NobleGasModel.orb(p) in NobleGasModel.p_orbitals:
                    interaction_matrix[p,q] = NobleGasModel.model_parameters['coulomb_p']

        return interaction_matrix


    def calculate_coulomb_energy(self, o1, o2, r12, NobleGasModel):
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
        if o1 == 's' and o2 in NobleGasModel.p_orbitals:
            coulomb_energy = np.dot(NobleGasModel.vec[o2], r12) / r12_length**3
        if o2 == 's' and o1 in NobleGasModel.p_orbitals:
            coulomb_energy = -1*np.dot(NobleGasModel.vec[o1], r12) / r12_length**3
        if o1 in NobleGasModel.p_orbitals and o2 in NobleGasModel.p_orbitals:
            coulomb_energy = ( np.dot(NobleGasModel.vec[o1], NobleGasModel.vec[o2]) / r12_length**3
                   - 3.0 * np.dot(NobleGasModel.vec[o1], r12) * np.dot(NobleGasModel.vec[o2], r12) / r12_length**5 )
        return coulomb_energy



    def calculate_potential_vector(self, NobleGasModel):#maybe make noblegasparamws part of self.

        '''Returns the electron-ion potential energy vector for an input list of atomic coordinates.'''

        potential_vector = np.zeros(self.ndof)
        for p in range(self.ndof):
            potential_vector[p] = 0.0
            for atom_i,r_i in enumerate(self.atomic_coordinates):
                r_pi = self.atomic_coordinates[NobleGasModel.atom(p)] - r_i
                if atom_i != NobleGasModel.atom(p):
                    potential_vector[p] += (self.calculate_pseudopotential_energy(NobleGasModel.orb(p), r_pi, NobleGasModel) - NobleGasModel.ionic_charge * self.calculate_coulomb_energy(NobleGasModel.orb(p), 's', r_pi, NobleGasModel) )
        return potential_vector




    def calculate_chi_tensor(self, NobleGasModel):
        '''Returns the chi tensor for an input list of atomic coordinates'''

        chi_tensor = np.zeros((self.ndof,self.ndof,self.ndof))
        for p in range(self.ndof):
            for orb_q in NobleGasModel.orbital_types:
                q = NobleGasModel.ao_index(NobleGasModel.atom(p), orb_q) # p & q on same atom
                for orb_r in NobleGasModel.orbital_types:
                    r = NobleGasModel.ao_index(NobleGasModel.atom(p), orb_r) # p & r on same atom
                    chi_tensor[p,q,r] = self.chi_on_atom(NobleGasModel.orb(p), NobleGasModel.orb(q), NobleGasModel.orb(r), NobleGasModel)

        return chi_tensor



    def chi_on_atom(self, o1, o2, o3, NobleGasModel):
        '''Returns the value of the chi tensor for 3 orbital indices on the same atom.'''
        if o1 == o2 and o3 == 's':
            return 1.0
        if o1 == o3 and o3 in NobleGasModel.p_orbitals and o2 == 's':
            return NobleGasModel.model_parameters['dipole']
        if o2 == o3 and o3 in NobleGasModel.p_orbitals and o1 == 's':
            return NobleGasModel.model_parameters['dipole']
        return 0.0


    def calculate_hopping_energy(self, o1, o2, r12, NobleGasModel):
        '''Returns the hopping matrix element for a pair of orbitals of type o1 & o2 separated by a vector r12.'''
        r12_rescaled = r12 / NobleGasModel.model_parameters['r_hop']
        r12_length = np.linalg.norm(r12_rescaled)
        ans = np.exp( 1.0 - r12_length**2 )
        if o1 == 's' and o2 == 's':
            ans *= NobleGasModel.model_parameters['t_ss']
        if o1 == 's' and o2 in NobleGasModel.p_orbitals:
            ans *= np.dot(NobleGasModel.vec[o2], r12_rescaled) * NobleGasModel.model_parameters['t_sp']
        if o2 == 's' and o1 in NobleGasModel.p_orbitals:
            ans *= -1*np.dot(NobleGasModel.vec[o1], r12_rescaled)* NobleGasModel.model_parameters['t_sp']
        if o1 in NobleGasModel.p_orbitals and o2 in NobleGasModel.p_orbitals:
            ans *= ( (r12_length**2) * np.dot(NobleGasModel.vec[o1], NobleGasModel.vec[o2]) * NobleGasModel.model_parameters['t_pp2']
                     - np.dot(NobleGasModel.vec[o1], r12_rescaled) * np.dot(NobleGasModel.vec[o2], r12_rescaled)
                     * ( NobleGasModel.model_parameters['t_pp1'] + NobleGasModel.model_parameters['t_pp2'] ) )
        return ans


    def calculate_hamiltonian_matrix(self, NobleGasModel):
        '''Returns the 1-body Hamiltonian matrix for an input list of atomic coordinates.'''

        hamiltonian_matrix = np.zeros((self.ndof,self.ndof))
        potential_vector = self.calculate_potential_vector(NobleGasModel)
        for p in range(self.ndof):
            for q in range(self.ndof):
                if NobleGasModel.atom(p) != NobleGasModel.atom(q):
                    r_pq = self.atomic_coordinates[NobleGasModel.atom(p)] - self.atomic_coordinates[NobleGasModel.atom(q)]
                    hamiltonian_matrix[p,q] = self.calculate_hopping_energy(NobleGasModel.orb(p), NobleGasModel.orb(q), r_pq, NobleGasModel)
                if NobleGasModel.atom(p) == NobleGasModel.atom(q):
                    if p == q and NobleGasModel.orb(p) == 's':
                        hamiltonian_matrix[p,q] += NobleGasModel.model_parameters['energy_s']
                    if p == q and NobleGasModel.orb(p) in NobleGasModel.p_orbitals:
                        hamiltonian_matrix[p,q] += NobleGasModel.model_parameters['energy_p']
                    for orb_r in NobleGasModel.orbital_types:
                        r = NobleGasModel.ao_index(NobleGasModel.atom(p), orb_r)
                        hamiltonian_matrix[p,q] += ( self.chi_on_atom(NobleGasModel.orb(p), NobleGasModel.orb(q), orb_r, NobleGasModel) * potential_vector[r] )

        return hamiltonian_matrix



    def calculate_pseudopotential_energy(self, o, r, NobleGasModel):
        '''Returns the energy of a pseudopotential between a multipole of type o and an atom separated by a vector r.'''
        ans = NobleGasModel.model_parameters['v_pseudo']
        r_rescaled = r / NobleGasModel.model_parameters['r_pseudo']
        r_length = np.linalg.norm(r_rescaled)
        ans *= np.exp( 1.0 - r_length**2 )
        if o in NobleGasModel.p_orbitals:
            ans *= -2.0 * np.dot(NobleGasModel.vec[o], r_rescaled)
        return ans


    def calculate_atomic_density_matrix(self, NobleGasModel):
        '''Returns a trial 1-electron density matrix for an input list of atomic coordinates.'''

        density_matrix = np.zeros((self.ndof,self.ndof))
        for p in range(self.ndof):
            density_matrix[p,p] = NobleGasModel.orbital_occupations[NobleGasModel.orb(p)]
        return density_matrix


    def calculate_fock_matrix(self):
        '''Returns the Fock matrix defined by the input Hamiltonian, interaction, & density matrices.'''
        fock_matrix = self.hamiltonian_matrix.copy()
        fock_matrix += 2.0*np.einsum('pqt,rsu,tu,rs', self.chi_tensor, self.chi_tensor, self.interaction_matrix, self.density_matrix, optimize=True)
        fock_matrix -= np.einsum('rqt,psu,tu,rs', self.chi_tensor, self.chi_tensor, self.interaction_matrix, self.density_matrix, optimize=True)
        return fock_matrix


    def calculate_density_matrix(self, NobleGasModel):
        '''Returns the 1-electron density matrix defined by the input Fock matrix.'''
    
        num_occ = (NobleGasModel.ionic_charge//2)*np.size(self.fock_matrix,0)//NobleGasModel.orbitals_per_atom
        orbital_energy, orbital_matrix = np.linalg.eigh(self.fock_matrix)
        occupied_matrix = orbital_matrix[:,:num_occ]
        density_matrix = occupied_matrix @ occupied_matrix.T

        return density_matrix

    def scf_cycle(self, NobleGasModel, max_scf_iterations = 200, mixing_fraction = 0.25, convergence_tolerance = 1e-7):
        '''Returns converged density & Fock matrices defined by the input Hamiltonian, interaction, & density matrices.'''
        old_density_matrix = self.density_matrix.copy()
        for iteration in range(max_scf_iterations):
            print(F'iteration is {iteration}')
            self.fock_matrix = self.calculate_fock_matrix()
            self.density_matrix = self.calculate_density_matrix(NobleGasModel)

            error_norm = np.linalg.norm(old_density_matrix - self.density_matrix)
            print(F'Error norm is {error_norm}')
            if error_norm < convergence_tolerance:
                return self.density_matrix, self.fock_matrix

            old_density_matrix = (mixing_fraction * self.density_matrix + (1.0 - mixing_fraction) * old_density_matrix)
        print("WARNING: SCF cycle didn't converge")
        return self.density_matrix, self.fock_matrix


    def calculate_energy_ion(self, NobleGasModel):
        '''Returns the ionic contribution to the total energy for an input list of atomic coordinates.'''
        energy_ion = 0.0
        for i, r_i in enumerate(self.atomic_coordinates):
            for j, r_j in enumerate(self.atomic_coordinates):
                if i < j:
                    energy_ion += (NobleGasModel.ionic_charge**2)*self.calculate_coulomb_energy('s', 's', r_i - r_j, NobleGasModel)
        return energy_ion


    def calculate_energy_scf(self):
        '''Returns the Hartree-Fock total energy defined by the input Hamiltonian, Fock, & density matrices.'''
        energy_scf = np.einsum('pq,pq',self.hamiltonian_matrix + self.fock_matrix,self.density_matrix)
        return energy_scf



