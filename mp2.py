import numpy as np
from hartree_fock import HartreeFock
from noble_gas_model import NobleGasModl

class MP2(HartreeFock):
    def __init__(self, NobleGasModel, atomic_coordinates):
        super().__init__(NobleGasModel, atomic_coordinates)
        self.orbitals_per_atom = NobleGasModel.orbitals_per_atom
        self.ionic_charge = NobleGasModel.ionic_charge
        self.energy_mp2 = 0.0


    def partition_orbitals(self):
        '''Returns a list with the occupied/virtual energies & orbitals defined by the input Fock matrix.'''
        num_occ = (self.ionic_charge // 2) * np.size(self.fock_matrix,0) // self.orbitals_per_atom
        orbital_energy, orbital_matrix = np.linalg.eigh(self.fock_matrix)
        occupied_energy = orbital_energy[:num_occ]
        virtual_energy = orbital_energy[num_occ:]
        occupied_matrix = orbital_matrix[:, :num_occ]
        virtual_matrix = orbital_matrix[:, num_occ:]

        return occupied_energy, virtual_energy, occupied_matrix, virtual_matrix

    def transform_interaction_tensor(self, occupied_matrix, virtual_matrix):
        '''Returns a transformed V tensor defined by the input occupied, virtual, & interaction matrices.'''
        chi2_tensor = np.einsum('qa,ri,qrp', virtual_matrix, occupied_matrix, self.chi_tensor, optimize=True)
        interaction_tensor = np.einsum('aip,pq,bjq->aibj', chi2_tensor, self.interaction_matrix, chi2_tensor, optimize=True)
        return interaction_tensor

    def calculate_energy_mp2(self):
        '''Returns the MP2 contribution to the total energy defined by the input Fock & interaction matrices.'''
        E_occ, E_virt, occupied_matrix, virtual_matrix = self.partition_orbitals()
        V_tilde = self.transform_interaction_tensor(occupied_matrix, virtual_matrix)

        energy_mp2 = 0.0
        num_occ = len(E_occ)
        num_virt = len(E_virt)
        for a in range(num_virt):
            for b in range(num_virt):
                for i in range(num_occ):
                    for j in range(num_occ):
                        energy_mp2 -= ((2.0 * V_tilde[a, i, b, j]**2 - V_tilde[a, i, b, j] * V_tilde[a, j, b, i])
                         / (E_virt[a] + E_virt[b] - E_occ[i] - E_occ[j]))
        return energy_mp2
    
    def __str__(self):
        return 'MP2 correction : ' + str(self.energy_mp2)
