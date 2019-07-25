import numpy as np

class MP2XTREME:
    def __init__(self, myHF, myMolecule):
        self.hf_energy = myHF.hf_energy
        self.occupied_energy = myHF.occupied_energy
        self.interaction_tensor = myHF.interaction_tensor
        self.myMolecule = myMolecule
        self.mp2_energy = 0.0
        self.mp2_correction = 0.0

    def partition_orbitals(self, fock_matrix):
        '''Returns a list with the occupied/virtual energies & orbitals defined by the input Fock matrix.'''
        num_occ = (self.myMolecule.ionic_charge // 2) * np.size(fock_matrix,
                                                0) // self.myMolecule.orbitals_per_atom
        orbital_energy, orbital_matrix = np.linalg.eigh(fock_matrix)
        occupied_energy = orbital_energy[:num_occ]
        virtual_energy = orbital_energy[num_occ:]
        occupied_matrix = orbital_matrix[:, :num_occ]
        virtual_matrix = orbital_matrix[:, num_occ:]

        return occupied_energy, virtual_energy, occupied_matrix, virtual_matrix

    def transform_interaction_tensor(self, occupied_matrix, virtual_matrix,
                                    interaction_matrix, chi_tensor):
        '''Returns a transformed V tensor defined by the input occupied, virtual, & interaction matrices.'''
        chi2_tensor = np.einsum('qa,ri,qrp',
                                virtual_matrix,
                                occupied_matrix,
                                chi_tensor,
                                optimize=True)
        interaction_tensor = np.einsum('aip,pq,bjq->aibj',
                                    chi2_tensor,
                                    interaction_matrix,
                                    chi2_tensor,
                                    optimize=True)
        return interaction_tensor

    def calculate_energy_mp2(self, fock_matrix, interaction_matrix, chi_tensor):
        '''Returns the MP2 contribution to the total energy defined by the input Fock & interaction matrices.'''
        E_occ, E_virt, occupied_matrix, virtual_matrix = self.partition_orbitals(
            fock_matrix)
        V_tilde = self.transform_interaction_tensor(occupied_matrix, virtual_matrix,
                                            interaction_matrix, chi_tensor)

        energy_mp2 = 0.0
        num_occ = len(E_occ)
        num_virt = len(E_virt)
        for a in range(num_virt):
            for b in range(num_virt):
                for i in range(num_occ):
                    for j in range(num_occ):
                        energy_mp2 -= (
                            (2.0 * V_tilde[a, i, b, j]**2 -
                            V_tilde[a, i, b, j] * V_tilde[a, j, b, i]) /
                            (E_virt[a] + E_virt[b] - E_occ[i] - E_occ[j]))
        self.mp2_correction = energy_mp2
        self.mp2_energy = self.hf_energy + energy_mp2
        return energy_mp2
