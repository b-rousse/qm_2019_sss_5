import numpy as np
from hartree_fock import HartreeFock
from noble_gas_model import NobleGasModl
"""
This class takes an instance of HartreeFock and NobelGas and is able to calculate the MP2 energy correction on top of the HF energy
"""
class MP2NoHF:
    def __init__(self, HartreeFock, NobleGasModel):#here I should only pass what I need for MP2.
        self.hf_energy = HartreeFock.energy_scf
        self.NobleGasModel = NobleGasModel
        self.mp2_energy = 0.0
        self.fock_matrix = HartreeFock.fock_matrix
        self.interaction_matrix = HartreeFock.interaction_matrix
        self.chi_tensor = HartreeFock.chi_tensor

    """
    This function partitions the converged molecular orbitals of the HF calculation into occupied and virtual orbitals,
    and does the same for their orbital energies.

    Parameters
    ----------
    fock_matrix : np.array
        A np.array of size (num_ao, num_ao). This is the final converged fock matrix obtained from the HF calculation
    
    Returns
    -------
    occupied_energy : np.array
        A np.array of size (num_occ) which contains the energies of the occupied orbitals, ordered in increasing energy.
    virtual_energy : np.array
        A np.array of size (num_ao - num_occ) which contains the energies of the virtual orbitals, ordered in increasing energy.
    occupied_matrix : np.array
        A np.array of size (num_ao, num_occ) where each column is a canonical occupied molecular orbital.
    virtual_matrix : np.array
        A np.array of size (num_ao, num_ao - num_occ) where each column is a canonical virtual molecular orbital.
    """
    def partition_orbitals(self):
        '''Returns a list with the occupied/virtual energies & orbitals defined by the input Fock matrix.'''
        num_occ = (self.NobleGasModel.ionic_charge // 2) * np.size(self.fock_matrix, 0) // self.NobleGasModel.orbitals_per_atom
        orbital_energy, orbital_matrix = np.linalg.eigh(self.fock_matrix)
        occupied_energy = orbital_energy[:num_occ]
        virtual_energy = orbital_energy[num_occ:]
        occupied_matrix = orbital_matrix[:, :num_occ]
        virtual_matrix = orbital_matrix[:, num_occ:]

        return occupied_energy, virtual_energy, occupied_matrix, virtual_matrix

    """
    This function rotates the interaction tensor into the basis of canonical molecular orbitals.

    Parameters
    ----------
    occupied_matrix : np.array
        A np.array of size (num_ao, num_occ) where each column is a canonical occupied molecular orbital.
    virtual_matrix : np.array
        A np.array of size (num_ao, num_ao - num_occ) where each column is a canonical virtual molecular orbital.
    interaction_matrix : np.array
        A np.array of size (num_ao, num_ao) representing the interaction between each atomic orbital in the AO basis
    chi_tensor : np.array
        A np.array of size (num_ao, num_ao, num_ao) representing the three electron integrals, I guess. In AO basis.

    Returns
    -------
    interaction_tensor : np.array
        A np.array of size (num_ao, num_ao, num_ao, num_ao). The 2 electron integrals in the canonical basis.
    """ 
    def transform_interaction_tensor(self, occupied_matrix, virtual_matrix, interaction_matrix, chi_tensor):
        '''Returns a transformed V tensor defined by the input occupied, virtual, & interaction matrices.'''
        chi2_tensor = np.einsum('qa,ri,qrp', virtual_matrix, occupied_matrix, chi_tensor, optimize=True)
        interaction_tensor = np.einsum('aip,pq,bjq->aibj', chi2_tensor, interaction_matrix, chi2_tensor, optimize=True)
        return interaction_tensor

    """
    Returns the MP2 contribution to the total energy defined by the input Fock & interaction matrices.

    Parameters
    ----------
    fock_matrix : np.array
        A np.array of size (num_ao, num_ao). This is the final converged fock matrix obtained from the HF calculation in AO basis.
    interaction_matrix : np.array
        A np.array of size (num_ao, num_ao) representing the interaction between each atomic orbital in the AO basis
    chi_tensor : np.array
        A np.array of size (num_ao, num_ao, num_ao) representing the three electron integrals, I guess. In AO basis.

    Returns
    -------
    energy_mp2 : float
        MP2 contribution to the total energy.
    """ 
    def calculate_energy_mp2(self, fock_matrix = None, interaction_matrix = None, chi_tensor = None):
        if fock_matrix == None:
            fock_matrix = self.fock_matrix
        if interaction_matrix == None:
            interaction_matrix = self.interaction_matrix
        if chi_tensor == None:
            chi_tensor = self.chi_tensor
        E_occ, E_virt, occupied_matrix, virtual_matrix = self.partition_orbitals()
        V_tilde = self.transform_interaction_tensor(occupied_matrix, virtual_matrix, interaction_matrix, chi_tensor)

        energy_mp2 = 0.0
        num_occ = len(E_occ)
        num_virt = len(E_virt)
        for a in range(num_virt):
            for b in range(num_virt):
                for i in range(num_occ):
                    for j in range(num_occ):
                        energy_mp2 -= (
                            (2.0 * V_tilde[a, i, b, j]**2 - V_tilde[a, i, b, j] * V_tilde[a, j, b, i])
                             / (E_virt[a] + E_virt[b] - E_occ[i] - E_occ[j]))
        self.mp2_correction = energy_mp2
        self.mp2_energy = self.hf_energy + energy_mp2
        return energy_mp2
