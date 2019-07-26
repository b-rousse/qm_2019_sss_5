import numpy as np
import mp2 as mp2
import hartree_fock as hf
import noble_gas_model as noble_gas_model

if __name__ == "__main__":
    NobleGasModel = noble_gas_model.NobleGasModl()
    atomic_coordinates = np.array([[0.0,0.0,0.0], [3.0,4.0,5.0]])
    hartree_fock_instance = hf.HartreeFock(NobleGasModel, atomic_coordinates)
    hartree_fock_instance.density_matrix = hartree_fock_instance.calculate_atomic_density_matrix(NobleGasModel)
    hartree_fock_instance.density_matrix, hartree_fock_instance.fock_matrix = hartree_fock_instance.scf_cycle(NobleGasModel)
    energy_scf = hartree_fock_instance.calculate_energy_scf()
    energy_ion = hartree_fock_instance.calculate_energy_ion(NobleGasModel)
    print(F'The SCF energy is  {energy_scf} and the ion energy is {energy_ion} ')
    #mp2_instance = mp2.MP2(hartree_fock_instance)
    #print(F'The MP2 energy is {mp2.MP2.calculate_energy_mp2}')