import numpy as np
import mp2_no_hf as mp2_no_hf
import hartree_fock as hf
import noble_gas_model as noble_gas_model

if __name__ == "__main__":
    NobleGasModel = noble_gas_model.NobleGasModl()
    atomic_coordinates = np.array([[0.0,0.0,0.0], [3.0,4.0,5.0]])
    hartree_fock_instance = hf.HartreeFock(NobleGasModel, atomic_coordinates)
    hartree_fock_instance.density_matrix = hartree_fock_instance.calculate_atomic_density_matrix(NobleGasModel)
    hartree_fock_instance.density_matrix, hartree_fock_instance.fock_matrix = hartree_fock_instance.scf_cycle(NobleGasModel)
    hartree_fock_instance.energy_scf = hartree_fock_instance.calculate_energy_scf()
    hartree_fock_instance.energy_ion = hartree_fock_instance.calculate_energy_ion(NobleGasModel)
    print(F'The SCF energy is  {hartree_fock_instance.energy_scf} and the ion energy is {hartree_fock_instance.energy_ion} ')
    #mp2_instance = mp2.MP2(NobleGasModel,atomic_coordinates)
    #mp2_instance.mp2_energy = mp2_instance.calculate_energy_mp2()
    #print(F'The MP2 energy is {mp2_instance.mp2_energy}')
    mp2_standalone_instance = mp2_no_hf.MP2NoHF(hartree_fock_instance, NobleGasModel)
    mp2_iso_energy = mp2_standalone_instance.calculate_energy_mp2()
