import numpy as np
import mp2 as mp2
import hartree_fock as hf
import noble_gas_model as noble_gas_model

if __name__ == "__main__":
    NobleGasModel = noble_gas_model.NobleGasModl()
    atomic_coordinates = np.array([[0.0,0.0,0.0], [3.0,4.0,5.0]])
    mp2_instance = mp2.MP2(NobleGasModel,atomic_coordinates)
    mp2_instance.density_matrix = mp2_instance.calculate_atomic_density_matrix(NobleGasModel)
    mp2_instance.density_matrix, mp2_instance.fock_matrix = mp2_instance.scf_cycle(NobleGasModel)
    mp2_instance.energy_scf = mp2_instance.calculate_energy_scf()
    mp2_instance.energy_ion = mp2_instance.calculate_energy_ion(NobleGasModel)
    print(F'The SCF energy is  {mp2_instance.energy_scf} and the ion energy is {mp2_instance.energy_ion} ')
    mp2_instance.mp2_energy = mp2_instance.calculate_energy_mp2()
    print(F'The MP2 energy is {mp2_instance.mp2_energy}')
