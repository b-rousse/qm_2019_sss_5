import numpy as np
import mp2 as mp2
import hartree_fock as hf
import noble_gas_model as noble_gas_model

if __name__ == "__main__":
    noble_gas_instance = noble_gas_model.NobleGasModel()
    atomic_coordinates = np.array([[0.0,0.0,0.0], [3.0,4.0,5.0]])
    hartree_fock_instance = hf.HartreeFock(noble_gas_instance, atomic_coordinates)
    hartree_fock_instance.density_matrix = hartree_fock_instance.calculate_atomic_density_matrix(NobleGasModel)
    hartree_fock_instance.density_matrix = hartree_fock_instance.canc

