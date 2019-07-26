
![](qm5.png)


## Description
qm5 Simulates a cluster of Argon atoms using quantum mechanics (QM). The codes uses simple semiempirical quantum mechanics. Since Argon is a noble gas, mostly London dispersion forces predominate. The lowest energy dipole transition is from the occupied $3p$ states to the unoccupied 4s steate, including these 4 atomic orbitals per atom.  


## Installation

Use the conda manager

```bash
conda install qm5-molssi -c conda-forge
```

## Usage

The main file is `routine.py` that performs as follows:

```python

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

```

qm5 works with 3 main files: `noble_gas_model.py`, `HartreeFock.py`, and `MP2.py`. Each of these files is a class that has attributes and methods associated with the class. All of them work together to produce the Hartree Fock energy with the MP2 correct    ion for a Noble Gas, e.g. Argon.



## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
