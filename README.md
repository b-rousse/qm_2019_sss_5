
![](qm5.png)


## Description
qm5 Simulates a cluster of Argon atoms using quantum mechanics (QM). The codes uses simple semiempirical quantum mechanics. Since Argon is a noble gas, mostly London dispersion forces predominate. The lowest energy dipole transition is from the occupied $3p$ states to the unoccupied 4s steate, including these 4 atomic orbitals per atom.  

$$ \begin{bmatrix} x & \dot{x} & \theta & \dot{\theta} & L & m & M \end{bmatrix} $$


## Installation

Use the conda manager

```bash
conda install qm5-molssi -c conda-forge
```

## Usage

qm5 works with 3 main files: `noble_gas_model`, `HartreeFock`, and `MP2`. Each of these files is a class that has attributes and methods associated with the class. All of them work together to produce the Hartree Fock energy with the MP2 correction for a Noble Gas, e.g. Argon. 


```python
import qm5

my_argon = qm5.Noble_Gas()
my_hf = qm5.HartreeFock(my_argon, user_defined_atomic_coordinates)
my_MP2 = qm5.MP2(my_hf)


```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
