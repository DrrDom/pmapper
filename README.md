# Pmapper - 3D pharmacophore signatures and fingerprints

Pmapper is a Python module to generate 3D pharmacophore signatures and fingerprints.
Signatures uniquely encode 3D pharmacophores with hashes suitable for fast identification of identical pharmacophores.

## Dependency

`rdkit >= 2017.09`  
`networkx >= 1.11`

## Installation
```text
pip install pmapper
```

## Examples

### Load modules
```python
from pmapper.pharmacophore import Pharmacophore as P
from rdkit import Chem
from rdkit.Chem import AllChem
from pprint import pprint
```
### Create pharmacophore from a single conformer using default feature definitions
```python
# load a molecule from SMILES and generate 3D coordinates
mol = Chem.MolFromSmiles('C1CC(=O)NC(=O)C1N2C(=O)C3=CC=CC=C3C2=O')  # talidomide
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)

# create pharmacophore
p = P()
p.load_from_mol(mol)
```
### Get 3D pharmacophore signature
```python
# get 3D pharmacophore signature
sig = p.get_signature_md5()
print(sig)
```
Output:
```text
98504647beeb143ae50bb6b7798ca0f0
```
### Get 3D pharmacophore signature with non-zero tolerance
```python
sig = p.get_signature_md5(tol=5)
print(sig)
```
Output:
```text
bc54806ba01bf59736a7b62b017d6e1d
```
### Create pharmacophores for a multiple conformer compound
```python
from pmapper.utils import load_multi_conf_mol

# create multiple conformer molecule
AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=1024)

ps = load_multi_conf_mol(mol)

sig = [p.get_signature_md5() for p in ps]

pprint(sig)  # identical signatures occur
```
Output:
```text
['d5f5f9d65e39cb8605f1fa9db5b2fbb0',
 '6204791002d1e343b2bde323149fa780',
 'abfabd8a4fcf5719ed6bf2c71a60852c',
 'dfe9f17d30210cb94b8dd7acf77feae9',
 'abfabd8a4fcf5719ed6bf2c71a60852c',
 'e739fb5f9985ce0c65a16da41da4a33f',
 '2297ddf0e437b7fc32077f75e3924dcd',
 'e739fb5f9985ce0c65a16da41da4a33f',
 '182a00bd9057abd0c455947d9cfa457c',
 '68f226d474808e60ab1256245f64c2b7']
```
Identical hashes should correspond to pharmacophores with low RMSD. Pharmacophores #2 and #4 have identical hash `abfabd8a4fcf5719ed6bf2c71a60852c`. Let's check RMSD.
```python
from pmapper.utils import get_rms
for i in range(len(ps)):
    print("rmsd bewteen 2 and %i pharmacophore:" % i, round(get_rms(ps[2], ps[i]), 2))
```
Output
```text
rmsd bewteen 2 and 0 pharmacophore: 0.63
rmsd bewteen 2 and 1 pharmacophore: 0.99
rmsd bewteen 2 and 2 pharmacophore: 0.0
rmsd bewteen 2 and 3 pharmacophore: 0.41
rmsd bewteen 2 and 4 pharmacophore: 0.18
rmsd bewteen 2 and 5 pharmacophore: 0.19
rmsd bewteen 2 and 6 pharmacophore: 1.15
rmsd bewteen 2 and 7 pharmacophore: 0.32
rmsd bewteen 2 and 8 pharmacophore: 0.69
rmsd bewteen 2 and 9 pharmacophore: 0.36
```
They really have RMSD < binning step (1A by default). However, other pharmacophores with distinct hashes also have low RMSD to #2. Identical hashes guarantee low RMSD between corresponding pharmacophores, but not vice versa.

### Pharmacophore match
Create a two-point pharmacophore model and match with a pharmacophore of a molecule (both pharmacophores should have identical binning steps)
```python
q = P()
q.load_from_feature_coords([('a', (3.17, -0.23, 0.24)), ('D', (-2.51, -1.28, -1.14))])
p.fit_model(q)
```
Output
```text
(0, 1)
```
If they do not match `None` will be returned

### Generate 3D pharmacophore fingerprint
```python
# generate 3D pharmacophore fingerprint which takes into account stereoconfiguration
b = p.get_fp(min_features=4, max_features=4)   # set of activated bits
print(b)
```
Output (a set of activated bit numbers):
```text
{259, 1671, 521, 143, 912, 402, 278, 406, 1562, 1692, 1835, 173, 558, 1070, 942, 1202, 1845, 823, 1476, 197, 968, 1355, 845, 1741, 1364, 87, 1881, 987, 1515, 378, 628, 1141, 1401, 1146, 2043}
```
Change settings:
```python
b = p.get_fp(min_features=4, max_features=4, nbits=4096, activate_bits=2)
print(b)
```
Output (a set of activated bit numbers):
```text
{389, 518, 2821, 1416, 2952, 395, 3339, 511, 3342, 1937, 1042, 2710, 1817, 1690, 3482, 3737, 286, 1824, 1700, 804, 1318, 2729, 3114, 812, 556, 175, 3763, 2356, 3124, 1077, 1975, 3384, 1081, 185, 65, 1223, 713, 1356, 1998, 1487, 2131, 85, 3670, 1877, 3030, 2395, 1116, 2141, 1885, 347, 2404, 1382, 1257, 3049, 2795, 3691, 2541, 1646, 2283, 241, 113, 3698, 756, 2548, 4086, 2293, 1528, 2802, 127}
```
### Save/load pharmacophore
```python
p.save_to_pma('filename.pma')
```
Output is a text file having json format.
```python
p = P()
p.load_from_pma('filename.pma')
```
### Support other formats
Pharmacophores can be saved/loaded from LigandScout pml-files. Also pharmacophores can be read from xyz-files.

### Caching
Pharmacophores can be created with enabled `cache` argument. This will speed up all futher repeated calls to retrive hash, fingerprints or descriptors.
```python
p = P(cache=True)
```

## Speed tests
Generation of pharmacophore signatures (hashes) is a CPU-bound task. The computation speed depends on the number of features in pharmacophores.  
Tests were run on 500 compounds (a random subset from Drugbank). Up to 50 conformers were generated for each compound. Up to 100 pharmacophores having a particular number of features were chosen randomly from the whole number of 25000 pharmacophores to generate pharmacophore signatures. Cache was disabled but enabled cache would not affect calculation time for the first time function calls.  
Laptop configuration:
- Intel(R) Core(TM) i7-5500U CPU @ 2.40GHz
- 12 GB RAM
- calculation was run in 1 thread (the module is thread safe and calculations can be parallelized)

```text
pharmacophore generation: 19.21 s
total number of pharmacophores: 25000

pharmacophore hash generation:
50 pharmacophores having 2 features: 0.00 s; time per pharmacophore: 0.00000 s
100 pharmacophores having 3 features: 0.01 s; time per pharmacophore: 0.00010 s
100 pharmacophores having 4 features: 0.01 s; time per pharmacophore: 0.00010 s
100 pharmacophores having 5 features: 0.04 s; time per pharmacophore: 0.00040 s
100 pharmacophores having 6 features: 0.12 s; time per pharmacophore: 0.00120 s
100 pharmacophores having 7 features: 0.24 s; time per pharmacophore: 0.00240 s
100 pharmacophores having 8 features: 0.51 s; time per pharmacophore: 0.00510 s
100 pharmacophores having 9 features: 0.94 s; time per pharmacophore: 0.00940 s
100 pharmacophores having 10 features: 1.86 s; time per pharmacophore: 0.01860 s
100 pharmacophores having 11 features: 3.02 s; time per pharmacophore: 0.03020 s
100 pharmacophores having 12 features: 4.17 s; time per pharmacophore: 0.04170 s
100 pharmacophores having 13 features: 7.04 s; time per pharmacophore: 0.07040 s
100 pharmacophores having 14 features: 9.29 s; time per pharmacophore: 0.09290 s
100 pharmacophores having 15 features: 12.94 s; time per pharmacophore: 0.12940 s
100 pharmacophores having 16 features: 17.79 s; time per pharmacophore: 0.17790 s
100 pharmacophores having 17 features: 23.58 s; time per pharmacophore: 0.23580 s
100 pharmacophores having 18 features: 33.83 s; time per pharmacophore: 0.33830 s
100 pharmacophores having 19 features: 40.43 s; time per pharmacophore: 0.40430 s
100 pharmacophores having 20 features: 58.30 s; time per pharmacophore: 0.58300 s
```

## Documentation
Mode documentation can be found here - https://pmapper.readthedocs.io/en/latest/

## Citation
Ligand-Based Pharmacophore Modeling Using Novel 3D Pharmacophore Signatures  
Alina Kutlushina, Aigul Khakimova, Timur Madzhidov, Pavel Polishchuk  
*Molecules* **2018**, 23(12), 3094  
https://doi.org/10.3390/molecules23123094

## License
BSD-3 clause
