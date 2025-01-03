# Pmapper - 3D pharmacophore signatures and fingerprints

Pmapper is a Python module to generate 3D pharmacophore signatures and fingerprints.
Signatures uniquely encode 3D pharmacophores with hashes suitable for fast identification of identical pharmacophores.

## Dependency

`rdkit >= 2017.09`  
`networkx >= 2`

## Installation
```text
pip install pmapper
```

## Changelog
**1.0.0**
- added functionality to calculate 3D pharmacophore descriptors for molecules with exclusion of single atoms (for the purpose of model interpretation)
- added convenience function get_feature_ids
- added function add_feature to manually edit/construct a pharmacophore
- added save/load of pharmit pharmacophore models


- IMPORTANT: changed the hashing procedure to make it more stable (pickle dependency was removed). This breaks compatibility with previously generated md5 hashes with `get_signature_md5`, `iterate_pharm` and `iterate_pharm1` functions, all other functionality was not affected. 

**1.0.1**
- `fit_model` function can return rms by request

**1.0.2**
- `fit_model` function now returns a dict of mapped feature ids

**1.0.3**
- add `get_subpharmacophore` function
- fix `get_mirror_pharmacophore` function to use the same bin step and cached args as for the source pharmacophore instance  

**1.0.4**
- fix installation of dependency `networkx`
- add citations on examples of `pmapper` descriptors used for machine learning

**1.1**
- change SMARTS pattern to avoid matching positively charge nitrogen atoms as H-bond acceptors

**1.1.1**
- fix missing aromatic features when reading LigandScout models

**1.1.2**
- add treatment of exclusion volume features in LigandScout files

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
Pharmacophores can be created with enabled `cached` argument. This will speed up all futher repeated calls to retrive hash, fingerprints or descriptors.
```python
p = P(cached=True)
```

## Speed tests
Generation of pharmacophore signatures (hashes) is a CPU-bound task. The computation speed depends on the number of features in pharmacophores.  
Tests were run on a random subset of compounds from Drugbank. Up to 50 conformers were generated for each compound.   
Laptop configuration:
- Intel(R) Core(TM) i7-5500U CPU @ 2.40GHz
- 12 GB RAM
- calculation was run in 1 thread (the module is thread safe and calculations can be parallelized)

To run the speed test use `pmapper_speed_test` command line tool

```text
========== Reading of conformers of molecules ==========
329 molecules were read in 0.0134 s

========== Creation of pharmacophores (with enabled caching) ==========
1938 pharmacophores were created in 3.17065 s

========== First calculation of hashes ==========
2 pharmacophores with 0 features - 0.00014s or 7e-05s per pharmacophore
2 pharmacophores with 1 features - 0.0001s or 5e-05s per pharmacophore
12 pharmacophores with 2 features - 0.00042s or 3e-05s per pharmacophore
44 pharmacophores with 3 features - 0.00212s or 5e-05s per pharmacophore
100 pharmacophores with 4 features - 0.00933s or 9e-05s per pharmacophore
103 pharmacophores with 5 features - 0.05155s or 0.0005s per pharmacophore
105 pharmacophores with 6 features - 0.10857s or 0.00103s per pharmacophore
109 pharmacophores with 7 features - 0.25322s or 0.00232s per pharmacophore
117 pharmacophores with 8 features - 0.59508s or 0.00509s per pharmacophore
101 pharmacophores with 9 features - 0.8795s or 0.00871s per pharmacophore
105 pharmacophores with 10 features - 1.61349s or 0.01537s per pharmacophore
100 pharmacophores with 11 features - 2.24937s or 0.02249s per pharmacophore
103 pharmacophores with 12 features - 3.53308s or 0.0343s per pharmacophore
117 pharmacophores with 13 features - 6.49837s or 0.05554s per pharmacophore
103 pharmacophores with 14 features - 7.54796s or 0.07328s per pharmacophore
142 pharmacophores with 15 features - 14.92654s or 0.10512s per pharmacophore
104 pharmacophores with 16 features - 13.86378s or 0.13331s per pharmacophore
100 pharmacophores with 17 features - 17.94023s or 0.1794s per pharmacophore
120 pharmacophores with 18 features - 28.01455s or 0.23345s per pharmacophore
136 pharmacophores with 19 features - 42.53481s or 0.31276s per pharmacophore
113 pharmacophores with 20 features - 45.88228s or 0.40604s per pharmacophore

========== Second calculation of hashes of the same pharmacophores ==========
2 pharmacophores with 0 features - 5e-05s or 2e-05s per pharmacophore
2 pharmacophores with 1 features - 3e-05s or 1e-05s per pharmacophore
12 pharmacophores with 2 features - 0.00012s or 1e-05s per pharmacophore
44 pharmacophores with 3 features - 0.00041s or 1e-05s per pharmacophore
100 pharmacophores with 4 features - 0.00089s or 1e-05s per pharmacophore
103 pharmacophores with 5 features - 0.00166s or 2e-05s per pharmacophore
105 pharmacophores with 6 features - 0.00316s or 3e-05s per pharmacophore
109 pharmacophores with 7 features - 0.00707s or 6e-05s per pharmacophore
117 pharmacophores with 8 features - 0.0166s or 0.00014s per pharmacophore
101 pharmacophores with 9 features - 0.02005s or 0.0002s per pharmacophore
105 pharmacophores with 10 features - 0.03527s or 0.00034s per pharmacophore
100 pharmacophores with 11 features - 0.05271s or 0.00053s per pharmacophore
103 pharmacophores with 12 features - 0.08097s or 0.00079s per pharmacophore
117 pharmacophores with 13 features - 0.13274s or 0.00113s per pharmacophore
103 pharmacophores with 14 features - 0.1588s or 0.00154s per pharmacophore
142 pharmacophores with 15 features - 0.32687s or 0.0023s per pharmacophore
104 pharmacophores with 16 features - 0.29255s or 0.00281s per pharmacophore
100 pharmacophores with 17 features - 0.38286s or 0.00383s per pharmacophore
120 pharmacophores with 18 features - 0.61327s or 0.00511s per pharmacophore
136 pharmacophores with 19 features - 0.93486s or 0.00687s per pharmacophore
113 pharmacophores with 20 features - 0.94041s or 0.00832s per pharmacophore
```

## Documentation
More documentation can be found here - https://pmapper.readthedocs.io/en/latest/

## Citation
Ligand-Based Pharmacophore Modeling Using Novel 3D Pharmacophore Signatures  
Alina Kutlushina, Aigul Khakimova, Timur Madzhidov, Pavel Polishchuk  
*Molecules* **2018**, 23(12), 3094  
https://doi.org/10.3390/molecules23123094

##### Further publications

###### MD pharmacophores
Virtual Screening Using Pharmacophore Models Retrieved from Molecular Dynamic Simulations  
Pavel Polishchuk, Alina Kutlushina, Dayana Bashirova, Olena Mokshyna, Timur Madzhidov  
*Int. J. Mol. Sci.* **2019**, 20(23), 5834  
https://doi.org/10.3390/ijms20235834

###### Pmapper descriptors in machine learning
QSAR Modeling Based on Conformation Ensembles Using a Multi-Instance Learning Approach  
Zankov, D. V.; Matveieva, M.; Nikonenko, A. V.; Nugmanov, R. I.; Baskin, I. I.; Varnek, A.; Polishchuk, P.; Madzhidov, T. I.  
*J. Chem. Inf. Model.* **2021**, 61 (10), 4913-4923
https://doi.org/10.1021/acs.jcim.1c00692  
  
Multi-Instance Learning Approach to the Modeling of Enantioselectivity of Conformationally Flexible Organic Catalysts  
Zankov, D.; Madzhidov, T.; Polishchuk, P.; Sidorov, P.; Varnek, A.  
*J. Chem. Inf. Model.* **2023**, 63 (21), 6629-6641  
https://doi.org/10.1021/acs.jcim.3c00393  


## License
BSD-3 clause
