# Pmapper - 3D pharmacophore signatures and fingerprints

Pmapper is a Python module to generate 3D pharmacophore signatures and fingerprints.
Signatures uniquely encode 3D pharmacophores with hashes suitable for fast identification of identical pharmacophores.

## Dependency

`rdkit >= 2017.09`
`networkx >= 1.11`

## Examples

### Load modules
```python
from rdkit import Chem
from rdkit.Chem import AllChem, ChemicalFeatures
from pharmacophore import Pharmacophore as P, read_smarts_feature_file, load_multi_conf_mol
from pprint import pprint
```
### Create pharmacophore from a single comformer using feature definition from SMARTS file
```
# load a molecule from SMILES and generate 3D coordinates
mol = Chem.MolFromSmiles('C1CC(=O)NC(=O)C1N2C(=O)C3=CC=CC=C3C2=O')  # talidomide
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)

# load pharmacophore feature definitions from SMARTS file
smarts = read_smarts_feature_file('smarts_features.txt')

# create pharmacophore
p = P()
p.load_from_smarts(mol, smarts)
```
### Get 3D pharmacophore signature
```python
# get 3D pharmacophore signature
sig = p.get_signature_md5()
print(sig)
```
Output:
```
f2e16f52f6f6ca6e97fc5844bfd35d36
```
### Get 3D pharmacophore signature with non-zero tolerance
```python
sig = p.get_signature_md5(tol=5)
print(sig)
```
Output:
```
fb535302db2e5d624aa979b6e8dfbdf2
```
### Create pharmacophore from a single comformer using RDKit feature factory
```python
# load pharmacophore using RDKit factory and get 3D pharmacophore signature
factory = ChemicalFeatures.BuildFeatureFactory('smarts_features.fdef')

p.load_from_feature_factory(mol, factory)
sig = p.get_signature_md5()
print(sig)
```
Output:
```
f2e16f52f6f6ca6e97fc5844bfd35d36
```
### Create pharmacophores for a multiple conformer compound
```python
# create multiple conformer molecule
AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=1024)

ps = load_multi_conf_mol(mol, smarts_features=smarts)

sig = [p.get_signature_md5() for p in ps]

pprint(sorted(sig))  # identical signatures occur
```
Output:
```
['13d168458ab1f251157f2422efcce312',
 '13d168458ab1f251157f2422efcce312',
 '182a4cfa756fe8b7f736a7f7ac0e8e0a',
 '182a4cfa756fe8b7f736a7f7ac0e8e0a',
 '4234e9d249874a5009f1e312dd885d80',
 'ab273dd083c4f2e3424ba917b121b846',
 'b6ec58553d2984bd398b4520bd1545cc',
 'bfc43365b2657d08b6bb888e4d8ec71b',
 'f5ca8e406dae31182e2b06fde7452b75',
 'fc4a85e818fc0b3f034a7af42fa5ca69']
 ```
### Generate 3D pharmacophore fingerprint
```python
# generate 3D pharmacophore fingerprint which takes into account stereoconfiguration
b = p.get_fp(min_features=4, max_features=4)   # set of activated bits
print(b)
```
Output (a set of activated bit numbers):
```
{1922, 1795, 779, 1040, 528, 920, 154, 1437, 287, 1313, 1447, 1961, 941, 690, 1203, 65, 1346, 709, 1486, 1366, 2006, 1750, 1016, 346, 603, 1116, 354, 995, 228, 2024, 1900, 1524, 888, 2043}
```
Change settings:
```python
b = p.get_fp(min_features=4, max_features=4, nbits=1024, activate_bits=2)
print(b)
```
Output (a set of activated bit numbers):
```
{897, 514, 259, 389, 520, 264, 143, 16, 529, 656, 787, 660, 24, 285, 157, 32, 673, 550, 683, 173, 301, 558, 45, 945, 177, 692, 950, 443, 444, 61, 960, 961, 448, 321, 709, 197, 587, 460, 77, 718, 720, 80, 339, 596, 723, 470, 980, 345, 601, 476, 354, 614, 743, 1003, 875, 494, 367, 497, 114, 1012, 244, 630, 377, 762, 507, 508, 1021}
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
### Support LigandScout pml-files
LigandScout models saved as pml-files can be read using `p.load_ls_model`. Also a pharmacophore can be stored in this format in order to export to LigandScout (`p.save_ls_model`).

## Citation
...

## License
BSD-3 clause