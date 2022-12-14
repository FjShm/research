# Create Straight Chains


## Description

This program creates an initial LAMMPS data file with arbitrary length and arbitrary number of straight chains.
(straight chain = chain without side chain)
Length of each bond and each angle can be specified. Dihedral angles are not supported (random).
Duplicate particles are also checked so that no other particles exist within the specified range.


## How to use
```bash
python create_straight_chains.py [-in path/to/input-yaml-file (default=input.yaml)] [-out name_of_exported_LAMMPS_data_file (default=init.data)]
```

## input file
|key|format|default value|explanation|
|---|---|---|---|
|`mass`|`dict`|{1: 68.12}|mass of each types of beads.<br>key: bead type<br>value: mass|
|`bond_length`|`dict`|{1: 5.}|bond length of each types of bonds.<br>key: bond type<br>value: length|
|`angle_degree`|`dict`|{1: 120.}|angle of each types of angles.<br>key: angle type<br>value: angle|
|`cell_length`|`float`|120.|Length of one side of the simulation box.|
|`N`|`int`|24|# of particles in a single chain.|
|`M`|`int`|50|# of chains.|
|`origin`|`str`|""|Where in the simulation box is the origin O of the coordinates. If "center", the body center of the simulation box is the origin O. If not, The corner of the simulation box is the origin O.|
|`sigma_bond`|`float`|1.0|The bond length is calculated randomly according to a normal distribution with the value set by `bond_length` as μ and `sigma_bond` as σ.|
|`sigma_angle`|`float`|1.0|Same as `sigma_bond`|
|`overlap_threashold`|`float`|1e-10|Threshold at which particles are judged to overlap. If there are other particles within this distance of a particle when generating coordinates, the coordinate generation is redone.|
|`special_bonds`|`int`|0|Beads specified in special_bonds are not included in the overlap calculation.|
|`type`|`dict`|{1: 1}|Types of each beads.<br>key: id of bead in chain (1, 2, ..., N)<br>value: type of bead<br>If the number of types in the input file is less than `N`, types are expanded repeatedly.<br>e.g. (N=5)<br>1: 1<br>2: 2<br>This is same as below.<br>1: 1<br>2: 2<br>3: 1<br>4: 2<br>5: 1|
|`btype`|`dict`|{1: 1}|Types of each bonds.<br>key: id of bond in chain (1, 2, ..., N-1)<br>value: type of bond<br>The expansion rule is same as `type`.|
|`atype`|`dict`|None|Types of each angles.<br>key: id of angle in chain (1, 2, ..., N-2)<br>value: type of angle<br>The expansion rule is same as `type`.<br>If you don't set this parameter, "Angles" is not written in LAMMPS data file.|
|`dtype`|`dict`|None|Types of each dihedrals.<br>key: id of dihedral in chain (1, 2, ..., N-3)<br>value: type of dihedral<br>The expansion rule is same as `type`.<br>If you don't set this parameter, "Dihedrals" is not written in LAMMPS data file.|


## About program

### checking overlap
Every time one coordinate of a chain is generated, it is checked for overlap.
The rate of coordinate generation gradually slows down as the probability of overlap increases with the progress of coordinate generation.
[KDtree](https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KDTree.html) is used to check for duplicates.
If you want to set `overlap_threashold` larger than the bond length, set also `special_bonds` larger than 1.

### decide dihedral angle
Completely random decision. Rotation of the axes was performed using the [Rodrigues' Rotation Formula](https://mathworld.wolfram.com/RodriguesRotationFormula.html).


## Sample System

### Kremer-Grest model

- $N=100, M=100$
- 1 bond type
- No angles
- No dihedrals

```yaml
mass:
  1: 1.
bond_length:
  1: 1.
angle_degree: ~
cell_length: 20.
N: 100
M: 100
origin: "center"  # "center" or not
sigma_bond: 0.01
  #sigma_angle: 0.1
overlap_threashold: 0.1
special_bonds: 1

# < Structures > #

type:
  1: 1

btype:
  1: 1

#atype:
#  1: 1

#dtype:
#  1: 1

```
![KG](https://user-images.githubusercontent.com/100293098/207528049-d5ba5c23-16e8-4290-be07-a9e91c61a982.png)

WallTime: 0m 29s

### 1 bead type

- $N=24, M=512$
- 1 bond type
- 1 angle type
- 1 dihedral type

```yaml
special_bonds: 3
mass:
  1: 68.117
bond_length:
  1: 5.
angle_degree:
  1: 120.
cell_length: 150.
N: 24
M: 512
origin: "center"  # "center" or not
sigma_bond: 0.01
sigma_angle: 0.1
overlap_threashold: 2.

# < Structures > #

type:
  1: 1

btype:
  1: 1

atype:
  1: 1

dtype:
  1: 1

```

![1_bead_type](https://user-images.githubusercontent.com/100293098/207528380-789ff2ce-dc55-4910-a9d7-594d93e8e3bf.png)

WallTime: 1m 48s


### 2 bead types

- $N=49, M=512$
- 2 bond types
- 2 angle types
- 2 dihedral types

```yaml
mass:
  1: 28.05316
  2: 40.06386
bond_length:
  1: 2.5
  2: 3.1
angle_degree:
  1: 90.
  2: 180.
cell_length: 150.
N: 49
M: 512
origin: "center"  # "center" or not
sigma_bond: 0.01
sigma_angle: 0.1
overlap_threashold: 1.5
special_bonds: 3

# < Structures > #

type:
  1: 1
  2: 2

btype:
  1: 1
  2: 2

atype:
  1: 1
  2: 2

dtype:
  1: 1
  2: 2

```

![2_bead_types](https://user-images.githubusercontent.com/100293098/207529716-181e82a7-8c12-40b7-817a-ea09ead4fb33.png)

WallTime: 17m 51s



### 3 bead types

- $N=49, M=512$
- 4 bond types
- 4 angle types
- 4 dihedral types

```yaml
mass:
  1: 28.05316
  2: 40.06386
  3: 31.033925
bond_length:
  1: 2.5
  2: 3.1
  3: 2.5
  4: 3.1
angle_degree:
  1: 90.
  2: 180.
  3: 90.
  4: 90.
cell_length: 150.
N: 49
M: 512
origin: "center"  # "center" or not
sigma_bond: 0.01
sigma_angle: 0.1
overlap_threashold: 1.5
special_bonds: 3

# < Structures > #

type:
  1: 3
  2: 2
  3: 1
  4: 2
  5: 1
  6: 2
  7: 1
  8: 2
  9: 1
  10: 2
  11: 1
  12: 2
  13: 1
  14: 2
  15: 1
  16: 2
  17: 1
  18: 2
  19: 1
  20: 2
  21: 1
  22: 2
  23: 1
  24: 2
  25: 1
  26: 2
  27: 1
  28: 2
  29: 1
  30: 2
  31: 1
  32: 2
  33: 1
  34: 2
  35: 1
  36: 2
  37: 1
  38: 2
  39: 1
  40: 2
  41: 1
  42: 2
  43: 1
  44: 2
  45: 1
  46: 2
  47: 1
  48: 2
  49: 3

btype:
  1: 3
  2: 2
  3: 1
  4: 2
  5: 1
  6: 2
  7: 1
  8: 2
  9: 1
  10: 2
  11: 1
  12: 2
  13: 1
  14: 2
  15: 1
  16: 2
  17: 1
  18: 2
  19: 1
  20: 2
  21: 1
  22: 2
  23: 1
  24: 2
  25: 1
  26: 2
  27: 1
  28: 2
  29: 1
  30: 2
  31: 1
  32: 2
  33: 1
  34: 2
  35: 1
  36: 2
  37: 1
  38: 2
  39: 1
  40: 2
  41: 1
  42: 2
  43: 1
  44: 2
  45: 1
  46: 2
  47: 1
  48: 4

atype:
  1: 3
  2: 2
  3: 1
  4: 2
  5: 1
  6: 2
  7: 1
  8: 2
  9: 1
  10: 2
  11: 1
  12: 2
  13: 1
  14: 2
  15: 1
  16: 2
  17: 1
  18: 2
  19: 1
  20: 2
  21: 1
  22: 2
  23: 1
  24: 2
  25: 1
  26: 2
  27: 1
  28: 2
  29: 1
  30: 2
  31: 1
  32: 2
  33: 1
  34: 2
  35: 1
  36: 2
  37: 1
  38: 2
  39: 1
  40: 2
  41: 1
  42: 2
  43: 1
  44: 2
  45: 1
  46: 2
  47: 4

dtype:
  1: 3
  2: 2
  3: 1
  4: 2
  5: 1
  6: 2
  7: 1
  8: 2
  9: 1
  10: 2
  11: 1
  12: 2
  13: 1
  14: 2
  15: 1
  16: 2
  17: 1
  18: 2
  19: 1
  20: 2
  21: 1
  22: 2
  23: 1
  24: 2
  25: 1
  26: 2
  27: 1
  28: 2
  29: 1
  30: 2
  31: 1
  32: 2
  33: 1
  34: 2
  35: 1
  36: 2
  37: 1
  38: 2
  39: 1
  40: 2
  41: 1
  42: 2
  43: 1
  44: 2
  45: 1
  46: 4
```
![3types](https://user-images.githubusercontent.com/100293098/207531446-908429a9-f929-4280-aa13-2e0e68f3a906.png)

WallTime: 18m 3s

