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
|`type`|`dict`|{1: 1}|Types of each beads.<br>key: id of bead in chain (1, 2, ..., N)<br>value: type of bead<br>If the number of types in the input file is less than `N`, types are expanded repeatedly.<br>e.g. (N=5)<br>1: 1<br>2: 2<br>This is same as below.<br>1: 1<br>2: 2<br>3: 1<br>4: 2<br>5: 1|
|`btype`|`dict`|{1: 1}|Types of each bonds.<br>key: id of bond in chain (1, 2, ..., N-1)<br>value: type of bond<br>The expansion rule is same as `type`.|
|`atype`|`dict`|{1: 1}|Types of each angles.<br>key: id of angle in chain (1, 2, ..., N-2)<br>value: type of angle<br>The expansion rule is same as `type`.|
|`dtype`|`dict`|{1: 1}|Types of each dihedrals.<br>key: id of dihedral in chain (1, 2, ..., N-3)<br>value: type of dihedral<br>The expansion rule is same as `type`.|


## About program

### checking overlap
Every time one coordinate of a chain is generated, it is checked for overlap.
The rate of coordinate generation gradually slows down as the probability of overlap increases with the progress of coordinate generation.
[KDtree](https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KDTree.html) was used to check for duplicates.
Since all particles are checked, if `overlap_threashold` is larger than the bond length, coordinate generation will not proceed.

### decide dihedral angle
Completely random decision. Rotation of the axes was performed using the [Rodrigues' Rotation Formula](https://mathworld.wolfram.com/RodriguesRotationFormula.html).
