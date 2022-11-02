# write\_dump.h
2022.10.27 Fujii Shoma

---

## Desctiption
Provides a class that writes `*.lammpstrj` files.
Available by simply placing a header(`write_dump.h`) in the include directory.

### member variables
| variable name |type| data |
|---|---|---|
|`timestep`|int\*|timestep of the frame|
|`num_atoms`|int\*|# of atoms in the frame|
|`num_frames`|int\*|# of total frames|
|`cellbox_a`|Eigen::Vector3d\*|basis vector of the simulation box $\boldsymbol a$ (see [LAMMPS dump box bounds](https://docs.lammps.org/dump.html), [LAMMPS triclinic](https://docs.lammps.org/Howto_triclinic.html))|
|`cellbox_b`|Eigen::Vector3d\*|basis vector of the simulation box $\boldsymbol b$|
|`cellbox_c`|Eigen::Vector3d\*|basis vector of the simulation box $\boldsymbol c$|
|`cellbox_origin`|Eigen::Vector3d\*|Origin of the simulation box|
|`atoms_all_data`|Eigen::MatrixXd\*|All data listed in the ATOMS section. If "id" exists, the row index is sorted in the order of "id". If "id" does not exist, the order is the same as in the dump file.|
|`header_map`|std::map\<std::string, int\>\*|A map that maps header names (such as "id", "xu", ...) to their column numbers.|

### member functions
|function name|type|details|
|---|---|---|
|open(const std::string&)|void|Open the lammpstrj file with the path given as an argument.|
|clear(void)|void|-|
|write\_1frame(void)|void|Read through the dump file for one frame. If read\_all\_frames, described below, has already been executed, it only refers to the next frame's data from memory.|
|set\_by\_ReadDump(ReadDump::ReadDump&)<br>set\_by\_ReadDump(ReadDump::ReadDump&)|void|Copy all member variable's pointers form readdump to writedump.|
|set\_headers(const std::vector\<std::string\>&)<br>set\_headers(std::string...)|void|Set header names to write to dump.|

## How to use
### Install
Copy `write_dump.h` to your include directory, or same directory as code to compile.

### Construct
```
WriteDump::WriteDump wd;
wd.open("path/to/.lammpstrj");

// or
WriteDump::WriteDump wd("path/to/.lammpstrj");
```

### sample 1




## Restrictions
- `Eigen` library is needed
- `read_dump.h` is needed
- g++17 or later is needed
