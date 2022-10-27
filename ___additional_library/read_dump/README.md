# read\_dump.h
2022.10.27 Fujii Shoma

---

## Desctiption
Provides a class that reads `*.lammpstrj` files and holds their various data.
Available by simply placing a header(`read_dump.h`) in the include directory.

### member variables
| variable name |type| data |
|---|---|---|
|`timestep`|int|timestep of the frame|
|`num_atoms`|int|# of atoms in the frame|
|`num_frames`|int|# of total frames|
|`cellbox_a`|Eigen::Vector3d|basis vector of the simulation box $\boldsymbol a$ (see [LAMMPS dump box bounds](https://docs.lammps.org/dump.html), [LAMMPS triclinic](https://docs.lammps.org/Howto_triclinic.html))|
|`cellbox_b`|Eigen::Vector3d|basis vector of the simulation box $\boldsymbol b$|
|`cellbox_c`|Eigen::Vector3d|basis vector of the simulation box $\boldsymbol c$|
|`cellbox_origin`|Eigen::Vector3d|Origin of the simulation box|
|`atoms_all_data`|Eigen::MatrixXd\*|All data listed in the ATOMS section. If "id" exists, the row index is sorted in the order of "id". If "id" does not exist, the order is the same as in the dump file.|
|`header_map`|std::map\<std::string, int\>\*|A map that maps header names (such as "id", "xu", ...) to their column numbers.|

### member functions
|function name|type|details|
|---|---|---|
|fileopen(const std::string&)|void|Open the lammpstrj file with the path given as an argument.|
|clear(void)|void|-|
|join_3columns(<br>std::vector\<Eigen::Vector3d\>&,<br>std::string,<br>std::string,<br>std::string<br>)|void|The data of the columns matching the three specified header names are concatenated and returned.|
|read\_1frame(void)|void|Read through the dump file for one frame. If read\_all\_frames, described below, has already been executed, it only refers to the next frame's data from memory.|
|read\_all\_frames(void)|void|Load all data in the dump file.|
|header\_validation(...std::string)|void|Verify that the required column names are present in the header. Multiple column names (std::string) may be packed into a vector type and passed as an argument, or any number of column names (std::string) may be passed as an argument.|
|set\_wanted\_frames(const std::vector\<double\>&)<br>set\_wanted\_frames(const std::vector\<int\>&)|void|Set the timestep of the frame to be read or the position of the frame in all frames (0-1).|
|search\_nearest\_timestep(void)|int|The position of the frame within all frames is given as a number between 0 and 1, and the timestep closest to that value is returned.|
|check\_if\_wanted\_frame(void)|bool|Check to see if the frame being read is one of the frames to be read. If set\_wanted\_frames has been executed in advance, it returns true/false according to the set contents. If not, only true is returned.|

## How to use
### Install
Copy `read_dump.h` to your include directory, or same directory as code to compile.

### Construct
```
ReadDump::ReadDump rd;
rd.open("path/to/.lammpstrj");

// or
ReadDump::ReadDump rd("path/to/.lammpstrj");
```

### sample 1
Read the dump file one frame at a time and display the coordinates of atom id=100.
```c++
ReadDump::ReadDump rd("dump.u.lammpstrj");

while(rd.read_1frame()){
  std::cout << rd.timestep << std::endl;
  rd.header_validation("id", "mol", "xu", "yu", "zu");
  
  // display coordination of id=100 (2 ways) //
  int id = 100;
  
  // method 1: direct reference
  std::cout << "method 1\n";
  int xu = rd.header_map->at("xu");
  int yu = rd.header_map->at("yu");
  int zu = rd.header_map->at("zu");
  std::cout << rd.atoms_all_data->coeff(id-1, xu) << std::endl;
  std::cout << rd.atoms_all_data->coeff(id-1, yu) << std::endl;
  std::cout << rd.atoms_all_data->coeff(id-1, zu) << std::endl;
  
  // method 2: use join_3columns()
  std::cout << "method 2\n";
  std::vector<Eigen::Vector3d> coordinate;
  rd.join_3columns(coordinate, "xu", "yu", "zu");
  std::cout << coordinate[id-1] << std::endl;
}
```

### sample 2
Load all the data in the dump file first and display the coordinates of atom id=100.
```c++
ReadDump::ReadDump rd("dump.u.lammpstrj");
rd.read_all_frames();

while(rd.read_1frame()){
  std::cout << rd.timestep << std::endl;
  rd.header_validation("id", "mol", "xu", "yu", "zu");
  
  // display coordination of id=100 (2 ways) //
  int id = 100;
  
  // method 1: direct reference
  std::cout << "method 1\n";
  int xu = rd.header_map->at("xu");
  int yu = rd.header_map->at("yu");
  int zu = rd.header_map->at("zu");
  std::cout << rd.atoms_all_data->coeff(id-1, xu) << std::endl;
  std::cout << rd.atoms_all_data->coeff(id-1, yu) << std::endl;
  std::cout << rd.atoms_all_data->coeff(id-1, zu) << std::endl;
  
  // method 2: use join_3columns()
  std::cout << "method 2\n";
  std::vector<Eigen::Vector3d> coordinate;
  rd.join_3columns(coordinate, "xu", "yu", "zu");
  std::cout << coordinate[id-1] << std::endl;
}
```

### sample 3-1
Pre-set the frames you want to load, **by specifying ratio**.
```c++
ReadDump::ReadDump rd("dump.u.lammpstrj");
rd.read_all_frames(); // MUST
std::vector<double> ratio = {0., 0.4, 1.0}; // In no particular order, duplicates are acceptable. Must be 0 <= ratio <= 1.
rd.set_want_frames(ratio);

while(rd.read_1frame()){
  if (!rd.check_if_wanted_frame()){
    std::cout << "\n!! skip to read timestep " << rd.timestep << std::endl;
    continue;
  }
  std::cout << rd.timestep << std::endl;
  
  // display coordination of id=100
  int id = 100;
  std::vector<Eigen::Vector3d> coordinate;
  rd.join_3columns(coordinate, "xu", "yu", "zu");
  std::cout << coordinate[id-1] << std::endl;
}
```

### sample 3-2
Pre-set the frames you want to load, **by specifying timesteps**.
```c++
ReadDump::ReadDump rd("dump.u.lammpstrj");
rd.read_all_frames(); // MUST
std::vector<int> timestep = {0, 20000000, 40000000}; // In no particular order, duplicates are acceptable. Timesteps that do not exist are ignored.
rd.set_want_frames(timestep);

while(rd.read_1frame()){
  if (!rd.check_if_wanted_frame()){
    std::cout << "\n!! skip to read timestep " << rd.timestep << std::endl;
    continue;
  }
  std::cout << rd.timestep << std::endl;
  
  // display coordination of id=100
  int id = 100;
  std::vector<Eigen::Vector3d> coordinate;
  rd.join_3columns(coordinate, "xu", "yu", "zu");
  std::cout << coordinate[id-1] << std::endl;
}
```

## Restrictions
- `Eigen` library is needed
- g++17 or later is needed to compile `sample/test.cpp`
