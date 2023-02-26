import os
import sys
import argparse
import yaml
import random
import numpy as np
import datetime
import pprint
import logging
from sklearn.neighbors import KDTree
from tqdm import tqdm

"""
# How to use

```bash
# if input filename is "input.yaml"
# and output LAMMPS data filename is "init.data"
python create_A-model.py

# See help
python create_straight_chains.py -h
```
"""


def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s:%(name)s] %(message)s"
    )
    arg = Arg()
    if not os.path.exists(arg.input):
        exit(f"Error: Input file does not exist.\n- {arg.input}")
    logging.info(f"Start to read the following file: {arg.input}")

    f = open(arg.input)
    inputs = yaml.safe_load(f)
    inputs = {} if inputs is None else inputs
    f.close()

    params = PARAMs(inputs.copy())
    params.validate()
    params.show()
    del inputs

    logging.info("Now creating position (this may take some time).")
    atoms = ATOMS(params)
    atoms.main()

    logging.info("Complete to create positions.")
    writer = Writer()
    logging.info("Preparing Info...")
    writer.abouts(params)
    logging.info("Preparing Cells...")
    writer.cells(params)
    logging.info("Preparing Masses...")
    writer.masses(params)
    logging.info("Preparing Atoms...")
    writer.atoms(params, atoms)
    logging.info("Preparing Bonds...")
    writer.bonds(params)
    logging.info("Preparing Angles...")
    writer.angles(params)
    logging.info("Preparing Dihedrals...")
    writer.dihedrals(params)
    logging.info("Now writing...")
    writer.write(arg.output)
    logging.info(f"All done!\n\n- exported: {arg.output}\n")


class DefaultParams:
    """
    mass [g/cm3]
    cite: https://www.chemistry.or.jp/know/atom_2022.pdf
    C = average(12.0096, 12.0116)
    H = average(1.00784, 1.00811)
    polyisoprene unit = C5H8
    """

    __mass: dict = {1: 68.12}  # cis-isoprene [g/cm3]
    __l: dict = {1: 5.0}
    __theta: dict = {1: 2.0 / 3.0 * np.pi}  # [rad]
    __cell_length: float = 120.0  # [Ang]
    __N: int = 24
    __M: int = 50
    __origin: str = ""
    # "center" or not.  Where does the
    # simulation box coincide with the origin
    # of the coordinates.

    __sigma_bond: float = 1.0
    __sigma_angle: float = 1.0
    # Bond lengths are generated according to
    # a normal distribution with bond_length_Ang
    # as the mean and sigma_bond as the standard
    # deviation.
    # sigma_angle is same as sigma_bond.

    __overlap_threa: float = 1e-10
    # Threshold at which particles are judged
    # to overlap each other.

    __bead_type: dict = {1: 1}
    __btype: dict = {1: 1}
    __atype: dict = None
    __dtype: dict = None

    __isAngle: bool = False
    __isDihedral: bool = False

    __special_bonds: int = 3

    @property
    def mass(self) -> dict:
        return self.__mass.copy()

    @property
    def l(self) -> dict:
        return self.__l.copy()

    @property
    def theta(self) -> dict:
        return self.__theta.copy()

    @property
    def cell_length(self) -> float:
        return self.__cell_length

    @property
    def N(self) -> int:
        return self.__N

    @property
    def M(self) -> int:
        return self.__M

    @property
    def origin(self) -> str:
        return self.__origin

    @property
    def sigma_bond(self) -> float:
        return self.__sigma_bond

    @property
    def sigma_angle(self) -> float:
        return self.__sigma_angle

    @property
    def overlap_threa(self) -> float:
        return self.__overlap_threa

    @property
    def bead_type(self) -> dict:
        return self.__bead_type.copy()

    @property
    def btype(self) -> dict:
        return self.__btype.copy()

    @property
    def atype(self) -> dict:
        return None if self.__atype is None else self.__atype.copy()

    @property
    def dtype(self) -> dict:
        return None if self.__dtype is None else self.__dtype.copy()

    @property
    def special_bonds(self) -> int:
        return self.__special_bonds

    @property
    def isAngle(self) -> bool:
        return self.__isAngle

    @property
    def isDihedral(self) -> bool:
        return self.__isDihedral

    @mass.setter
    def _mass(self, val: dict) -> None:
        if val is not None:
            self.__mass = val.copy() if type(val) == dict else {1: val}

    @l.setter
    def _l(self, val: dict) -> None:
        if val is not None:
            self.__l = val.copy() if type(val) == dict else {1: val}

    @theta.setter
    def _thetafdeg(self, val: dict) -> None:
        if val is not None:
            if type(val) == dict:
                for k in val.keys():
                    val[k] = np.deg2rad(val[k])
                self.__theta = val.copy()
            else:
                self.__theta = {1: val}

    @cell_length.setter
    def _cell_length(self, val: float) -> None:
        if val is not None:
            self.__cell_length = val

    @N.setter
    def _N(self, val: int) -> None:
        if val is not None:
            self.__N = val

    @M.setter
    def _M(self, val: int) -> None:
        if val is not None:
            self.__M = val

    @origin.setter
    def _origin(self, val: str) -> None:
        if val is not None:
            self.__origin = val

    @sigma_bond.setter
    def _sigma_bond(self, val: float) -> None:
        if val is not None:
            self.__sigma_bond = val

    @sigma_angle.setter
    def _sigma_angle(self, val: float) -> None:
        if val is not None:
            self.__sigma_angle = val

    @overlap_threa.setter
    def _overlap_threa(self, val: float) -> None:
        if val is not None:
            self.__overlap_threa = val

    @bead_type.setter
    def _bead_type(self, val: dict) -> None:
        if val is not None:
            if type(val) != dict:
                exit("The value of 'type' must be dict.")
            self.__bead_type = dict((k, v) for k, v in sorted(val.items()))

    @btype.setter
    def _btype(self, val: dict) -> None:
        if val is not None:
            if type(val) != dict:
                exit("The value of 'btype' must be dict.")
            self.__btype = dict((k, v) for k, v in sorted(val.items()))

    @atype.setter
    def _atype(self, val: dict) -> None:
        if val is not None:
            if type(val) != dict:
                exit("The value of 'atype' must be dict.")
            self.__isAngle = True
            self.__atype = dict((k, v) for k, v in sorted(val.items()))

    @dtype.setter
    def _dtype(self, val: dict) -> None:
        if val is not None:
            if type(val) != dict:
                exit("The value of 'dtype' must be dict.")
            self.__isDihedral = True
            self.__dtype = dict((k, v) for k, v in sorted(val.items()))

    @special_bonds.setter
    def _special_bonds(self, val: int) -> None:
        if val is not None:
            if val < 0:
                exit(f"special_bonds must be over 0.\ninvalid value: {val}")
            self.__special_bonds = val


class PARAMs(DefaultParams):
    def __init__(self, inputs: dict) -> None:
        for key, val in inputs.items():
            if key == "mass":
                self._mass = val
            elif key == "bond_length":
                self._l = val
            elif key == "angle_degree":
                self._thetafdeg = val
            elif key == "cell_length":
                self._cell_length = val
            elif key == "N":
                self._N = val
            elif key == "M":
                self._M = val
            elif key == "origin":
                self._origin = val
            elif key == "sigma_bond":
                self._sigma_bond = val
            elif key == "sigma_angle":
                self._sigma_angle = val
            elif key == "overlap_threashold":
                self._overlap_threa = val
            elif key == "type":
                self._bead_type = val
            elif key == "btype":
                self._btype = val
            elif key == "atype":
                self._atype = val
            elif key == "dtype":
                self._dtype = val
            elif key == "special_bonds":
                self._special_bonds = val
        self._expand_types()

    def __expand_types(self, xtype: dict, Nmax: int) -> dict:
        if xtype is None:
            return xtype
        nxtype = len(xtype)
        if nxtype < Nmax:
            for i in range(nxtype + 1, Nmax + 1):
                refid = i % nxtype
                refid = nxtype if refid == 0 else refid
                xtype[i] = xtype[refid]
        return xtype

    def _expand_types(self) -> None:
        self._bead_type = self.__expand_types(self.bead_type, self.N)
        self._btype = self.__expand_types(self.btype, self.N - 1)
        self._atype = self.__expand_types(self.atype, self.N - 2)
        self._dtype = self.__expand_types(self.dtype, self.N - 3)

    def validate(self) -> None:
        # bead type
        types_of_mass = sorted(set(self.mass.keys()))
        types = sorted(set(self.bead_type.values()))
        if types_of_mass != types:
            logging.error("invalid bead type\n")
            exit(f"types of mass: {types_of_mass}\ntype: {types}")

    def show(self) -> None:
        pprint.pprint(vars(self))


class ATOMS:
    atom_coordinations: np.ndarray
    mirror_indexes: np.array

    def __init__(self, params: PARAMs) -> None:
        self.__N = params.N
        self.__M = params.M
        self.__theta = params.theta
        self.__l = params.l
        self.__cell_length = params.cell_length
        self.__sigma_bond = params.sigma_bond
        self.__sigma_angle = params.sigma_angle
        self.__overlap_threa = params.overlap_threa
        self.__origin = params.origin
        self.__btype = params.btype
        self.__atype = params.atype
        self.__special_bonds = params.special_bonds
        self.atom_coordinations = np.full(
            (self.__N * self.__M, 3), (self.__cell_length + self.__overlap_threa) + 1.0
        )
        self.mirror_indexes = np.empty((self.__N * self.__M, 3), dtype=int)

    def main(self) -> None:
        for m in tqdm(range(self.__M)):
            while True:
                chain = self.__generate_1_chain()
                chain, mirror = self.__fix_periodic(chain)
                X = self.atom_coordinations[: m * self.__N, :].copy()
                if not self.__overlapped(chain, X):
                    self.atom_coordinations[
                        m * self.__N : (m + 1) * self.__N, :
                    ] = chain
                    self.mirror_indexes[m * self.__N : (m + 1) * self.__N, :] = mirror
                    break
        if self.__origin == "center":
            self.atom_coordinations -= 0.5 * self.__cell_length

    def __generate_next_bond_vector(
        self, vec_before: np.array, bead_id: int
    ) -> np.array:
        next_vec = self.__rect2spherical(vec_before)
        # decide bond length
        next_vec[0] = np.random.normal(
            loc=self.__l[self.__btype[bead_id - 1]], scale=self.__sigma_bond, size=1
        )[0]
        # decide angle
        dtheta = (
            random.randrange(-1, 2, 2)
            * np.random.normal(
                loc=np.pi - self.__theta[self.__atype[bead_id - 2]],
                scale=self.__sigma_angle,
                size=1,
            )[0]
        ) if self.__atype is not None else random.uniform(-np.pi, np.pi)
        next_vec[1] += dtheta
        # rotate according to Rodrigues' rotation formula
        next_vec = self.__spherical2rect(next_vec)
        phi = random.uniform(-np.pi, np.pi)
        n = vec_before / np.linalg.norm(vec_before)
        next_vec = self.__rotation(phi, n, next_vec)
        return next_vec

    def __generate_1st_bond_vector(self) -> np.array:
        r = np.random.normal(
            loc=self.__l[self.__btype[1]], scale=self.__sigma_bond, size=1
        )[0]
        t = random.uniform(0.0, np.pi)
        p = random.uniform(-np.pi, np.pi)
        return self.__spherical2rect(np.array([r, t, p]))

    def __spherical2rect(self, vec: np.array) -> np.array:
        x = vec[0] * np.sin(vec[1]) * np.cos(vec[2])
        y = vec[0] * np.sin(vec[1]) * np.sin(vec[2])
        z = vec[0] * np.cos(vec[1])
        return np.array([x, y, z])

    def __rect2spherical(self, vec: np.array) -> np.array:
        r = np.linalg.norm(vec)
        t = np.arccos(vec[2] / r)
        p = np.arccos(vec[0] / (r * np.sin(t)))
        return np.array([r, t, p])

    def __generate_1_chain(self) -> np.ndarray:
        r = np.empty((self.__N, 3))
        r[0, :] = np.random.rand(3) * self.__cell_length
        bond_vec = self.__generate_1st_bond_vector()
        r[1, :] = r[0, :] + bond_vec
        for n in range(2, self.__N):
            bond_vec = self.__generate_next_bond_vector(bond_vec, n + 1)
            r[n, :] = r[n - 1, :] + bond_vec
        return r

    def __append_periodic_skin(self, arr: np.ndarray) -> np.ndarray:
        if len(arr) == 0:
            return arr
        for i in range(3):
            rng = (
                arr[:, i] <= self.__overlap_threa,
                self.__cell_length - self.__overlap_threa <= arr[:, i],
            )
            skin = [arr[rng[0]].copy(), arr[rng[1]].copy()]
            skin[0][:, i] += (self.__cell_length,)
            skin[1][:, i] -= (self.__cell_length,)
            arr = np.array(list(arr) + list(skin[0]) + list(skin[1]))
        return arr

    def __overlapped(self, query: np.ndarray, X: np.ndarray) -> bool:
        X = self.__append_periodic_skin(X.copy())
        query = self.__append_periodic_skin(query.copy())
        tree = KDTree(np.array(list(X) + list(query)))
        ind = tree.query_radius(query, r=self.__overlap_threa)

        # check overlapped bead except for special_bonds
        ind -= len(X)
        flgs = np.empty(len(ind), dtype=bool)
        for i in range(len(ind)):
            flgs[i] = (ind[i] < 0).sum() > 0 or (
                np.abs(ind[i] - i) > self.__special_bonds
            ).sum() > 0
        return flgs.sum() > 0

    def __fix_periodic(self, positions: np.ndarray) -> None:
        mod = positions % self.__cell_length
        q = np.array(positions // self.__cell_length, dtype=int)
        return mod, q

    def __rotation(self, rad: float, n: np.array, vec: np.array) -> np.array:
        cos = np.cos(rad)
        sin = np.sin(rad)
        return vec * cos + (1.0 - cos) * np.dot(vec, n) * n + sin * np.cross(n, vec)


class Writer:
    def __init__(self) -> None:
        today = datetime.datetime.now().strftime("%d %b %Y")
        fname = os.path.basename(__file__)
        self.__contents = f"LAMMPS data file via {fname}, version {today}\n\n"

    def write(self, path: str) -> None:
        f = open(path, mode="w")
        f.write(self.__contents)
        f.close()

    def abouts(self, params: PARAMs) -> None:
        N = params.N
        M = params.M
        nbead_type = len(set(params.bead_type.values()))
        nbtype = len(set(params.btype.values()))
        self.__contents += f"{N * M} atoms\n"
        self.__contents += f"{nbead_type} atom types\n"
        self.__contents += f"{(N - 1) * M} bonds\n"
        self.__contents += f"{nbtype} bond types\n"
        if params.isAngle:
            natype = len(set(params.atype.values()))
            self.__contents += f"{(N - 2) * M} angles\n"
            self.__contents += f"{natype} angle types\n"
        if params.isDihedral:
            ndtype = len(set(params.dtype.values()))
            self.__contents += f"{(N - 3) * M} dihedrals\n"
            self.__contents += f"{ndtype} dihedral types\n\n"

    def cells(self, params: PARAMs) -> None:
        origin = np.array([0.0, 0.0, 0.0])
        if params.origin == "center":
            origin -= params.cell_length * 0.5
        self.__contents += f"{origin[0]} {origin[0] + params.cell_length} xlo xhi\n"
        self.__contents += f"{origin[1]} {origin[1] + params.cell_length} ylo yhi\n"
        self.__contents += f"{origin[2]} {origin[2] + params.cell_length} zlo zhi\n\n"

    def masses(self, params: PARAMs) -> None:
        self.__contents += "Masses\n\n"
        for k, v in params.mass.items():
            self.__contents += f"{k} {v}\n"
        self.__contents += "\n"

    def atoms(self, params: PARAMs, atoms: ATOMS) -> None:
        r = atoms.atom_coordinations
        mir = atoms.mirror_indexes
        self.__contents += "Atoms # molecular\n\n"
        for m in range(params.M):
            for n in range(params.N):
                idx = m * params.N + n
                self.__contents += f"{idx + 1} "
                self.__contents += f"{m + 1} {params.bead_type[n+1]} "
                self.__contents += f"{r[idx, 0]} {r[idx, 1]} {r[idx, 2]} "
                self.__contents += f"{mir[idx, 0]} {mir[idx, 1]} {mir[idx, 2]}\n"
        self.__contents += "\n\n"

    def bonds(self, params: PARAMs) -> None:
        self.__contents += "Bonds\n\n"
        bN = params.N - 1
        for m in range(params.M):
            for n in range(bN):
                b_idx = m * bN + n
                bd_idx = m * params.N + n
                self.__contents += f"{b_idx + 1} {params.btype[n+1]} "
                self.__contents += f"{bd_idx + 1} {bd_idx + 2}\n"
        self.__contents += "\n\n"

    def angles(self, params: PARAMs) -> None:
        if not params.isAngle:
            return
        self.__contents += "Angles\n\n"
        aN = params.N - 2
        for m in range(params.M):
            for n in range(aN):
                a_idx = m * aN + n
                bd_idx = m * params.N + n
                self.__contents += f"{a_idx + 1} {params.atype[n+1]} "
                self.__contents += f"{bd_idx + 1} {bd_idx + 2} {bd_idx + 3}\n"
        self.__contents += "\n\n"

    def dihedrals(self, params: PARAMs) -> None:
        if not params.isDihedral:
            return
        self.__contents += "Dihedrals\n\n"
        dN = params.N - 3
        for m in range(params.M):
            for n in range(dN):
                d_idx = m * dN + n
                bd_idx = m * params.N + n
                self.__contents += f"{d_idx + 1} {params.dtype[n+1]} "
                self.__contents += (
                    f"{bd_idx + 1} {bd_idx + 2} {bd_idx + 3} {bd_idx + 4}\n"
                )
        self.__contents += "\n\n"

    def clear(self) -> None:
        self.__init__()


def Arg():
    parser = argparse.ArgumentParser(description="Create initial position of A-model.")
    parser.add_argument(
        "-in",
        "--input",
        default="input.yaml",
        type=str,
        help="PATH to YAML format input file",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="init.data",
        type=str,
        help="output LAMMPS data file name",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
