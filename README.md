# research

|directory|explanation|
|---|---|
|`IBI`|Used when calculating CG potential.|
|`R2_n`||
|`compare_with_ideal_chain`||
|`crystallinity`|Calculate the crystallinity from lammpstrj file.|
|`grow_system`||
|`ideal_chain`||
|`rotate_dump`|Correct lammpstrj file when extensional deformation is performed using the UEF method. The direction of elongation should be fixed in the z-axis direction.|


## `R2_n`
Used to check equilibrium. 

$$
\begin{align}
  \frac{\langle R^2(n)\rangle}{n} \\
  R^2(n)=|\boldsymbol r_{i+n}-\boldsymbol r_{i}|^2 \\
  n:\text{ the number of bonds}
 \end{align}
$$


## `compare_with_ideal_chain`
A simple program to compare statistics with those in an ideal chain. Only the following statistics are incorporated.

$$
\begin{align}
  \frac{\langle R^4(N)\rangle}{\langle R^2(N)\rangle^2} \\
  \frac{\langle R^2(N)\rangle}{\langle R_{\rm G}^2\rangle} \\
  N:\text{ the number of beads per chain}
\end{align}
$$


## `crystallinity`
Used to calculate orientation and crystallinity. Since the presentation varies from paper to paper, A, B,... Since the method of representation varies from paper to paper, several were created as in A, B,.... Currently, only A calculates meaningful values. (R4.10.17)

`crystallinity-A` even outputs the visualized lammpstrj files.

Required files(crystallinity-A):
- \*.lammpstrj (without `rotate_dump`)
- rotation.txt

![image](https://user-images.githubusercontent.com/100293098/196086239-32259b5b-1ecb-4060-aefe-1b905113d471.png)



## `grow_system`
A summary of the steps taken to create an equilibrium state for the 100% pure AB model. Various LAMMPS input files and initial placement creation scripts need to be rewritten when creating other coarse-grained models.


## `ideal_chain`
Create a lammpstrj file in an ideal chain (random walk) of arbitrary length and number of It is just a random number file, but it can be used to check the statistics of the ideal chain.

## `rotate_dump`
The lammpstrj file resulting from elongation simulation using the UEF method is simply difficult to see because the elongation direction in the Cartesian coordinate system changes each time the cell is switched if visualized as is. This program rewrites the lammpstrj file so that the elongation direction always faces the same direction.

Required files:
- \*.lammpstrj (dump command)
- rotation.txt (fix command)
