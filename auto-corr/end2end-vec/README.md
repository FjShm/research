# end-to-end vector auto correlation

## Desctiption
$$
\langle \boldsymbol R^{(M)}(t_0+t')\boldsymbol R^{(M)}(t_0) \rangle
$$

## Variables
|variable|mean|
|---|---|
|$\boldsymbol R$|end-to-end vecgtor|


## Input
|key|explanation|other info|
|---|---|---|
|`input_dump_path`|path to `*.lammpstrj` file|(MUST)|
|`output_path`|path to output file (new create)|(MUST)|
|`beads_per_chain`|same as $N$|(MUST)|
|`num_chain`|# of chains (same as $M$)|(MUST)|
|`dt_fs`|timestep dt, units: [fs]|(MUST)|
|`max_frames_to_average`|how many frames to use to do time average|arbitrary (default = 100)|
|`total_frames_dump`|# of frames "input_dump_path" has|arbitrary (default = -1)|


## Output
```
0 800 1.00
0.001 697 0.87125
    :
    :
```
column1: $t'$, column2: auto-correlation of $\boldsymbol R$, column3: standardization of column3

## Restrictions
- using `read_dump.h`
- using `eigen`
- length of all chains must be same

