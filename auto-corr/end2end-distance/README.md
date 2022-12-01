# end-to-end vector auto correlation

## Desctiption
$$
\begin{align}
    \frac{\langle R^{(M)}(t_0+t') R^{(M)}(t_0) \rangle}{\langle R^{(M)2}(t_0) \rangle} \\
    R ...{\rm end-to-end distance}
\end{align}
$$

## Variables
|variable|mean|
|---|---|
|$R$|end-to-end distance|


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
0 802.494 1
0.1 792.869 0.988006
        :
        :
```
column1: $t'$, column2: auto-correlation of $\boldsymbol R$, column3: standardization of column3

## Restrictions
- using `read_dump.h`
- using `eigen`
- length of all chains must be same

