# mode-p $\boldsymbol X_{\rm p}$ auto-correlation

## Desctiption
$$
\begin{align}
    C(t)&=\frac{\langle\boldsymbol X_{\rm p}^{(m)}(t)\boldsymbol X_{\rm p}^{(m)}(0)\rangle_{m,t}}{\langle\boldsymbol X_{\rm p}^{(m)}(0)^2\rangle_{m,t}} \\
    \boldsymbol X_{\rm p}^{(m)}(t)&=\sum_{n=1}^N \sqrt{\frac{2}{N}}\cos{\left(\frac{(n-1/2)p\pi}{N}\right)\boldsymbol R_n^{(m)}(t)}
\end{align}
$$

## Variables
|variable|mean|
|---|---|
|$\boldsymbol R_n^{(m)}(t)$|position of bead $n$ (of $m$ at $t$)|
|$m$|# of chain|
|$p$|mode (1, 2, ...)|

## Input
|key|explanation|other info|
|---|---|---|
|`input_dump_path`|path to `*.lammpstrj` file|(MUST)|
|`output_path`|path to output file (new create)|(MUST)|
|`beads_per_chain`|same as $N$|(MUST)|
|`num_chain`|# of chains (same as $M$)|(MUST)|
|`dt_fs`|timestep dt, units: [fs]|(MUST)|
|`mode-p`|mode p|(MUST)|
|`max_frames_to_average`|how many frames to use to do time average|arbitrary (default = 100)|
|`total_frames_dump`|# of frames "input_dump_path" has|arbitrary (default = -1)|

## Restrictions
- using `read_dump.h`
- using `eigen`
- length of all chains must be same
