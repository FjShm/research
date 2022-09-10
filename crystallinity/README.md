# crystallinity


## crystallinity-A

cite: https://doi.org/10/1021/acsapm.2c00147

$$
\begin{align}
\boldsymbol d_n&=\frac{\boldsymbol r_{n+1}-\boldsymbol r_n}{|r_{n+1}-\boldsymbol r_n|} \\
\lambda&=\frac{1}{2k+1}\left|{\sum_{i=-k}^k \boldsymbol d_{n+i}}\right|
\end{align}
$$

---

## crystallinity-B

$$
\begin{align}
{\boldsymbol d}_n&=\frac{{\boldsymbol r}_{n+1}-{\boldsymbol r}_n}{|r_{n+1}-{\boldsymbol r}_n|} \\
{\boldsymbol a}'_\alpha&=\sum_{i=-k}^k {\boldsymbol d}_{\alpha+i} \\
{\boldsymbol a}_{\alpha}&=\frac{{\boldsymbol a}'_\alpha}{|{\boldsymbol a}'_\alpha|} \\
C_\alpha&=|{\boldsymbol a}_\alpha\cdot{\boldsymbol \varepsilon}_{\rm T}| \\
&{\rm where}\ {\boldsymbol\varepsilon}_{\rm T}{\rm :\ elongational\ vector}
\end{align}
$$

Average $C_\alpha$ over all chains.
NOTE: Before do this program, dump file must be rotated by `rotate_dump`.

