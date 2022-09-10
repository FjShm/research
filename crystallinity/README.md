# crystallinity


## crystallinity-A

cite: https://doi.org/10/1021/acsapm.2c00147

$$
\boldsymbol d_n=\frac{\boldsymbol r_{n+1}-\boldsynbol r_n}{|r_{n+1}-\boldsynbol r_n|} \\
\lambda=\frac{1}{2k+1}\left|\sum_{i=-k}^k \boldsymbol d_{n+i}\right|
$$

---

## crystallinity-B

$$
\boldsymbol d_n=\frac{\boldsymbol r_{n+1}-\boldsynbol r_n}{|r_{n+1}-\boldsynbol r_n|} \\
\boldsymbol a'_\alpha=\sum_{i=-k}^k \boldsymbol d_{\alpha+i} \\
\boldsymbol a_\alpha=\frac{\boldsymbol a'_\alpha}{|\boldsymbol a'_\alpha|} \\
C_\alpha=|\boldsymbol a_\alpha\cdot\boldsymbol \varepsilon_{\rm T}| \\
{\rm where}\ \boldsymbol\varepsilon_{\rm T}{\rm :\ elongational\ vector}
$$

Average $C_\alpha$ over all chains.
NOTE: Before do this program, dump file must be rotated by `rotate_dump`.

