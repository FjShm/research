# Calculation of $G(t),\ G',\ G''$

## $G$
### definition
$$
\begin{equation}
G(t)=\frac{V}{k_{\rm B}T}\langle\Pi_{xy}(t)\Pi_{xy}(0)\rangle
\end{equation}
$$

Take the average for the three directions.

$$
\begin{align}
    G(t)&=\frac{1}{3}\frac{V}{k_{\rm B}T}\left[
        \langle\Pi_{xy}(t)\Pi_{xy}(0)\rangle
        +\langle\Pi_{xz}(t)\Pi_{xz}(0)\rangle
        +\langle\Pi_{yz}(t)\Pi_{yz}(0)\rangle
    \right]\\
    &=\frac{1}{3}\frac{V}{k_{\rm B}T}X\\
    &=\frac{1}{3}\frac{V}{k_{\rm B}T}[(1-r)X+rX]\\
&\left(
\begin{matrix}
    X=\langle\Pi_{xy}(t)\Pi_{xy}(0)\rangle
      +\langle\Pi_{xz}(t)\Pi_{xz}(0)\rangle
      +\langle\Pi_{yz}(t)\Pi_{yz}(0)\rangle\\
    0\le r\le 1
\end{matrix}
\right)
\end{align}
$$


Since the system is isotropic, the normal stress difference should relate to the off-diagonal parts of the stress tensor such as

$$
\frac{1}{4}\left[
\langle N_{xy}(t)N_{xy}(0)\rangle
+\langle N_{xz}(t)N_{xz}(0)\rangle
+\langle N_{yz}(t)N_{yz}(0)\rangle
\right]=
\langle\Pi_{xy}(t)\Pi_{xy}(0)\rangle
+\langle\Pi_{xz}(t)\Pi_{xz}(0)\rangle
+\langle\Pi_{yz}(t)\Pi_{yz}(0)\rangle\ \ (=X)
$$

Substitute this expression for $rX$ above. If $r=2/5$,

$$
\begin{align}
    G(t)&=\frac{1}{3}(1-r)\frac{V}{k_{\rm B}T}\left[
    \langle\Pi_{xy}(t)\Pi_{xy}(0)\rangle
    +\langle\Pi_{xz}(t)\Pi_{xz}(0)\rangle
    +\langle\Pi_{yz}(t)\Pi_{yz}(0)\rangle
    \right]\\
    &\ \ \ \ +\frac{1}{3}r\frac{V}{k_{\rm B}T}\cdot\frac{1}{4}\left[
    \langle N_{xy}(t)N_{xy}(0)\rangle
    +\langle N_{xz}(t)N_{xz}(0)\rangle
    +\langle N_{yz}(t)N_{yz}(0)\rangle
    \right]\\
    &=\frac{1}{5}\frac{V}{k_{\rm B}T}\left[
    \langle\Pi_{xy}(t)\Pi_{xy}(0)\rangle
    +\langle\Pi_{xz}(t)\Pi_{xz}(0)\rangle
    +\langle\Pi_{yz}(t)\Pi_{yz}(0)\rangle
    \right]\\
    &\ \ \ \ +\frac{1}{30}\frac{V}{k_{\rm B}T}\left[
    \langle N_{xy}(t)N_{xy}(0)\rangle
    +\langle N_{xz}(t)N_{xz}(0)\rangle
    +\langle N_{yz}(t)N_{yz}(0)\rangle
    \right]
\end{align}
$$
