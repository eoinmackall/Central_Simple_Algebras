Let $X$ be a smooth connected projective complex variety.
For $\ell\in H^2(X,\mathbb Z)$, let $L$ be a topological
complex line bundle with $c_1(L)=\ell$.

Define the integral twisted Hodge subgroup
$$
\mathcal H_\ell =
\{v\in KU^0(X):
e^{-\ell/2}\operatorname{ch}(v)
\in\bigoplus_i H^{i,i}(X)\}.
$$

Let $D_\ell:KU^0(X)\to KU^0(X)$ be the involution
$v\mapsto [L]v^\vee$.

For $\varepsilon\in\{+1,-1\}$, let $GW^\varepsilon(X;L)$
be the Grothendieck group of topological complex vector
bundles equipped with nondegenerate $L$-valued
$\varepsilon$-symmetric bilinear forms, and let
$$
F_\varepsilon:GW^\varepsilon(X;L)\to KU^0(X)
$$
be the forgetful homomorphism.

Define
$$
d_{\mathrm{sd}}(X,\ell)=
\gcd\{\operatorname{rk}(v):
v\in\mathcal H_\ell,\ D_\ell(v)=v\}
$$
and
$$
d_\varepsilon(X,\ell)=
\gcd\{\operatorname{rk}(v):
v\in\mathcal H_\ell\cap\operatorname{im}(F_\varepsilon)\}.
$$

Conjecture: There exist $X$, $\ell$, and $\varepsilon$
as above, with
$$
\ell\notin 2H^2(X,\mathbb Z)+c_1(\operatorname{Pic}(X)),
$$
such that
$$
d_\varepsilon(X,\ell)=2d_{\mathrm{sd}}(X,\ell)>0.
$$
