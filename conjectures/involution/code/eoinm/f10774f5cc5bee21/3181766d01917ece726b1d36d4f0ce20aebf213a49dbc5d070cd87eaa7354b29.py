import json
from fractions import Fraction as F

# H^*(X,Q) for X = S x E, S a K3 surface, E an elliptic curve.
# u = p^*lambda (deg 2), w = q^*eta (deg 2).  Relations: u^3 = 0 (H^6(S)=0, dim S = 2),
# w^2 = 0 (H^4(E)=0, dim E = 1), and all classes of cohomological degree > 6 vanish.
# Element = dict {(a,b): F} for u^a w^b, degree 2(a+b); truncate at degree 6 and at a<=2,b<=1.
def mul(x, y):
    r = {}
    for (a, b), c in x.items():
        for (d, e), f in y.items():
            A, B = a + d, b + e
            if A > 2 or B > 1 or A + B > 3:
                continue
            r[(A, B)] = r.get((A, B), F(0)) + c * f
    return {k: v for k, v in r.items() if v != 0}
def scale(x, t): return {k: v * t for k, v in x.items() if v * t != 0}
def add(x, y):
    r = dict(x)
    for k, v in y.items(): r[k] = r.get(k, F(0)) + v
    return {k: v for k, v in r.items() if v != 0}
one, u, w = {(0, 0): F(1)}, {(1, 0): F(1)}, {(0, 1): F(1)}
def expo(x, nmax=12):
    r, term = dict(one), dict(one)
    for n in range(1, nmax + 1):
        term = scale(mul(term, x), F(1, n))
        if not term: break
        r = add(r, term)
    return r
def bydeg(x):
    d = {}
    for (a, b), c in x.items(): d.setdefault(a + b, {})[f"u^{a}w^{b}"] = str(c)
    return {k: d[k] for k in sorted(d)}

ell = add(u, scale(w, F(2)))          # c_1(L)
a   = add(u, w)                       # c_1(A) = p^*lambda + q^*eta
b   = w                               # c_1(B) = q^*eta
assert add(a, b) == ell               # c_1(A) + c_1(B) = c_1(L)  (hypothesis of lemma 7b91e0)

inv = expo(scale(ell, F(-1, 2)))      # e^{-ell/2}
ch  = add(expo(a), expo(b))           # ch([A] + [B])
res = mul(inv, ch)
target = add(scale(one, F(2)), {(2, 0): F(1, 4)})       # 2 + p^*(lambda^2/4)
assert res == target, (res, target)

# independent route: e^{-ell/2}(e^{a} + e^{b}) = e^{u/2} + e^{-u/2} = 2cosh(u/2)
alt = add(expo(scale(u, F(1, 2))), expo(scale(u, F(-1, 2))))
assert alt == target, alt

# relation checks used above
assert mul(w, w) == {} and mul(mul(u, u), u) == {}

# also record that the WRONG candidate [O]+[L] leaves a degree-2 (i.e. H^2) part: -ell/2 + ... = 0?
resOL = mul(inv, add(one, expo(ell)))
assert resOL == add(scale(one, F(2)), scale(mul(ell, ell), F(1, 4)))   # = 2 + ell^2/4 = 2+u^2/4+uw
assert resOL != target                                                 # uw-term present, not in H^{i,i} generally

out = {"c1_additivity": True,
       "twisted_char_v_by_degrees": bydeg(res),
       "equals_2_plus_lambda2_over_4": res == target,
       "cosh_route_agrees": alt == target,
       "u3_zero": mul(mul(u, u), u) == {}, "w2_zero": mul(w, w) == {},
       "OL_candidate_by_degrees": bydeg(resOL),
       "OL_differs_from_target": resOL != target,
       "degrees_present_in_res": sorted({k[0] + k[1] for k in res}),
       "bidegrees_present_in_res": sorted({(2 * k[0] + 0 * k[1] + k[1], 2 * k[1] + k[0]) for k in res})}
print("RESULT: " + json.dumps(out, sort_keys=True))
