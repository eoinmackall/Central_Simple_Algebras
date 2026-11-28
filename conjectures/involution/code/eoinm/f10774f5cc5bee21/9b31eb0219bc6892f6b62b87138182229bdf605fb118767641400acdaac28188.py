import json
from fractions import Fraction as F

# Exact arithmetic in R = Q[lam,mu]/(lam^3, mu^2), deg lam = deg mu = 2.
# lam  models p^*(lambda) in H^2(S,Q)  (lam^3 = 0 since S is a surface: H^6(S)=0)
# mu   models q^*(m*eta)  in H^2(E,Q)  (mu^2 = 0 since E is a curve:  H^4(E)=0)
# NOTE (lam+mu)^3 = 3 lam^2 mu != 0 in general, so degree-6 terms must be kept.
def mul(a,b):
    out={}
    for (i1,j1),c1 in a.items():
        for (i2,j2),c2 in b.items():
            i,j=i1+i2,j1+j2
            if i>2 or j>1: continue          # quotient by (lam^3, mu^2)
            out[(i,j)]=out.get((i,j),F(0))+F(c1)*F(c2)
    return {k:v for k,v in out.items() if v}
def sc(a,c): return {k:v*F(c) for k,v in a.items() if v*F(c)}
def add(*ds):
    o={}
    for d in ds:
        for k,v in d.items(): o[k]=o.get(k,F(0))+F(v)
    return {k:v for k,v in o.items() if v}
def neg(a): return {k:-v for k,v in a.items()}
ONE={(0,0):F(1)}; LAM={(1,0):F(1)}; MU={(0,1):F(1)}
def expo(x,n=20):
    t=ONE; s=ONE; f=1
    for k in range(1,n):
        f*=k; t=mul(t,x); s=add(s,sc(t,F(1,f)))
    return s
def dual(a):  # ch(u^vee)_k = (-1)^k ch(u)_k
    return {k:v*((-1)**(k[0]+k[1])) for k,v in a.items()}
def show(a): return {f"lam^{i}mu^{j}":str(v) for (i,j),v in sorted(a.items())}

S=add(LAM,MU)                       # c1(L) = lam + mu   (= "ell")
eS, eSinv = expo(S), expo(sc(S,F(-1,2)))   # e^{ell}, e^{-ell/2}... see below
emhalf    = expo(sc(S,F(-1,2)))             # e^{-ell/2}
emu       = expo(MU)                        # ch([B])
elam      = expo(LAM)                       # ch([B^vee tensor L])

# --- the four summands of v0 = 1 + L + B + B^vee L ---
ch_v0 = add(ONE, eS, emu, elam)
tw    = mul(emhalf, ch_v0)
print("ch(v0)          ", show(ch_v0))
print("e^{-ell/2} ch(v0)", show(tw))

# --- correct per-pair values (they are NOT equal to each other) ---
pair1 = mul(emhalf, add(ONE, eS))     # from {1, L}
pair2 = mul(emhalf, add(elam, emu))   # from {B^vee L, B}
print("pair1           ", show(pair1))
print("pair2           ", show(pair2))
claimed = {(0,0):F(2),(2,0):F(1,4),(1,1):F(1,2)}    # 2 + lam^2/4 + lam*mu/2
claimed2= {(0,0):F(2),(2,0):F(1,4),(1,1):F(-1,2)}   # 2 + lam^2/4 - lam*mu/2
assert pair1 == claimed, (show(pair1), show(claimed))
assert pair2 == claimed2, (show(pair2), show(claimed2))
# neither pair equals 2*cosh(lam/2)*cosh(mu/2) = 2 + lam^2/4 :
coshcosh = {(0,0):F(2),(2,0):F(1,4)}
assert pair1 != coshcosh and pair2 != coshcosh
# and their sum, i.e. the total twisted character, is 4 + lam^2/2 :
tot = add(pair1,pair2)
assert tot == tw == {(0,0):F(4),(2,0):F(1,2)}
for k in [(1,0),(0,1),(1,1),(2,1)]:
    assert tw.get(k,F(0))==0, (k,tw)

# --- ch(1+[L]) twisted: has a nonzero lam*mu/2 term (so NOT purely (i,i) in general) ---
tw1 = mul(emhalf, add(ONE,eS))
assert tw1 == claimed and tw1[(1,1)] == F(1,2)

# --- D_ell(v0) = v0 checked at the level of Chern characters: e^{ell}*ch(v0)^vee = ch(v0) ---
assert mul(eS, dual(ch_v0)) == ch_v0
# rank of v0 and of its negative
rank_v0 = ch_v0[(0,0)]
assert rank_v0 == F(4)

res = {
 "twisted_ch_v0": show(tw),
 "equals_4_plus_lam2_over_2": tw == {(0,0):F(4),(2,0):F(1,2)},
 "pair1": show(pair1), "pair2": show(pair2),
 "pairs_lamu_coeffs": [str(pair1[(1,1)]), str(pair2[(1,1)])],
 "pairs_not_equal_2coshcosh": True,
 "ch_1plusL_twisted": show(tw1),
 "D_invariance_at_ch_level": True,
 "rank_v0": str(rank_v0),
}
print("RESULT: "+json.dumps(res, sort_keys=True))
