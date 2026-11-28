from fractions import Fraction
from itertools import combinations_with_replacement
from math import factorial
import json, random

# ---------- Part A: symbolic check of ch(v^vee)_(i) = (-1)^i ch(v)_(i) and of the degree-2 product ----------
N = 6  # keep cohomology degrees h^0..h^N on CP^n, n >= 3 (h^j != 0 for j <= n)

def trunc(p):
    p = list(p) + [Fraction(0)]*(N+1-len(p))
    return p[:N+1]

def pmul(a,b):
    c = [Fraction(0)]*(2*N+1)
    for i,x in enumerate(a):
        if x == 0: continue
        for j,y in enumerate(b):
            c[i+j] += x*y
    return trunc(c)

def ech(m):  # ch(O(m)) = exp(m h) on CP^n
    return [Fraction(m**j, factorial(j)) for j in range(N+1)]

def ch_sum(exps):
    return trunc([sum(ech(a)[j] for a in exps) for j in range(N+1)])

# check duality identity on all virtual sums of line bundles tested
fams = []
for r in range(4):
    for tup in combinations_with_replacement(range(4), r):
        fams.append(tup)
bad = 0
for a in fams:
    for b in fams:
        v = trunc([x-y for x,y in zip(ch_sum(a), ch_sum(b))])
        vve = trunc([x-y for x,y in zip(ch_sum([-x for x in a]), ch_sum([-x for x in b]))])
        for i in range(N+1):
            if vve[i] != Fraction((-1)**i)*v[i]:
                bad += 1
assert bad == 0

# check the degree-2 identity (exp(ell)*ch(v^vee))_(1) = r*ell - c1(v)  (coefficient of h), ell = k h
for k in range(0,5):
    for a in fams:
        for b in fams:
            v = trunc([x-y for x,y in zip(ch_sum(a), ch_sum(b))])
            vve = trunc([x-y for x,y in zip(ch_sum([-x for x in a]), ch_sum([-x for x in b]))])
            rhs = pmul(ech(k), vve)
            rank = len(a)-len(b)
            assert rhs[1] == Fraction(rank*k) - v[1], (k, a, b)

# ---------- Part B: full search -- D_ell-invariance implies 2c1 = r*ell, and odd k forces even rank ----------
info = {}
for k in range(1,6):
    ranks = set(); nfound = 0; c1ok = True
    for a in fams:
        for b in fams:
            v = trunc([x-y for x,y in zip(ch_sum(a), ch_sum(b))])
            vve = trunc([x-y for x,y in zip(ch_sum([-x for x in a]), ch_sum([-x for x in b]))])
            if v == pmul(ech(k), vve):
                rank = len(a)-len(b); nfound += 1; ranks.add(rank)
                if 2*v[1] != Fraction(rank*k): c1ok = False
    assert c1ok
    if k % 2 == 1:
        assert all(rr % 2 == 0 for rr in ranks), (k, sorted(ranks))
    info[str(k)] = {"nsolutions": nfound, "ranks": sorted(ranks), "c1_identity": c1ok}

# ---------- Part C: the lattice / Bezout step, in Z^n with a random change of basis (coordinate-free look) ----------
random.seed(12345)
tested = 0
for trial in range(3000):
    n = random.randint(1,5)
    U = [[random.randint(-3,3) for _ in range(n)] for _ in range(n)]
    def apply(M, x):  # change of basis matrix, entries of the class in the new frame
        return [sum(M[i][j]*x[j] for j in range(n)) for i in range(n)]
    w = [random.randint(-6,6) for _ in range(n)]
    r = 2*random.randint(-4,4)+1            # odd rank
    ell = apply(U, [2*wi for wi in w])      # ell = 2w in the lattice
    c   = apply(U, [r*wi for wi in w])      # so that 2c = r*ell
    # Bezout reconstruction, purely from the equation r*ell = 2c in a torsion-free abelian group
    from math import gcd
    g = gcd(abs(r), 2)
    assert g == 1
    s = pow(r, -1, 2)
    assert (r*s - 1) % 2 == 0
    t = (r*s - 1)//2                        # r*s - 2*t = 1
    rec = [2*(s*c[i] - t*ell[i]) for i in range(n)]
    assert rec == ell, (r, s, t, ell, c, rec)
    tested += 1
assert tested == 3000

print("RESULT: " + json.dumps({"partA_dual_and_degree2": "ok",
                               "partB_by_k": info,
                               "partC_bezout_lattice_trials": tested,
                               "partC": "ok"}, sort_keys=True))
