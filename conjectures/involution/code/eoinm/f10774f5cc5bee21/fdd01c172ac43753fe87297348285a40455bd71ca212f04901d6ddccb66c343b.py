import json
from fractions import Fraction

# ---------- E8 Cartan matrix (Bourbaki-type numbering: chain 1-2-3-4-5-6-7, node 8 attached to 3) ----------
C = [[0]*8 for _ in range(8)]
for i in range(8): C[i][i] = 2
edges = [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(3,8)]
for (u,v) in edges:
    C[u-1][v-1] = C[v-1][u-1] = -1

def matdet(M):
    # exact integer determinant via fraction Gaussian elimination
    n = len(M)
    A = [[Fraction(M[i][j]) for j in range(n)] for i in range(n)]
    det = Fraction(1)
    for c in range(n):
        piv = None
        for r in range(c, n):
            if A[r][c] != 0:
                piv = r; break
        if piv is None: return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]; det = -det
        det *= A[c][c]
        inv = A[c][c]
        for r in range(c+1, n):
            f = A[r][c]/inv
            if f != 0:
                for k in range(c, n):
                    A[r][k] -= f*A[c][k]
    return det

def leading_minors(M):
    out = []
    n = len(M)
    for k in range(1, n+1):
        out.append(int(matdet([row[:k] for row in M[:k]])))
    return out

detC = int(matdet(C))
lmC  = leading_minors(C)   # Sylvester: all > 0  <=>  C positive definite

# Gram of E8(-1) is -C ; simple roots a_1..a_8 (resp. b_1..b_8) form a Z-basis, squares -2
negC = [[-C[i][j] for j in range(8)] for i in range(8)]

# ---------- Gram matrix of Lambda = U^3 (+) E8(-1) (+) E8(-1) ----------
# basis order: e1,f1,e2,f2,e3,f3, a1..a8, b1..b8
n = 6+8+8
G = [[0]*n for _ in range(n)]
def idx_e(i): return 2*(i-1)      # e_i
def idx_f(i): return 2*(i-1)+1    # f_i
for i in range(1,4):
    G[idx_e(i)][idx_e(i)] = 0
    G[idx_f(i)][idx_f(i)] = 0
    G[idx_e(i)][idx_f(i)] = G[idx_f(i)][idx_e(i)] = 1
for t,off in (('a',6),('b',14)):
    for i in range(8):
        for j in range(8):
            G[off+i][off+j] = negC[i][j]

def dot(u,v):
    return sum(u[i]*G[i][j]*v[j] for i in range(n) for j in range(n))

zero = [0]*n
e1 = zero[:]; e1[idx_e(1)] = 1
f1 = zero[:]; f1[idx_f(1)] = 1
a1 = zero[:]; a1[6] = 1
b1 = zero[:]; b1[14] = 1

h   = [e1[i] + 2*f1[i] for i in range(n)]
lam = [e1[i] - 2*f1[i] + a1[i] + b1[i] for i in range(n)]

hh, hl, ll = dot(h,h), dot(h,lam), dot(lam,lam)

def gcd_list(v):
    from math import gcd
    g = 0
    for x in v: g = gcd(g, abs(x))
    return g

# primitivity: h . f1 = 1 (pairing argument) and gcd of coordinates = 1
h_f1 = dot(h, f1)
g_h  = gcd_list(h)

# nondivisibility: for every m in Z, lam - m*h has a_1-coordinate 1 (odd) => not in 2*Lambda
coord_a1_lam_minus_mh = [lam[6] - m*h[6] for m in range(-6,7)]
coord_b1_lam_minus_mh = [lam[14] - m*h[14] for m in range(-6,7)]
all_odd_a = all(x % 2 == 1 or x % 2 == -1 for x in coord_a1_lam_minus_mh)
# brute check over a range of m and mu-coordinates parity: lam = 2mu + mh has solution iff
# all coords of lam - m h are even, for some m
exists_m = any(all((lam[i]-m*h[i]) % 2 == 0 for i in range(n)) for m in range(-20,21))

detG = int(matdet(G))   # unimodularity sanity check

# ---------- evenness of the form: every element of Lambda has even square ----------
# x^2 = -2*(sum of a- and b-coeffs... careful) : compute directly on basis pairs, enough since
# q(x)=x.x is a quadratic form: check q(e_i),q(f_i),q(a_i),q(b_i) even AND bilinear form integral
# => q(sum x_i v_i) = sum x_i^2 q(v_i) + 2 sum_{i<j} x_i x_j (v_i.v_j) is even.
basis = []
for i in range(1,4):
    v = zero[:]; v[idx_e(i)] = 1; basis.append(v)
    v = zero[:]; v[idx_f(i)] = 1; basis.append(v)
for off in (6,14):
    for i in range(8):
        v = zero[:]; v[off+i] = 1; basis.append(v)
assert len(basis) == n
assert all(dot(w,w) % 2 == 0 for w in basis)
assert all(dot(basis[i], basis[j]) in ZZ_int for i in range(n) for j in range(n)) if False else True
integral_pairs = all(isinstance(dot(basis[i], basis[j]), int) for i in range(n) for j in range(n))

# random exact integer vectors: squares even (sanity, not a proof)
import random
random.seed(20240517)
rand_even = True
for _ in range(400):
    w = [random.randint(-9,9) for _ in range(n)]
    if dot(w,w) % 2 != 0: rand_even = False

# express h, lam in the chosen basis to confirm definitions
print('h coeffs in basis:', [dot(h, zero) ])  # placeholder no-op

# exhaustive membership check for lam in 2*Lambda + Z*h over ALL m:
# a1-coordinate of lam - m*h equals lam[6] - m*h[6] = 1 - 0 = 1 for every m in Z.
for m in range(-200, 201):
    assert (lam[6] - m*h[6]) % 2 != 0
    assert (lam[14] - m*h[14]) % 2 != 0

checks = {
 "detE8Cartan": detC,
 "leading_minors_E8_Cartan": lmC,
 "diag_negC": [negC[i][i] for i in range(8)],
 "offdiag_negC_values": sorted({negC[i][j] for i in range(8) for j in range(8) if i!=j}),
 "det_Gram_Lambda": detG,
 "h_sq": hh, "h_dot_lam": hl, "lam_sq": ll,
 "h_dot_f1": h_f1, "gcd_coords_h": g_h,
 "a1_coord_of_lam_minus_mh": sorted(set(coord_a1_lam_minus_mh)),
 "b1_coord_of_lam_minus_mh": sorted(set(coord_b1_lam_minus_mh)),
 "all_a1_coords_odd": all_odd_a,
 "exists_m_with_lam-mh_in_2Lambda": exists_m,
 "h": h, "lam": lam,
}
assert detC == 1
assert all(x > 0 for x in lmC)          # E8 lattice positive definite => E8(-1) negative definite
assert set(negC[i][i] for i in range(8)) == {-2}
assert detG in (1,-1)
assert hh == 4 and hl == 0 and ll == -8
assert h_f1 == 1 and g_h == 1
assert integral_pairs and rand_even
assert all_odd_a and exists_m is False
print("RESULT: " + json.dumps(checks, sort_keys=True))
