import json, hashlib

K.<z> = CyclotomicField(8)
assert z.multiplicative_order() == 8          # z = exp(pi i/4)
r = [z^(2*a+1) for a in range(4)]
for x in r:
    assert x^4 == -1                          # r_a^4 = -1

# partition p -> the two coordinate pairs used
pairs = {0: ((0,1),(2,3)), 1: ((0,2),(1,3)), 2: ((0,3),(1,2))}

order = [(0,0,0),(0,0,1),(0,0,2),(0,0,3),(0,1,0),(0,1,1),(0,1,2),(0,2,0),(0,2,1),(0,2,2),
         (1,0,0),(1,0,1),(1,0,2),(1,1,0),(1,1,1),(1,1,2),
         (2,0,0),(2,0,1),(2,0,2),(2,1,0)]
assert len(order) == 20 and len(set(order)) == 20
n = 20

R.<x0,x1,x2,x3> = PolynomialRing(K)
Fq = x0^4 + x1^4 + x2^4 + x3^4

def mat(p,a,b):
    (i,j),(k,l) = pairs[p]
    M = matrix(K,2,4,[0]*8)
    M[0,i] = 1; M[0,j] = r[a]
    M[1,k] = 1; M[1,l] = r[b]
    return M

def forms(p,a,b):
    (i,j),(k,l) = pairs[p]
    v = [K(0)]*4
    L1 = None
    row = [K(0)]*4; row[i] = 1; row[j] = r[a]
    L1 = sum(row[t]*R.gens()[t] for t in range(4))
    row2 = [K(0)]*4; row2[k] = 1; row2[l] = r[b]
    L2 = sum(row2[t]*R.gens()[t] for t in range(4))
    return L1, L2

# 1) each chosen line lies on S: substitute the parametrisation into F
T.<u,v> = PolynomialRing(K)
for (p,a,b) in order:
    (i,j),(k,l) = pairs[p]
    xs = [K(0)]*4
    xs[j] = u;  xs[i] = -r[a]*u
    xs[l] = v;  xs[k] = -r[b]*v
    sub = sum((xs[t]*T.gen(0)*0 + (xs[t]))^4*(1) for t in range(4))
    sub = T.zero()
    for t in range(4):
        xt = xs[t]
        sub = sub + xt^4
    assert sub == 0, (p,a,b)
    # the two linear forms are independent (so it really is a P^1)
    assert mat(p,a,b).rank() == 2

# 2) pairwise incidence: rank of the four spanning linear forms
meet = []
skew = []
ranks_ok = True
for I in range(n):
    for J in range(I+1, n):
        M = mat(*order[I]).stack(mat(*order[J]))
        rk = M.rank()
        if rk <= 2: ranks_ok = False     # would mean the two lines coincide
        if rk == 3: meet.append((I+1, J+1))
        elif rk == 4: skew.append((I+1, J+1))
        else: ranks_ok = False
assert ranks_ok
assert len(meet) + len(skew) == n*(n-1)//2

cert = [(1,2),(1,3),(1,4),(1,5),(1,8),(1,11),(1,15),(1,20),(2,3),(2,4),(2,6),(2,9),(2,12),(2,16),(2,19),
        (3,4),(3,7),(3,10),(3,13),(3,18),(4,14),(4,17),(5,6),(5,7),(5,8),(5,14),(5,19),
        (6,7),(6,9),(6,11),(6,15),(6,18),(7,10),(7,12),(7,16),(7,17),(8,9),(8,10),(8,13),(8,18),
        (9,10),(9,14),(9,17),(10,11),(10,15),(10,20),(11,12),(11,13),(11,14),(11,17),
        (12,13),(12,15),(12,18),(12,20),(13,16),(13,19),(14,15),(14,16),(14,18),(14,20),
        (15,16),(15,19),(17,18),(17,19),(17,20),(18,19)]
cert = sorted(set(tuple(sorted(t)) for t in cert))
computed = sorted(meet)
assert cert == computed, (len(cert), len(computed), set(cert)^set(computed))

# 3) intersection matrix: -2 on the diagonal (adjunction), 1 meeting, 0 skew
M20 = matrix(ZZ, n, n)
for i in range(n): M20[i,i] = -2
for (I,J) in computed:
    M20[I-1,J-1] = 1; M20[J-1,I-1] = 1
d = M20.det()
snf = M20.smith_form()
if isinstance(snf, tuple): snf = snf[0]
snf = matrix(ZZ, snf.nrows(), snf.ncols(), [[ZZ(snf[i,j]) for j in range(snf.ncols())] for i in range(snf.nrows())])
diag = [int(snf[i,i]) for i in range(n)]
assert M20.is_symmetric()
rk = M20.change_ring(QQ).rank()
assert d != 0 and rk == n
assert prod(diag) == abs(d)

# 4) cross-check: number of meeting pairs inside/between partitions
def clsp(t): return t[0]
within = sum(1 for (I,J) in computed if order[I-1][0]==order[J-1][0])
cross  = len(computed) - within

def py(o):
    if isinstance(o, dict): return {str(k): py(v) for k,v in sorted(o.items())}
    if isinstance(o, (list, tuple)): return [py(t) for t in o]
    if isinstance(o, (bool,)): return o
    if isinstance(o, int): return int(o)
    if isinstance(o, str): return o
    return int(o)

out = {
 "n_lines": int(n),
 "det": int(d),
 "snf_diag": [int(t) for t in diag],
 "rank_over_QQ": int(rk),
 "n_meet": int(len(computed)),
 "n_meet_hash": hashlib.sha256(str(computed).encode()).hexdigest()[:16],
 "within_same_partition": int(within),
 "cross_partition": int(cross),
 "lines_on_S": "asserted",
 "all_pairs_rank_ge_3": True,
 "matrix_offdiag_entries": sorted(set(int(M20[i,j]) for i in range(n) for j in range(n) if i!=j)),
}
print("RESULT: " + json.dumps(py(out), sort_keys=True))
