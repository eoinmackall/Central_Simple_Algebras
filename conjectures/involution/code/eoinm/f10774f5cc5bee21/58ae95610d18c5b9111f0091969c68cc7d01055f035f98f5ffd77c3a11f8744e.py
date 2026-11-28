import json
from fractions import Fraction as F

# Elements: dict (i,j) -> coeff, i = power of 'a', j = power of 'ell', truncated at total degree <= MAXD.
MAXD = 6

def mul(P, Q, maxd=MAXD):
    R = {}
    for (i,j),c in P.items():
        for (k,l),d in Q.items():
            if i+k <= maxd and j+l <= maxd and (i+k)+(j+l) <= maxd:
                R[(i+k,j+l)] = R.get((i+k,j+l), F(0)) + c*d
    return {t:c for t,c in R.items() if c != 0}

def add(P,Q):
    R = dict(P)
    for t,c in Q.items():
        R[t] = R.get(t, F(0)) + c
    return {t:c for t,c in R.items() if c != 0}

def sub(P,Q):
    R = dict(P)
    for t,c in Q.items():
        R[t] = R.get(t, F(0)) - c
    return {t:c for t,c in R.items() if c != 0}

def const(x):
    return {} if x == 0 else {(0,0): F(x)}

def var_a():  return {(1,0): F(1)}
def var_l():  return {(0,1): F(1)}

def smul(c,P):
    if c == 0: return {}
    return {t: c*d for t,d in P.items()}

def exp_of_linear(alpha, gamma, maxd=MAXD):
    # exp(alpha*a + gamma*ell) = sum_k (alpha*a+gamma*ell)^k / k!  , truncated
    L = {}
    if alpha: L[(1,0)] = F(alpha)
    if gamma: L[(0,1)] = F(gamma)
    term = const(1)   # L^k/k!
    total = const(1)
    for k in range(1, maxd+1):
        term = mul(term, L, maxd)              # L^k (unnormalized part kept exactly since each step multiplies by L)
        # careful: recompute properly below
    # recompute cleanly
    total = {}
    pow_k = const(1)
    fact = 1
    for k in range(0, maxd+1):
        if k > 0:
            pow_k = mul(pow_k, L, maxd)
            fact *= k
        if k == 0:
            pow_k = const(1); fact = 1
        total = add(total, smul(F(1, fact), pow_k))
        if not pow_k and k>0:
            pass
    return total

a = var_a(); ell = var_l()
half = F(1,2)

ch_A        = exp_of_linear(1, 0)                 # e^{a}
ch_AdualL   = exp_of_linear(-1, 1)                # e^{ell - a}
ch_v        = add(ch_A, ch_AdualL)
twist       = exp_of_linear(0, -half)             # e^{-ell/2}
lhs         = mul(twist, ch_v)

beta        = add(a, smul(-half, ell))            # a - ell/2
e_beta      = exp_of_linear(1, -half)
e_negbeta   = exp_of_linear(-1, half)
rhs         = add(e_beta, e_negbeta)

assert lhs == rhs, (sorted(lhs.items()), sorted(rhs.items()))

# Now express e^{beta} + e^{-beta} in the single graded variable beta (deg 2), up to degree 6 in beta (=12 cohomological)
# Truncating at total (a,ell)-degree 6 keeps beta^k for k<=6. Check e^b+e^-b = 2 + b^2 + b^4/12 + b^6/360
def poly1(coeffs):  # coeffs[k] for beta^k
    return {k:c for k,c in enumerate(coeffs) if c != 0}
def mul1(P,Q,maxk=6):
    R = {}
    for i,c in P.items():
        for j,d in Q.items():
            if i+j <= maxk: R[i+j] = R.get(i+j,F(0)) + c*d
    return {t:c for t,c in R.items() if c != 0}
def exp1(maxk=6):
    total = poly1([F(1)])
    for k in range(1, maxk+1):
        f = F(1,1)
        for m in range(2, k+1): f = f*F(1,m)
        total = add(total, poly1([F(0)]*k + [f]))
    return total
# (e^b + e^-b) computed as series, but here we must substitute b = a - ell/2 : instead verify via the two-variable dict:
# compare lhs against 2 + beta^2 + beta^4/12 + beta^6/360  (as two-variable polys, truncated at total degree 6)
b2 = mul(beta,beta); b4 = mul(b2,b2); b6 = mul(b4,b2)
expect = add(add(add(const(2), b2), smul(F(1,12), b4)), smul(F(1,360), b6))
assert lhs == expect, (sorted(lhs.items()), sorted(expect.items()))

# what survives on a threefold: cohomology vanishes above degree 6, i.e. monomials of total (a,ell)-degree > 3 vanish.
def keep_deg3(P):
    return {t:c for t,c in P.items() if sum(t) <= 3}
assert keep_deg3(lhs) == keep_deg3(add(const(2), b2)), (keep_deg3(lhs), keep_deg3(add(const(2),b2)))

result = {"identity_e_ell_over_2_ch_equals_ebeta_e_minus_beta": True,
          "full_expansion_up_to_degree6": {str(sum(k)): str(v) for k,v in sorted(lhs.items(), key=lambda x: (sum(x[0]), x[0]))},
          "truncated_to_threefold_H0_plus_H4": {"constant": "2", "beta^2": "1"},
          "difference_ebeta+e-minus_beta_minus(2+beta2)_starts_at_degree": str(min(sum(t) for t in sub(lhs, add(const(2),b2)))) }
print("RESULT: " + json.dumps(result, sort_keys=True))
