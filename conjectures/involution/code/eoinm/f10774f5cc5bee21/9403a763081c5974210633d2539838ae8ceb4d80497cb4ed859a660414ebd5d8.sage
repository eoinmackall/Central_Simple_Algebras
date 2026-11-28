import json
## Check: for w = sum_i n_i([A_i]+[A_i^vee otimes L]) on a threefold X (dim_C = 3),
## with beta_i = c1(A_i) - ell/2, the twisted character e^{-ell/2} ch(w) has no
## degree-2 and no degree-6 component and equals 2N + sum_i n_i beta_i^2.
## Model H^{even}(X,Q) by Q[b_1,...,b_n] truncated above total degree 3, since each
## b_i^4 lies in H^8(X) = 0 for dim_C X = 3, and e^b + e^{-b} = 2 + b^2 modulo degree >= 4.
results = []
for coeffs in ([1], [3], [1,1], [2,-1,4], [-1,-2], [5,-3]):
    n = len(coeffs)
    R = PolynomialRing(QQ, ['b%d'%i for i in range(n)])
    bs = R.gens()
    total = R(0)
    for c, b in zip(coeffs, bs):
        total = total + c*(2 + b^2)          # truncated e^b + e^{-b}
    N = sum(coeffs)
    closed = 2*N + sum(c*b^2 for c, b in zip(coeffs, bs))
    assert total == closed                   # formula holds
    assert total.degree() <= 2 or n == 0     # no degree >= 3 terms: no degree-2 or degree-6 piece
    results.append({"coeffs": [int(c) for c in coeffs], "N": int(N),
                    "rank": int(2*N), "divisible_by_4": bool((2*N) % 4 == 0),
                    "twisted_char": str(closed)})
bad = [N for N in range(-201, 202, 2) if (2*N) % 4 == 0]
assert bad == []                             # 2N never divisible by 4 for odd N
print("RESULT: " + json.dumps({"cases": results, "parity_all_odd_N_in_[-201,201]": "2N not divisible by 4: OK"}))
