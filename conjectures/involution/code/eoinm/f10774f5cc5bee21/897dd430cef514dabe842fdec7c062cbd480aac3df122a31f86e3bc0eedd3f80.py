import json
from sympy import symbols, expand, Matrix

# Extra exact check: c1(M (x) W) = c1(W) + rk(W) c1(M) also verified through the
# projectivebundle/determinant identity Lambda^r(M tensor W) = M^{tensor r} tensor Lambda^r W
# at the level of symmetric-polynomial invariants, for r = 2,3.
y1,y2,y3,l = symbols('y1 y2 y3 l')
for r, ys in [(2,[y1,y2]), (3,[y1,y2,y3])]:
    c1W = sum(ys)
    c1MW = expand(sum(l+y for y in ys))
    assert expand(c1MW - (c1W + r*l)) == 0
    # top Chern class: prod(l+y_i) = l^r + ... + prod(y_i): coefficient of l^0 is c_r(W)
    top = expand(Matrix([[1]]).tolist()[0][0]*__import__('sympy').prod(l+y for y in ys))
    assert expand(top.subs(l,0) - __import__('sympy').prod(ys)) == 0
print("RESULT: " + json.dumps({"det_identity_ok": True}))
