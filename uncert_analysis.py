import numpy as np
import sympy as sp

# type each parameter here
a = 2603.8918197861994
b = 3449.085703325699
d = 1.3848654292403553
e = -1.0608020168978702
f = 16042.928201531528
g = 2395.4402913385798
Astatic = 16289.764999999996

params = [a,b,d,e,f,g,Astatic]

# type each uncertainty here
da = 4.4643667520845804
db = 4.310164651965979
dd = 0.023946905539360915
de = 0.01874955064022889
df = 53.79008995334951
dg = 50.74465644587862
dAstatic = 129.54215829219518

uncerts = [da,db,dd,de,df,dg,dAstatic]

def S(z,y,x,w,v,u,t):
    # note that this result is in degrees per second per bit.
    return (1/np.sqrt(z**2+y**2)* np.arccos((u*np.cos(np.sqrt(x**2+w**2))+v) / (t))) * 180/np.pi

# a list of symbols to represent the independent variables
symbols = a_s,bs,ds,es,fs,gs,As = sp.symbols('a b d e f g A')

# the equation for sensitivity that uses our indep variable symbols:
S_sym = (1/((a_s**2 + bs**2))**(0.5)* sp.acos((gs*sp.cos((ds**2 + es**2)**(0.5)) + fs) / (As))) * 180/np.pi

def grad_func(a,b,d,e,f,g,A):
    # returns the vector of numerical partial derivatives (i.e., the gradient)
    sym_partial_derivs = [sp.diff(S_sym, var) for var in symbols]
    grad = []
    for func in sym_partial_derivs:
        func = sp.lambdify(symbols, func, 'numpy')
        grad.append(func(a,b,d,e,f,g,A))
    return grad

gradient = grad_func(*params)

variance = 0 # initialize it at 0, then sum it over all the variables and uncerts
for i in range(len(gradient)):
    variance += (gradient[i] * uncerts[i])**2

print('Sensitivity:', S(*params)*1000, 'mdps/bit')
print('Uncertainty in Sensitivity:', np.sqrt(variance)*1000, 'mdps/bit')
