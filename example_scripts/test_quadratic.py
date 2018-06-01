import warpy
import numpy as np

def test_quadratic():
    """Check that u(x,t)=x(L-x)(1+t/2) is exactly reproduced."""
    
    def u_exact(x, t):
        return x*(L-x)*(1 + 0.5*t)

    def I(x):
        return u_exact(x, 0)

    def V(x):
        return 0.5*u_exact(x, 0)

    def f(x, t):
        return 2*(1 + 0.5*t)*c**2

    L = 2.5
    c = 1.5
    C = 0.75
    Nx = 6 # Very coarse mesh for this exact test
    dt = C*(L/Nx)/c
    T = 18
    
    def assert_no_error(u, x, t, n):
        u_e = u_exact(x, t[n])
        diff = np.abs(u - u_e).max()
        tol = 1E-13
        print n,t[n],diff, tol
        assert diff < tol

    sol = warpy.Solver(I=I, V=V, f=f, c=c,U_0=None,U_L=None,
                       L=L, dt=dt, C=C, T=T,
                       user_action=assert_no_error)

    sol.solve()


if __name__ == "__main__":

    test_quadratic()
