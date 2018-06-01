import numpy as np
from warpy_grid import *

class Solver(object):

    def __init__(self):
        

        # Attributes
        # equation attributes
        self.c = kwargs.get("c") # c= c(x) wave speed
        self.f = kwargs.get("f") # f= f(x,t) source term
        
        # Mesh attributes
        self.X0 = kwargs.get("X0")
        self.X1 = kwargs.get("X1")
        self.T0 = kwargs.get("T0")
        self.T1 = kwargs.get("T1")

        self.grid = self.construct_grid()

        
    def construct_grid(self):
        
        # --- Compute time and space mesh ---
        Nt = int(round((self.T1-self.T0)/dt))
        self.t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time

        # determine required spacing in time
        # Find max(c) using a fake mesh and adapt dx to C and dt
        if isinstance(c, (float,int)):
            c_max = c
        elif callable(c):
            c_max = max([c(x_) for x_ in np.linspace(0, L, 101)])
        dx = dt*c_max/(stability_safety_factor*C)
        Nx = int(round(L/dx))
        x = np.linspace(0, L, Nx+1) # Mesh points in space



        

        
    def solve(
            I, V, f, c, U_0, U_L, L, dt, C, T,
            user_action=None, version='scalar',
            stability_safety_factor=1.0):


        
        """Solve u_tt=(c^2*u_x)_x + f on (0,L)x(0,T]."""

        # --- Compute time and space mesh ---
        Nt = int(round(T/dt))
        t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time


        # Find max(c) using a fake mesh and adapt dx to C and dt
        if isinstance(c, (float,int)):
            c_max = c
        elif callable(c):
            c_max = max([c(x_) for x_ in np.linspace(0, L, 101)])
        dx = dt*c_max/(stability_safety_factor*C)
        Nx = int(round(L/dx))
        x = np.linspace(0, L, Nx+1) # Mesh points in space


        
        # Make sure dx and dt are compatible with x and t
        dx = x[1] - x[0]
        dt = t[1] - t[0]
        # Make c(x) available as array
        if isinstance(c, (float,int)):
            c = np.zeros(x.shape) + c
        elif callable(c):
            # Call c(x) and fill array c
            c_ = np.zeros(x.shape)
        for i in range(Nx+1):
            c_[i] = c(x[i])
            c = c_
            
        q = c**2
        C2 = (dt/dx)**2; dt2 = dt*dt # Help variables in the scheme
        # --- Wrap user-given f, I, V, U_0, U_L if None or 0 ---
        if f is None or f == 0:
            f = (lambda x, t: 0) if version == 'scalar' else \
                lambda x, t: np.zeros(x.shape)
        if I is None or I == 0:
            I = (lambda x: 0) if version == 'scalar' else \
                lambda x: np.zeros(x.shape)
        if V is None or V == 0:
            V = (lambda x: 0) if version == 'scalar' else \
                lambda x: np.zeros(x.shape)
        if U_0 is not None:
            if isinstance(U_0, (float,int)) and U_0 == 0:
                U_0 = lambda t: 0
        if U_L is not None:
            if isinstance(U_L, (float,int)) and U_L == 0:
                U_L = lambda t: 0

                
        # --- Allocate memory for solutions ---
        
        u = np.zeros(Nx+1) # Solution array at new time level
        u_n = np.zeros(Nx+1) # Solution at 1 time level back
        u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
        import time; t0 = time.clock() # CPU time measurement
        # --- Valid indices for space and time mesh ---
        Ix = range(0, Nx+1)
        It = range(0, Nt+1)
        # --- Load initial condition into u_n ---
        for i in range(0,Nx+1):
            u_n[i] = I(x[i])
        if user_action is not None:
            user_action(u_n, x, t, 0)

        # --- Special formula for the first step ---
        for i in Ix[1:-1]:
            u[i] = u_n[i] + dt*V(x[i]) + \
                   0.5*C2*(0.5*(q[i] + q[i+1])*(u_n[i+1] - u_n[i]) - \
                           0.5*(q[i] + q[i-1])*(u_n[i] - u_n[i-1])) + \
            0.5*dt2*f(x[i], t[0])
        i = Ix[0]
        if U_0 is None:
            # Set boundary values (x=0: i-1 -> i+1 since u[i-1]=u[i+1]
            # when du/dn = 0, on x=L: i+1 -> i-1 since u[i+1]=u[i-1])
            ip1 = i+1
            im1 = ip1 # i-1 -> i+1
            u[i] = u_n[i] + dt*V(x[i]) + \
                   0.5*C2*(0.5*(q[i] + q[ip1])*(u_n[ip1] - u_n[i]) - \
                           0.5*(q[i] + q[im1])*(u_n[i] - u_n[im1])) + \
            0.5*dt2*f(x[i], t[0])
        else:
            u[i] = U_0(dt)
            
        i = Ix[-1]
        if U_L is None:
            im1 = i-1
            ip1 = im1 # i+1 -> i-1
            u[i] = u_n[i] + dt*V(x[i]) + \
                   0.5*C2*(0.5*(q[i] + q[ip1])*(u_n[ip1] - u_n[i]) - \
                           0.5*(q[i] + q[im1])*(u_n[i] - u_n[im1])) + \
            0.5*dt2*f(x[i], t[0])
        else:
            u[i] = U_L(dt)
            
        if user_action is not None:
            user_action(u, x, t, 1)
        # Update data structures for next step
        #u_nm1[:] = u_n; u_n[:] = u # safe, but slower
        u_nm1, u_n, u = u_n, u, u_nm1
        
        # --- Time loop ---
        for n in It[1:-1]:
            # Update all inner points
            if version == 'scalar':
                for i in Ix[1:-1]:
                    u[i] = - u_nm1[i] + 2*u_n[i] + \
                           C2*(0.5*(q[i] + q[i+1])*(u_n[i+1] - u_n[i]) - \
                               0.5*(q[i] + q[i-1])*(u_n[i] - u_n[i-1])) + \
                    dt2*f(x[i], t[n])
            elif version == 'vectorized':
                u[1:-1] = - u_nm1[1:-1] + 2*u_n[1:-1] + \
                          C2*(0.5*(q[1:-1] + q[2:])*(u_n[2:] - u_n[1:-1]) -
                              0.5*(q[1:-1] + q[:-2])*(u_n[1:-1] - u_n[:-2])) + \
                dt2*f(x[1:-1], t[n])

                
            else:
                raise ValueError('version=%s' % version)

            
            # Insert boundary conditions
            i = Ix[0]
            if U_0 is None:
                # Set boundary values
                # x=0: i-1 -> i+1 since u[i-1]=u[i+1] when du/dn=0
                # x=L: i+1 -> i-1 since u[i+1]=u[i-1] when du/dn=0
                ip1 = i+1
                im1 = ip1
                u[i] = - u_nm1[i] + 2*u_n[i] + \
                       C2*(0.5*(q[i] + q[ip1])*(u_n[ip1] - u_n[i]) - \
                           0.5*(q[i] + q[im1])*(u_n[i] - u_n[im1])) + \
                           dt2*f(x[i], t[n])
            else:
                u[i] = U_0(t[n+1])
                
            i = Ix[-1]
            if U_L is None:
                im1 = i-1
                ip1 = im1
                u[i] = - u_nm1[i] + 2*u_n[i] + \
                       C2*(0.5*(q[i] + q[ip1])*(u_n[ip1] - u_n[i]) - \
                           0.5*(q[i] + q[im1])*(u_n[i] - u_n[im1])) + \
                           dt2*f(x[i], t[n])
            else:
                u[i] = U_L(t[n+1])
                
            if user_action is not None:
                if user_action(u, x, t, n+1):
                    break
                
            # Update data structures for next step
            u_nm1, u_n, u = u_n, u, u_nm1
            
        cpu_time = time.clock() - t0
        return cpu_time, hashed_input


