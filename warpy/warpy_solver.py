import numpy as np
import matplotlib.pyplot as plt
from warpy_grid import *

class Solver(object):

    def __init__(self,*args,**kwargs):

        # Attributes
        
        # equation attributes
        self.I = kwargs.get("I")
        self.V = kwargs.get("V")
        self.f = kwargs.get("f")
        self.c = kwargs.get("c")
        self.U_0 = kwargs.get("U_0")
        self.U_L = kwargs.get("U_L")
        self.L = kwargs.get("L")
        self.dt = kwargs.get("dt")
        self.C = kwargs.get("C")
        self.T = kwargs.get("T")
        self.user_action = kwargs.get("user_action") 
        self.version = kwargs.get("version")
        self.stability_safety_factor = kwargs.get("stability_safety_factor") 
            
        # Defaults
        if (self.version is None):
            self.version = 'scalar'
        if (self.stability_safety_factor is None):
            self.stability_safety_factor = 1.0
            
        # Attributes

        
        # Mesh attributes
        self.X0 = kwargs.get("X0")
        self.X1 = kwargs.get("X1")
        self.T0 = kwargs.get("T0")
        self.T1 = kwargs.get("T1")

        #self.grid = self.construct_grid()

        
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



    def solve(self):
        
        """Solve u_tt=c^2*u_xx + f on (0,L)x(0,T]."""
        
        Nt = int(round(self.T/self.dt))
        t = np.linspace(0, Nt*self.dt, Nt+1)   # Mesh points in time
        dx = self.dt*self.c/float(self.C)
        Nx = int(round(self.L/dx))
        x = np.linspace(0, self.L, Nx+1)       # Mesh points in space
        C2 = self.C**2                         # Help variable in the scheme
        # Make sure dx and dt are compatible with x and t
        dx = x[1] - x[0]
        self.dt = t[1] - t[0]
        
        if self.f is None or self.f == 0 :
            self.f = lambda x, t: 0
        if self.V is None or self.V == 0:
            self.V = lambda x: 0

        u     = np.zeros(Nx+1)   # Solution array at new time level
        u_n   = np.zeros(Nx+1)   # Solution at 1 time level back
        u_nm1 = np.zeros(Nx+1)   # Solution at 2 time levels back
        
        import time;  t0 = time.clock()  # Measure CPU time

        # Load initial condition into u_n
        for i in range(0,Nx+1):
            u_n[i] = self.I(x[i])
            
        if self.user_action is not None:
            self.user_action(u_n, x, t, 0)

        # Special formula for first time step
        n = 0
        for i in range(1, Nx):
            u[i] = u_n[i] + self.dt*self.V(x[i]) + \
                   0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                   0.5*self.dt**2*self.f(x[i], t[n])
        u[0] = 0;  u[Nx] = 0
        
        if self.user_action is not None:
            self.user_action(u, x, t, 1)

        # Switch variables before next step
        u_nm1[:] = u_n;  u_n[:] = u

        for n in range(1, Nt):
            # Update all inner points at time t[n+1]
            for i in range(1, Nx):
                u[i] = - u_nm1[i] + 2*u_n[i] + \
                       C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                       self.dt**2*self.f(x[i], t[n])

            # Insert boundary conditions
            u[0] = 0;  u[Nx] = 0
            if self.user_action is not None:
                if self.user_action(u, x, t, n+1):
                    break

            # Switch variables before next step
            u_nm1[:] = u_n;  u_n[:] = u

        cpu_time = time.clock() - t0
        return u, x, t, cpu_time

        
        
    def solve2(self):

        
        """Solve u_tt=(c^2*u_x)_x + f on (0,L)x(0,T]."""

        # --- Compute time and space mesh ---
        Nt = int(round(self.T/self.dt))
        t = np.linspace(0, Nt*self.dt, Nt+1) # Mesh points in time


        # Find max(c) using a fake mesh and adapt dx to C and dt
        if isinstance(self.c, (float,int)):
            c_max = self.c
        elif callable(self.c):
            c_max = max([self.c(x_) for x_ in np.linspace(0, self.L, 101)])
        dx = self.dt*c_max/(self.stability_safety_factor*self.C)
        Nx = int(round(self.L/dx))
        x = np.linspace(0, self.L, Nx+1) # Mesh points in space

        
        
        # Make sure dx and dt are compatible with x and t
        dx = x[1] - x[0]
        self.dt = t[1] - t[0]
        
        # Make c(x) available as array
        if isinstance(self.c, (float,int)):
            self.c = np.zeros(x.shape) + self.c
        elif callable(self.c):
            # Call c(x) and fill array c
            c_ = np.zeros(x.shape)
            for i in range(Nx+1):
                c_[i] = self.c(x[i])
            self.c = c_
            
        q = self.c**2
        C2 = (self.dt/dx)**2
        C2 = 0.75**2
        dt2 = self.dt * self.dt # Help variables in the scheme
        
        # --- Wrap user-given f, I, V, U_0, U_L if None or 0 ---
        if self.f is None or self.f == 0:
            self.f = (lambda x, t: 0) if self.version == 'scalar' else \
                lambda x, t: np.zeros(x.shape)
        if self.I is None or self.I == 0:
            self.I = (lambda x: 0) if self.version == 'scalar' else \
                lambda x: np.zeros(x.shape)
        if self.V is None or self.V == 0:
            self.V = (lambda x: 0) if self.version == 'scalar' else \
                lambda x: np.zeros(x.shape)
        if self.U_0 is not None:
            if isinstance(self.U_0, (float,int)) and self.U_0 == 0:
                self.U_0 = lambda t: 0
        if self.U_L is not None:
            if isinstance(self.U_L, (float,int)) and self.U_L == 0:
                self.U_L = lambda t: 0

                
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
            u_n[i] = self.I(x[i])
        if self.user_action is not None:
            self.user_action(u_n, x, t, 0)

        plt.plot(x,u_n)
        plt.show()
            
        # --- Special formula for the first step ---
        for i in Ix[1:-1]:
            u[i] = u_n[i] + self.dt*self.V(x[i]) + \
                   0.5*C2*(0.5*(q[i] + q[i+1])*(u_n[i+1] - u_n[i]) - \
                           0.5*(q[i] + q[i-1])*(u_n[i] - u_n[i-1])) + \
                           0.5 * dt2 * self.f(x[i], t[0])
            
        i = Ix[0]
        print u[i],u[i-1]

        if self.U_0 is None:
            # Set boundary values (x=0: i-1 -> i+1 since u[i-1]=u[i+1]
            # when du/dn = 0, on x=L: i+1 -> i-1 since u[i+1]=u[i-1])
            ip1 = i+1
            im1 = ip1 # i-1 -> i+1
            u[i] = u_n[i] + self.dt*self.V(x[i]) + \
                   0.5*C2*(0.5*(q[i] + q[ip1])*(u_n[ip1] - u_n[i]) - \
                           0.5*(q[i] + q[im1])*(u_n[i] - u_n[im1])) + \
                           0.5*dt2*self.f(x[i], t[0])
            print u[i],u[i-1]
        else:
            u[i] = self.U_0(self.dt)
            
        i = Ix[-1]
        if self.U_L is None:
            im1 = i-1
            ip1 = im1 # i+1 -> i-1
            u[i] = u_n[i] + self.dt*self.V(x[i]) + \
                   0.5*C2*(0.5*(q[i] + q[ip1])*(u_n[ip1] - u_n[i]) - \
                           0.5*(q[i] + q[im1])*(u_n[i] - u_n[im1])) + \
                           0.5*dt2*self.f(x[i], t[0])
            print u[i],u[i-1]
        else:
            u[i] = self.U_L(self.dt)

            
        if self.user_action is not None:
            self.user_action(u, x, t, 1)
            
        # Update data structures for next step
        #u_nm1[:] = u_n; u_n[:] = u # safe, but slower
        u_nm1, u_n, u = u_n, u, u_nm1
        
        # --- Time loop ---
        for n in It[1:-1]:
            # Update all inner points
            if self.version == 'scalar':
                for i in Ix[1:-1]:
                    u[i] = - u_nm1[i] + 2*u_n[i] + \
                           C2*(0.5*(q[i] + q[i+1])*(u_n[i+1] - u_n[i]) - \
                               0.5*(q[i] + q[i-1])*(u_n[i] - u_n[i-1])) + \
                    dt2*self.f(x[i], t[n])
            elif self.version == 'vectorized':
                u[1:-1] = - u_nm1[1:-1] + 2*u_n[1:-1] + \
                          C2*(0.5*(q[1:-1] + q[2:])*(u_n[2:] - u_n[1:-1]) -
                              0.5*(q[1:-1] + q[:-2])*(u_n[1:-1] - u_n[:-2])) + \
                dt2*self.f(x[1:-1], t[n])

                
            else:
                raise ValueError('version=%s' % self.version)

            
            # Insert boundary conditions
            i = Ix[0]
            if self.U_0 is None:
                # Set boundary values
                # x=0: i-1 -> i+1 since u[i-1]=u[i+1] when du/dn=0
                # x=L: i+1 -> i-1 since u[i+1]=u[i-1] when du/dn=0
                ip1 = i+1
                im1 = ip1
                u[i] = - u_nm1[i] + 2*u_n[i] + \
                       C2*(0.5*(q[i] + q[ip1])*(u_n[ip1] - u_n[i]) - \
                           0.5*(q[i] + q[im1])*(u_n[i] - u_n[im1])) + \
                           dt2*self.f(x[i], t[n])
            else:
                u[i] = self.U_0(t[n+1])
                
            i = Ix[-1]
            if self.U_L is None:
                im1 = i-1
                ip1 = im1
                u[i] = - u_nm1[i] + 2*u_n[i] + \
                       C2*(0.5*(q[i] + q[ip1])*(u_n[ip1] - u_n[i]) - \
                           0.5*(q[i] + q[im1])*(u_n[i] - u_n[im1])) + \
                           dt2*self.f(x[i], t[n])
            else:
                u[i] = self.U_L(t[n+1])
                
            if self.user_action is not None:
                if self.user_action(u, x, t, n+1):
                    break
                
            # Update data structures for next step
            u_nm1, u_n, u = u_n, u, u_nm1
            
        cpu_time = time.clock() - t0
        return cpu_time, hashed_input


