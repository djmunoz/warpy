import numpy as np


class Grid(object):

    '''
    Holds data structures for a uniform mesh on a hypercube in
    space, plus a uniform mesh in time.
    ======== ==================================================
    Argument Explanation
    ======== ==================================================
    L List of 2-lists of min and max coordinates
    in each spatial direction.
    T Final time in time mesh.
    Nt Number of cells in time mesh.
    dt Time step. Either Nt or dt must be given.
    N List of number of cells in the spatial directions.
    d List of cell sizes in the spatial directions.
    Either N or d must be given.
    ======== ==================================================
    Users can access all the parameters mentioned above, plus
    ''x[i]'' and ''t'' for the coordinates in direction ''i''
    and the time coordinates, respectively.
    Examples:
    >>> from UniformFDMesh import Mesh
    >>>
    >>> # Simple space mesh
    >>> m = Mesh(L=[0,1], N=4)
    >>> print m.dump()
    space: [0,1] N=4 d=0.25
    >>>
    >>> # Simple time mesh
    >>> m = Mesh(T=4, dt=0.5)
    >>> print m.dump()
    time: [0,4] Nt=8 dt=0.5
    >>>
    >>> # 2D space mesh
    >>> m = Mesh(L=[[0,1], [-1,1]], d=[0.5, 1])
    >>> print m.dump()
    space: [0,1]x[-1,1] N=2x2 d=0.5,1
    >>>
    >>> # 2D space mesh and time mesh
    >>> m = Mesh(L=[[0,1], [-1,1]], d=[0.5, 1], Nt=10, T=3)
    >>> print m.dump()
    space: [0,1]x[-1,1] N=2x2 d=0.5,1 time: [0,3] Nt=10 dt=0.3
    '''


    
    def __init__(self,
                 L=None, T=None, t0=0,
                 N=None, d=None,
                 Nt=None, dt=None):

        
        if N is None and d is None:
            # No spatial mesh
            if Nt is None and dt is None:
                raise ValueError('Mesh constructor: either Nt or dt must be given')
            if T is None:
                raise ValueError('Mesh constructor: T must be given')
        if Nt is None and dt is None:
            if N is None and d is None:
                raise ValueError('Mesh constructor: either N or d must be given')
            if L is None:
                raise ValueError('Mesh constructor: L must be given')
            
        # Allow 1D interface without nested lists with one element
        if L is not None and isinstance(L[0], (float,int)):
            # Only an interval was given
            L = [L]
        if N is not None and isinstance(N, (float,int)):
            N = [N]
        if d is not None and isinstance(d, (float,int)):
            d = [d]

            
        # Set all attributes to None
        self.x = None
        self.t = None
        self.Nt = None
        self.dt = None
        self.N = None
        self.d = None
        self.t0 = t0

        
        if N is None and d is not None and L is not None:
            self.L = L
            if len(d) != len(L):
                raise ValueError('d has different size (no of space dim.) from L: %d vs %d', len(d), len(L))
            self.d = d
            self.N = [int(round(float(self.L[i][1] -
                                      self.L[i][0])/d[i]))
                      for i in range(len(d))]
        if d is None and N is not None and L is not None:
            self.L = L
            if len(N) != len(L):
                raise ValueError('N has different size (no of space dim.) from L: %d vs %d', len(N), len(L))
            self.N = N
            self.d = [float(self.L[i][1] - self.L[i][0])/N[i]
                      for i in range(len(N))]
            
        if Nt is None and dt is not None and T is not None:
            self.T = T
            self.dt = dt
            self.Nt = int(round(T/dt))
                    
        if dt is None and Nt is not None and T is not None:
            self.T = T
            self.Nt = Nt
            self.dt = T/float(Nt)
            
        if self.N is not None:
            self.x = [np.linspace(
                self.L[i][0], self.L[i][1], self.N[i]+1)
                      for i in range(len(self.L))]
        if Nt is not None:
            self.t = np.linspace(self.t0, self.T, self.Nt+1)

            
    def get_num_space_dim(self):
        return len(self.d) if self.d is not None else 0
    
    def has_space(self):
        return self.d is not None
    
    def has_time(self):
        return self.dt is not None
    
    def dump(self):
        s = ''
        if self.has_space():
            s += 'space: ' + \
                 'x'.join(['[%g,%g]' % (self.L[i][0], self.L[i][1])
                           for i in range(len(self.L))]) + ' N='
            s += 'x'.join([str(Ni) for Ni in self.N]) + ' d='
            s += ','.join([str(di) for di in self.d])
        if self.has_space() and self.has_time():
            s += ' '
        if self.has_time():
            s += 'time: ' + '[%g,%g]' % (self.t0, self.T) + \
                 ' Nt=%g' % self.Nt + ' dt=%g' % self.dt
        return s
