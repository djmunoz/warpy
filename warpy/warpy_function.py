class Function(object):
    """
    A scalar or vector function over a mesh (of class Mesh).
    ========== ===================================================
    Argument Explanation
    ========== ===================================================
    mesh Class Mesh object: spatial and/or temporal mesh.
    num_comp Number of components in function (1 for scalar).
    space_only True if the function is defined on the space mesh
    only (to save space). False if function has values
    in space and time.
    ========== ===================================================
    The indexing of ‘‘u‘‘, which holds the mesh point values of the
    function, depends on whether we have a space and/or time mesh.
    Examples:
    >>> from UniformFDMesh import Mesh, Function
    >>>
    >>> # Simple space mesh
    >>> m = Mesh(L=[0,1], N=4)
    >>> print m.dump()
    space: [0,1] N=4 d=0.25
    >>> f = Function(m)
    >>> f.indices
    [’x0’]
    >>> f.u.shape
    (5,)
    >>> f.u[4] # space point 4
    0.0
    >>>
    >>> # Simple time mesh for two components
    >>> m = Mesh(T=4, dt=0.5)
    >>> print m.dump()
    time: [0,4] Nt=8 dt=0.5
    >>> f = Function(m, num_comp=2)
    >>> f.indices
    [’time’, ’component’]
    >>> f.u.shape
    (9, 2)
    >>> f.u[3,1] # time point 3, comp=1 (2nd comp.)
    0.0
    >>>
    >>> # 2D space mesh
    >>> m = Mesh(L=[[0,1], [-1,1]], d=[0.5, 1])
    >>> print m.dump()
    space: [0,1]x[-1,1] N=2x2 d=0.5,1
    >>> f = Function(m)
    >>> f.indices
    [’x0’, ’x1’]
    >>> f.u.shape
    (3, 3)
    >>> f.u[1,2] # space point (1,2)
    0.0
    >>>
    >>> # 2D space mesh and time mesh
    >>> m = Mesh(L=[[0,1],[-1,1]], d=[0.5,1], Nt=10, T=3)
    >>> print m.dump()
    space: [0,1]x[-1,1] N=2x2 d=0.5,1 time: [0,3] Nt=10 dt=0.3
    >>> f = Function(m, num_comp=2, space_only=False)
    >>> f.indices
    [’time’, ’x0’, ’x1’, ’component’]
    >>> f.u.shape
    (11, 3, 3, 2)
    >>> f.u[2,1,2,0] # time step 2, space point (1,2), comp=0
    0.0
    >>> # Function with space data only
    >>> f = Function(m, num_comp=1, space_only=True)
    >>> f.indices
    [’x0’, ’x1’]
    >>> f.u.shape
    (3, 3)
    >>> f.u[1,2] # space point (1,2)
    0.0
    """
    
    def __init__(self, mesh, num_comp=1, space_only=True):
        self.mesh = mesh
        self.num_comp = num_comp
        self.indices = []
        # Create array(s) to store mesh point values
        if (self.mesh.has_space() and not self.mesh.has_time()) or \
           (self.mesh.has_space() and self.mesh.has_time() and \
            space_only):
        # Space mesh only
if num_comp == 1:
self.u = np.zeros(
[self.mesh.N[i] + 1
for i in range(len(self.mesh.N))])
self.indices = [
’x’+str(i) for i in range(len(self.mesh.N))]
else:
self.u = np.zeros(
[self.mesh.N[i] + 1
for i in range(len(self.mesh.N))] +
[num_comp])
self.indices = [
'x'+str(i)
for i in range(len(self.mesh.N))] +\
['component']
if not self.mesh.has_space() and self.mesh.has_time():
# Time mesh only
if num_comp == 1:
self.u = np.zeros(self.mesh.Nt+1)
self.indices = ['time']
else:
# Need num_comp entries per time step
self.u = np.zeros((self.mesh.Nt+1, num_comp))
self.indices = ['time', 'component']
if self.mesh.has_space() and self.mesh.has_time() \
and not space_only:
# Space-time mesh
size = [self.mesh.Nt+1] + \
[self.mesh.N[i]+1
for i in range(len(self.mesh.N))]
if num_comp > 1:
self.indices = ['time'] + \
['x'+str(i)
for i in range(len(self.mesh.N))] +\
['component']
size += [num_comp]
else:
self.indices = ['time'] + ['x'+str(i)
for i in range(len(self.mesh.N))]
self.u = np.zeros(size)
            
