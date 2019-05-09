# Solve 2d poisson equation.

import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import scipy.linalg
import matplotlib.pyplot as plt

def poisson(vort,dx=4400.):
    #print vort
    Nx = len(vort.coord('longitude').points)
    Ny = len(vort.coord('latitude').points)
    #print vort_red
    ####    CREATE DOMAIN   ####
    xmin = ymin = 0.
    xmax = Nx * dx; ymax = Ny * dx

    x = np.linspace(xmin,xmax,Nx)
    y = np.linspace(ymin,ymax,Ny)
    X, Y = np.meshgrid(x,y)
    h = dx
    ####   Create block matrix to solve the poisson equation ####
    block_size = Nx      # i.e. the size of B, I etc - the blocks of the matrix
    nblock = Ny  #i.e. the number of rwos of blocks
    Bdiag = -4 * np.eye(block_size)
    Bupper = np.eye(block_size, k=1)##
    Blower = np.eye(block_size, k=-1)
    B = Bdiag + Blower + Bupper
    
    
    # Neumann BCs: 
    B[0,1] = 2
    B[Nx - 1, Nx - 2] = 2
    
    
    #print B.shape
    #### A should be constructed using block matrix B and identity matrix
    #     B  I  0 ...
    # A = I  B  I ...
    #     0  I  B ...S
    #     ...........

    blst = [B] * nblock   # i.e. create all of the Bs that will form the block diagonal of A
    A = sp.sparse.block_diag(blst)   # Block diagonal of A

## If errors check the following first

    diags = [1] * (block_size * (nblock - 1))  # Identitiy matrices on either side of the main block diagonal
    Dupper = sp.sparse.diags(diags,block_size)
    Dlower = sp.sparse.diags(diags,-block_size)

    # Neuman 0 BCs:
    bc_data = [1] * block_size
    col_inds = [Nx + j for j in np.arange(Nx)]
    row_inds = [j for j in np.arange(Nx)]
    M = block_size * nblock
    top_bc = sp.sparse.csr_matrix((bc_data,(row_inds,col_inds)),shape=(M,M))
    
    col_inds = [(Ny-2)*Nx+j for j in np.arange(Nx)]
    row_inds = [(Ny-1)*Nx+j for j in np.arange(Nx)]
    bottom_bc = sp.sparse.csr_matrix((bc_data,(row_inds,col_inds)),shape=(M,M))
    
    Dupper = Dupper + top_bc
    Dlower = Dlower + bottom_bc
## 

    

    A += Dupper + Dlower # Create A
    

    bblock = vort.data * h**2 # RHS of poisson equation
    b = bblock.flatten() # Turn the matrix into an array so that Ax=b can be solved

### Include Neumann 0 BCs
    ## Top boundary - i = 0, j in (0:Nx-1)
    ## Bottom boundary - i=Ny-1, j in (0:Nx-1)
    #~ for j in np.arange(Nx):
        #~ if A[Nx+j,j] != 1 or A[(Ny-2)*Nx + j, (Ny-1)*Nx + j] != 1:
            #~ print 'Changing wrong value in A for BCs'
        #~ A[Nx+j,j] = 2
        #~ A[(Ny-2)*Nx + j, (Ny-1)*Nx + j] = 2
    #~ ## Left boundary, j = 0, i in (0:Ny-1)
    #~ ## Right boundary, j=Nx-1, i in (0:Ny-1)
    #~ for i in np.arange(Ny):
        #~ A[i*Nx + 1, i*Nx] = 2
        #~ A[(i+1)*Nx - 2, (i+1)*Nx - 1] = 2
    
    u = sp.sparse.linalg.spsolve(A, b)  # Au = b, solve for u
    T = u.reshape((nblock,block_size))

    print A
    print T
    
    psi = vort[:]
    psi.data = T
    
    return psi
    #Tfull = embed(T, Te=0)

    #g = exact_sol(x,y)
    #residue = abs(Tfull - g)
    #x = y = np.linspace(xmin, xmax, Nx)
    #X, Y = np.meshgrid(x,y)



#    plt.show()


def embed(T, Te=0):
    # Embeds T into a matrix giving two extrac colukmns and two extra rows. Boundary conditions currently set to 0
    N = T.shape[0] + 2
    Tfull = np.zeros((N,N))
    Tfull[0] = Te
    Tfull[1:-1, 1:-1] = T
    return Tfull

def f(x,y):
    sigma = 10.; x0=50; y0=50
    x = x[1:-1]; y = y[1:-1]
    X,Y = np.meshgrid(x,y)
    f = (1/(2*np.pi*sigma**6))*np.exp(-((X-x0)**2 + (Y-y0)**2)/(2*sigma**2))*((X-x0)**2+(Y-y0)**2-2*sigma**2)
    #f = 200 * (Y * (Y - 1) + X * (X - 1))
    #f = 100 * X * Y * (1 - X) * (1 - Y)
    return f

def exact_sol(x,y):
    sigma = 10.; x0=50; y0=50
    X,Y = np.meshgrid(x,y)
    g = (1/(2*np.pi*sigma**2))*np.exp(-((X-x0)**2 + (Y-y0)**2)/(2*sigma**2))
    #g = 100 * X * Y * (1 - X) * (1 - Y)
    return g


if __name__ == '__main__':
    main()
