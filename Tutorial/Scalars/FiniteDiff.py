import numpy as np
import numpy.linalg as nlg
import scipy.sparse as sp
try:
    from scipy.misc import factorial
except:
    from scipy.special import factorial

## FiniteDiff computes FD matrices
## ------
## Dr = FiniteDiff(x, Nx, kwargs)
## ------
#
# kwargs: [default]
#   Sparse = {[True],False}      # Whether or not to return sparse form
#   Uniform = {[True], False} # If domain is Uniform
#   DiffOrd = integer [2]       # Which derivative (i.e. 1 = 1st, 2 = 2nd)
#   Periodic = {True, [False]} # If domain is periodic
#
## Notes
#
#   Periodic option is only available for Uniform domains.
def FiniteDiff(x, n, Sparse=True, Uniform=True, DiffOrd = 1, Periodic=False):

    #FiniteDiff : Create a finite difference matrix of arbitrary order for an
    #             arbitrary grid.
    n = int(n)
    if len(x) == 3:
        # Using a length of 3 is shorthand. x = [a,b,c] is interpreted as
        # x = linspace(a,b,c). The advantage is that we don't actually need
        # to generate x, so we can save on memory when we want large grids.
        Nx = x[2]
    else:
        Nx = len(x)

    # Adjust for higher order derivatives
    n = n + (DiffOrd - 1)
    if DiffOrd == 0:
        Dx = np.eye((Nx,Nx))
        return Dx

    if Sparse:
        Dx = sp.lil_matrix((Nx, Nx))
    else:
        Dx = np.zeros([Nx, Nx])

    n2 = n//2
    cn2 = int(np.ceil(n/2.0))
    fn2 = int(np.floor(n/2.0))

    if Uniform:
        if len(x) == 3:
            # Using a length of 3 is shorthand. x = [a,b,c] is interpreted as
            # x = linspace(a,b,c). The advantage is that we don't actually need
            # to generate x, so we can save on memory when we want large grids.
            dx = (x[1] - x[0])/float((x[2]-1))
        else:
            dx = x[1] - x[0];
        
        # Deal with boundary issues
        if not(Periodic):
            for i in range(0,cn2):
                A = np.zeros([n+1,n+1])
                for j in range(0,n+1):
                    A[:,j] = np.power(((j-i)*dx)*np.ones([1,n+1]),range(0,n+1))/factorial(range(0,n+1))
                b = np.zeros(n+1)
                b[DiffOrd] = 1
                coeff = nlg.solve(A,b)
                coeff = coeff.conj().transpose()
                Dx[i, 0:n+1] = coeff
        
            for i in range(Nx-cn2,Nx):
                A = np.zeros([n+1,n+1])
                for j in range(Nx-n-1,Nx):
                    A[:,j-Nx+n+1] = np.power(((j-i)*dx)*np.ones([1,n+1]),range(0,n+1))/factorial(range(0,n+1))
                b = np.zeros(n+1)
                b[DiffOrd] = 1
                coeff = nlg.solve(A,b)
                coeff = coeff.conj().transpose()
                Dx[i, Nx-n-1:Nx] = coeff
        
        # Now do the internals.
        A = np.zeros([n+1,n+1])
        if n % 2 == 0: # If even...
            for j in range(-n//2-1, n//2):
                A[:,j+n//2+1] = np.power(((j+1)*dx)*np.ones([1,n+1]),range(0,n+1))/factorial(range(0,n+1))
            b = np.zeros(n+1)
            b[DiffOrd] = 1
            coeff = nlg.solve(A,b)
            coeff = coeff.conj().transpose()
           
            if not(Periodic):
                coeff = np.tile(coeff, [Nx, 1])
                Dx[n//2:Nx-n//2, :] = sp.spdiags(coeff.T, range(0,n+1), Nx-n, Nx).todense()
           
            if Periodic:
                Dx += sp.diags(coeff,range(-n2,n2+1),shape=(Nx,Nx)).todense()
                for j in range(1,n2+1):
                    Dx += sp.diags([coeff[n2-j]],Nx-j,shape=(Nx,Nx)).todense()
                    Dx += sp.diags([coeff[n2+j]],j-Nx,shape=(Nx,Nx)).todense()
        elif n % 2 == 1: # If odd...
            for j in range(-fn2-1,cn2):
                A[:,j+fn2+1] = \
                np.power(j*dx*np.ones([1,n+1]),range(0,n+1))\
                    /factorial(range(0,n+1))
            b = np.zeros(n+1)
            b[DiffOrd] = 1
            coeff = nlg.solve(A,b)
            coeff = coeff.conj().transpose()

            if not(Periodic):
                coeff = np.tile(coeff, [Nx, 1])
                Dx[cn2:Nx-fn2,:] = sp.spdiags(coeff.T, range(0,n+1), Nx-n, Nx).todense()

            if Periodic:
                Dx += sp.diags(coeff,range(-cn2,cn2),shape=(Nx,Nx)).todense()
                for j in range(1,cn2+1):
                    Dx += sp.diags([coeff[cn2-j]],Nx-j,shape=(Nx,Nx)).todense()
                for j in range(-1,-cn2,-1):
                    Dx += sp.diags([coeff[cn2-j]],-j-Nx,shape=(Nx,Nx)).todense()
                    
    else: # If not Uniform
        for i in range(0,Nx):
            if i < np.ceil(n/2.0):
                # Deal with 'left' boundary issues
                # Strict inequality since we index from 0
                A = np.zeros([n+1,n+1])
                for j in range(0,n+1):
                    dx = x[j]-x[i]
                    A[:,j] = (dx**range(0,n+1))/factorial(range(0,n+1))
                b = np.zeros(n+1)
                b[DiffOrd] = 1
                coeff = nlg.solve(A,b)
                coeff = coeff.conj().transpose()
                Dx[i, 0:n+1] = coeff
            elif i >= Nx - np.ceil(n/2.0):
                # Deal with 'right' boundary issues
                # Non-strict because of 0 indexing
                A = np.zeros([n+1,n+1])
                for j in range(Nx-n-1,Nx):
                    dx = x[j]-x[i]
                    A[:,j-Nx+n+1] = (dx**range(0,n+1))/factorial(range(0,n+1));
                b = np.zeros(n+1)
                b[DiffOrd] = 1
                coeff = nlg.solve(A,b)
                coeff = coeff.conj().transpose()
                Dx[i, Nx-n-1:Nx] = coeff
            else:
                # Deal with the internal pieces
                # If n is even, then just use a centred scheme
                if n % 2 == 0:
                    A = np.zeros([n+1,n+1])
                    for j in range(-n//2, n//2+1):
                        dx = x[i+j] - x[i]
                        A[:,j+n//2] = (dx**range(0,n+1))/factorial(range(0,n+1))
                    b = np.zeros(n+1)
                    b[DiffOrd] = 1;
                    coeff = nlg.solve(A,b)
                    coeff = coeff.conj().transpose()
                    Dx[i, i-n//2:i+n//2+1] = coeff
            
                # If n is odd, then bias to which side has the closest point.
                elif n % 2 == 1:
                    
                    FLn = int(np.floor(n/2.0))
                    CLn = int(np.ceil(n/2.0))
                    
                    if abs(x[i+CLn] - x[i]) <= abs(x[i-CLn] - x[i]):
                        # If the right end is closer
                        A = np.zeros([n+1,n+1])
                        for j in range(-FLn,CLn+1):
                            dx = x[i+j] - x[i]
                            A[:,j+FLn] = (dx**range(0,n+1))/factorial(range(0,n+1))
                        b = np.zeros(n+1)
                        b[DiffOrd] = 1
                        coeff = nlg.solve(A,b)
                        coeff = coeff.conj().transpose()
                        Dx[i, i-FLn:i+CLn+1] = coeff
                    
                    else:
                        # If the left end is closer
                        A = np.zeros([n+1,n+1])
                        for j in range(-CLn,CLn):
                            dx = x[i+j] - x[i]
                            A[:,j+CLn] = (dx**range(0,n+1))/factorial(range(0,n+1))
                        b = np.zeros(n+1)
                        b[DiffOrd] = 1
                        coeff = nlg.solve(A,b)
                        coeff = coeff.conj().transpose()
                        Dx[i, i-CLn:i+CLn] = coeff
    
    if Sparse:
        # If we're making a sparse matrix, convert it now into csr form.
        Dx = Dx.tocsr()
    
    return Dx
## ------
