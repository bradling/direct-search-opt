import numpy as np
import scipy as sp

def line_search(x, fx, s, a, b, k, v, kmax):
# This implementation quits line search after 1 success and 1 failure even if 
# we start searching right and the min is left. This should only happen on the 
# first iteration, as once the directions are updated with gram-schmidt we will
# start off searching in the correct direction.
    i = 0
    n = x.size
    jmax = kmax - k
    j = 0
    # search in each direction
    while (i < n) and (j < jmax):
        
        # these flags are used to determine when we have a success and failure 
        # in each direction. 
        success, failure = False, False
        while (j < jmax) and any([(success == False), (failure == False)]):
            # generate new point
            y = x + s[i]*v[i]
            fy = evalf(y); k += 1
            if fy < fx: 
            # if its better, its a success and increase si
                x, fx = y, fy
                s[i] *= a
                success = True
            else: 
            # if its a failure decrease step size
                failure = True              
                s[i] *= -b
            # endif        
            j += 1
        #endwhile
        i += 1
    #endwhile
    return x, fx, k
#end line_search


def rosenbrock(evalf, x0, s, a=3, b=0.5, kmax=1e5, smin=1e-6):
    # initialize search directions to basis vectors    
    v = np.eye(2)
    k = 0
    sMag = np.linalg.norm(s)
    x = x0
    n = x.size
    fx = evalf(x); k += 1
    
    # main loop - do until we hit max feval's or converged
    while (k < kmax) and (sMag > smin):
        # perform crude line search
        xnew, fxnew, k = line_search(x, fx, s, a, b, k, v, kmax)
        
        # redefine s as a vector of how far we moved in line search
        # the way the book updates s doesnt capture how far we moved in each 
        # direction, so this is modiefied from Rosenbrock as presented in BC.pdf
        s = xnew - x  
        x, fx = xnew, fxnew
        
        # now perform gram-schmidt orthogonalization
        # note this implementation is for 2D ONLY
        u = np.empty([2, 2])
        u[0] = s[0]*v[0] + s[1]*v[1]
        u[1] = s[1]+v[1]
        sMag = np.linalg.norm(u[0,:])
        v[0] = u[0] / sMag
        up = u[1] - np.dot(u[0], v[0]) * v[0]
        v[1] = up / np.linalg.norm(up)
        
        s = sMag * np.ones(2)

            
        # reset directions if we lose orthogonality
        if np.linalg.norm(up) < 1e-10:
            v = np.eye(2)
    #endwhile
            
    return x, fx, k, sMag
#end rosenbrock()


