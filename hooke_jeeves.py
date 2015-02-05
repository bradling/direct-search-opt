def hooke_jeeves(evalf, x0, s, a=1, r=0.5, kmax=1e5, smin=1e-6):
    import numpy as np
    import scipy as sp
    # Hooke and Jeeves  
    k = 0       # function evaulation counter
    n = x0.size
    
    # first step
    xb = x0
    fxb = evalf(xb); k += 1
    
    x, fx = pattern_search(xb, fxb, s); k += (2*n)
    # note: k incremented by 2n because pattern_search calls evalf() 2n times

    # keep reducing step size and continuing pattern search until success
    while (fx >= fxb) and (s > smin):  
        s = r*s                        # reduce step size
        x, fx = pattern_search(xb, fxb, s); k += (2*n)
    
    # if pattern search succeeded, enter main loop
    delta = x - xb
    xb = x; fxb = fx
    while (k < kmax) and (s > smin):
        # first take acceleration step
        xe = xb + a*delta
        fxe = evalf(xe); k += 1
        # pattern search around xe
        x, fx = pattern_search(xe, fxe, s); k += (2*n)
        if fx < fxb: 
        # patten serach succeeded; take the new point and lop
            delta = x-xe
            xb, fxb = x, fx
        else: 
        # pattern search about xe failed; pattern search around xb
            x, fx = pattern_search(xb, fxb, s); k += (2*n)
            if fx < fxb:
                delta = x-xb
                xb, fxb = x, fx
            else: 
            # patten search about xb failed; reduce s and try again
                s = r*s
            #endif
        #endif
    #endwhile
                
    return xb, fxb, k, s
# end hooke_jeeves()
    


# pattern search function for use with hooke_jeeves()
def pattern_search(x, fx, s):
    ii = 0
    n = x.size
    # loop through each dimension
    while ii < n:
        # define current basis vector
        d = np.zeros(n)
        d[ii] = 1
        
        # look at x +/- s*d, take lowest f value
        y = np.array([x, x + s*d, x-s*d])
        fVals = np.array([fx, evalf(y[1]), evalf(y[2]) ])
        idx = np.argmin(fVals)
        x = y[idx]; fx = fVals[idx]
        ii += 1
    return x, fx
# end pattern_search()
