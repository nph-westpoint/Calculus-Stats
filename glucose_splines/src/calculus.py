import numpy as np

def find_zeros(fp,f,t):
    """
    Function that uses discrete function, derivative, and time values to find 
    the local max and min's of a single variable function
    
    Input:
    ---------------------------------------------
    fp : np.array
        the derivative of f at each time, t -> f'(t)
    f : np.array
        the function value at each time, t -> f(t)
    t : np.array
        the independent variable, time

    Output:
    ---------------------------------------------
    Returns zeros : 2D numpy array
        A list of critical points of the function f.
        [[t0,f(t0)],[t1,f(t1)],...]

    """
    # reshape each nd array into a simple row vector
    fp = fp.reshape(len(fp),)
    f = f.reshape(len(f),)
    x0 = fp[0]
    i=0
    zeros=[]
    for x1 in fp[1:]:
        #test to see if the sign on the derivative changes
        if x0*x1<0:
            #if the sign of the derivative changes, then that means
            #you are at a local max or min => record it
            zeros.append([int((t[i]+t[i+1])/2),round((f[i]+f[i+1])/2,3)])
        i+=1
        x0=x1
    return np.array(zeros)