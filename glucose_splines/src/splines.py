import numpy as np
from collections.abc import Iterable


def cubic_spline(x,a):
    """
    Function to construct cubic spline interpolant S for the curve f defined
        by the discrete points (x0,a0), (x1,a1), ... (xn,an) such that:
            a0=f(x0), a1 = f(x1), ... an = f(xn);
            [the spline passes through each data pt]
            ai=S(xi) for all i, 
            [the derivatives for each section of the spline are equal]
            S_i'(xi)=S_i+1'(xi) for all i,
            [the second derivatives for each section of the spline are equal]
            S_i''(xi)=S_i+1''(xi) for all i.
        See Burden and Faires (2010) pp 144-150 for more details.
    
    Input:
    ---------------------------------------------
    x : in this case the x-values are the time values
    a : the measured value for the discrete data (glucose for example)
    
    Output:
    ---------------------------------------------
    Returns a list of parameters for a cubic function (a,b,c,d) for each 
    section of the curve over the time period. There will be n-1 cubic 
    functions that connect the data over the entire time period. 
        S_i(x) = a_i + b_i(x-xi) + c_i(x-xi)^2 + d_i(x-xi)^3
        
    """
    n = len(x)-1
    x = np.array(x)
    a = np.array(a)
    h = np.zeros(n)
    for i in range(n):
        h[i]=x[i+1]-x[i]
    alpha = np.zeros(n+1)
    for i in range(1,n):
        alpha[i] = 3/h[i]*(a[i+1]-a[i])-3/h[i-1]*(a[i]-a[i-1])
        
    l = np.zeros(n+1)
    mu = np.zeros(n+1)
    z = np.zeros(n+1)
    b = np.zeros(n+1)
    c = np.zeros(n+1)
    d = np.zeros(n+1)
    l[0] = 1
    
    for i in range(1,n):
        l[i] = 2*(x[i+1]-x[i-1])-h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]
        z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i]
    
    l[n-1] = 1
    
    for j in range(n-1,-1,-1):
        c[j] = z[j]-mu[j]*c[j+1]
        b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3
        d[j] = (c[j+1]-c[j])/(3*h[j])
    return(a[:-1],b[:-1],c[:-1],d[:-1])

def S(x,X,a,b,c,d):
    """
    This is the interpolant for the function based on its x-value. 
    
        S_i(x) = a_i + b_i(x-Xi) + c_i(x-Xi)^2 + d_i(x-Xi)^3
    
    The a,b,c,d values come from the previous cubic_spline function.
    
    Input:
    ---------------------------------------------
    x : the value at which you want an estimate of the function using the
        cubic spline.
    X : the x-values of the data for which the cubic spline parameters were
        computed.
    a : the computed parameters for the first parameter in the cubic spline.
    b : the computed parameters for the second parameter in the cubic spline.
    c : the computed parameters for the third parameter in the cubic spline.
    d : the computed parameters for the fourth parameter in the cubic spline.
    
    Output: This is a vectorized version of what we have done previously.
    """
    #ans = np.zeros(len(x))
    vector_idx = np.vectorize(lambda x,break_points: np.argmax(break_points>x),
                          excluded=['break_points'],
                         )
    a = np.array(a);b = np.array(b);c = np.array(c);d = np.array(d)
    X = np.array(X)
    
    #vectorize_idx - gives the index i for each value in x where x[i] is greater than
        # the time_break X; if the index is 0, then it is either too small or too large.
    try:
        idx = vector_idx(x,break_points=X) #gives the index of the left X-value
    except:
        return np.array([])  
    #use the last two values and the last set of cubic splines to calculate the largest value.
    max_value = (a[idx-1]+
                 b[idx-1]*(X[-1]-X[idx-2])+
                 c[idx-1]*(X[-1]-X[idx-2])**2+
                 d[idx-1]*(X[-1]-X[idx-2])**3
                 )
    conditions = [
        (idx==0) & (x<=X[0]),
        (idx==0) & (x>=X[-1]),
        (idx!=0),
    ]
    choices = [
        a[0],
        max_value,
        a[idx-1]+b[idx-1]*(x-X[idx-1])+c[idx-1]*(x-X[idx-1])**2+d[idx-1]*(x-X[idx-1])**3
    ]
    ans = np.select(conditions,choices,np.nan)
    return ans
