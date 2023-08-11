# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 12:30:25 2023

@author: grover.laporte
"""
import pandas as pd
import numpy as np
import plotly.express as px
import tkinter as tk
from tkinter.filedialog import askopenfilename


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
    
    """
    if x<X[0] or x>X[len(X)-1]:
        print("Not a good approximation for this value")
        return(None)
    a = np.array(a);b = np.array(b);c = np.array(c);d = np.array(d)
    X = np.array(X)
    #idx is the spline that is used, np.where returns multiple arguments
    idx = np.where(X[x<=X][0]==X)
    #want to use the index of the previous spline unless it is the first one
    if idx[0][0] != 0:
        idx[0][0] -= 1
    #use the parameters of the correct spline to estimate f(x) using S(x)
    ans = a[idx]+b[idx]*(x-X[idx])+c[idx]*(x-X[idx])**2+d[idx]*(x-X[idx])**3
    return(ans)

##############################################################################
########  Curve Object #######################################################
##############################################################################

class Curve (object):
    """
    A curve is the splines values for a set of discrete data. 
    
    Input:
        -------------------------------------------------------
        row - the row of data (in a pandas Series)
        col - the number of non-numerical columns in the data.
 
    Output:
        -------------------------------------------------------
        A dataframe with the calculus statistics complete
    """
    def __init__(self,row,col):
        self.name = str(row.iloc[:col])
        
        self.labels = list(row.index[col:])
        self.num_labels = np.array([int(c) for c in self.labels])
        self.data = row.iloc[col:]
        self.num_data = self.data.values
        self.spline_coeff()
        
        ## Times need to be cut on the left / right because of derivatives
        self.deltat = 0.1
        stop_time=int(list(row.index)[len(row.index)-1])-0.8
        self.t = np.arange(0.8,stop_time,self.deltat)
        
        ## Calculate the function values and the derivatives
        self.f = self.splines(self.t)
        self.f=self.f.reshape(len(self.f),)
        self.fp = self.fprime(self.t)
        self.fp=self.fp.reshape(len(self.fp),)
        self.fdp = self.fpp(self.t)
        self.fdp=self.fdp.reshape(len(self.fdp),)
        self.baseline = self.num_data[0]
        ## reshaping to make all of the data the same shape
        
    def __str__(self):
        res = ""
        try:
            for label in self.df.columns:
                res += label
            return res
        except:
            return str(self.name)
        
    
    def spline_coeff(self):
        # Use the helper function to record all of the parameter values.
        a,b,c,d = cubic_spline(self.num_labels,self.num_data)
        self.a = a
        self.b = b
        self.c = c
        self.d = d
    
    def splines(self,X):
        """
        The interpolating function. With this code you can pass in a vector or
            single value to estimate the function values using your spline.
        """
        try:
            res = []
            for x in X:
                res.append(S(x,self.num_labels,self.a,self.b,self.c,self.d))
            return np.array(res)
        except:
            return float(S(X,self.num_labels,self.a,self.b,self.c,self.d))

    def fprime(self,X):
        """
        Uses the 5-point Midpoint Formula from Burden/Faires to estimate the
            derivative of the splines. See pp 178. We use a large value for 
            h simply because we do not want to run into round off error if
            the measurements start getting small and the 5-point Midpoint 
            Formula is very accurate for h = 0.01 => on the order of 1E-8 
            which is more than enough.
        """
        h = 1E-2       
        return 1/(12*h)*(self.splines(X-2*h)-8*self.splines(X-h)+
                          8*self.splines(X+h)-self.splines(X+2*h))
    
    def fpp (self,X):
        """
        Same as above.
        """
        h = 1E-2
        return 1/(12*h)*(self.fprime(X-2*h)-8*self.fprime(X-h)+
                        8*self.fprime(X+h)-self.fprime(X+2*h))

    def plot_f(self):
        """
        Plots the curve and the data using plotly which is interactive.
        """
        fig = px.line(x=self.t,y=self.f,color_discrete_sequence=["blue"],
                      title = "Continuous Curve").update_layout(
                      xaxis_title="time", yaxis_title="value")
        fig.add_scatter(x=self.num_labels,y=self.num_data,name = "data",
                        mode = "markers",marker_color="salmon",
                        marker_size=8)
        return fig
    
    def auc(self,time):
        """
        This function allows the user to calculate any AUC value. Although
            we return only values for 180 and 240, there is the possibility to
            change this for an individual user.
        """
        tot = 0
        idx = np.where(self.t<=time)[0]
        for i in idx[:-1]:
            tot += (self.f[i]+self.f[i+1])/2*self.deltat
        AUC = float(tot)
        iAUC = AUC - time*self.baseline
        return AUC,iAUC
    
    def local_opt(self):
        """
        One of the two main functions used to calculate the statistics.
        """
        self.local_max = []
        self.local_min = []
        #t0 is the time at the maximum increase
        t0 = np.where(self.fp==self.fp.max())[0][0]
        #t1 is the time at the maximum decrease
        t1 = np.where(self.fp==self.fp.min())[0][0]
        #use t0, t1 to get the points for max_inc and max_dec
        self.max_inc = [self.t[t0],self.fp[t0]]
        self.max_dec = [self.t[t1],self.fp[t1]]
        #use the helper function to find all of the critical points
        self.critical_points = find_zeros(self.fp,self.f,self.t)
        self.total_critical_pts = len(self.critical_points)
        #use the auc function to find all of the needed auc and iauc values
        self.AUC180, self.iAUC180 = self.auc(180)
        self.AUC240, self.iAUC240 = self.auc(240)
        
        #classify the critical points as max or min based on the 2nd derivative
        for point in self.critical_points:
            if self.fpp(point[0])<=0:
                self.local_max.append(point)
            else:
                self.local_min.append(point)
        
        # turn the list into a numpy array (prob do not need this)
        self.local_max = np.array([np.array(e) for i,e in enumerate(self.local_max)])
        self.local_min = np.array([np.array(e) for i,e in enumerate(self.local_min)])

        
    def inflection_points(self):
        """
        Second of two main functions that calculate the calculus statistics.
        finds inflection points, absolute max, and inflection pt after 
            absolute max
        """
        #find all of the inflection points and put them in self.ips
        self.ips = find_zeros(self.fdp,self.f,self.t)
        #turn the list into a numpy array (not sure if this is needed)
        self.ips = np.array([np.array(e) for i, e in enumerate(self.ips)])
        ex = self.local_max
        
        #find the inflection point after the absolute maximum
        if len(ex)>0 and len(self.ips)>0:
            # putting the first and last points on the vector of local max pts
            # this is because the abs max can actually be an end point.
            ex = np.concatenate([ex,[[self.t[0],self.f[0]]]])
            ex = np.concatenate([ex,[[self.t[-1],self.f[-1]]]])
            
            #find the time of the absolute max
            max_time=ex[np.where(ex[:,1]==ex[:,1].max())[0],:][0]
            
            #determine if an inflection point comes after the max 
            #returns a boolean mask
            bool=self.ips[:,0]>max_time[0]
            self.abs_max = max_time[1]
            self.abs_max_time = max_time[0]
            
        else:
            #When the abs max occurs at the boundary figure out which
            self.abs_max_time = 0 if self.f[0]>self.f[-1] else self.t[-1]
            self.abs_max = self.f[0] if self.f[0]>self.f[-1] else self.f[-1]
        try:
            #if the boolean mask has any True's, then get the first one
            #recall the boolean mask is if an inflection pt occurs after max
            self.ip_after_max = self.ips[bool,:][0]
        except:
            #if there is not one after the max, just return None
            self.ip_after_max = [None,None]
        
        #do not need to return anything since it is not used, use the data
        #structure instead.
        return self.ip_after_max
    
    def calc_stats(self):
        """
        Use the functions in the object class and the helper functions to 
        create a data frame that we can pass out as the solution.
        """
        df = pd.DataFrame([],columns = ["name","baseline","auc180","iauc180","auc240","iauc240",
                                    "max_inc_time","max_inc_val","max_dec_time","max_dec_val",
                                    "abs_max_val","abs_max_time","inf_pt_val","inf_pt_time","tot_crit_pts",
                                    "local_max","local_min"])
        self.local_opt()
        self.inflection_points()
        df.loc[0,"name"]=self.name
        df.loc[0,"baseline"]=self.baseline
        df.loc[0,"auc180"]=self.AUC180
        df.loc[0,"auc240"]=self.AUC240
        df.loc[0,"iauc180"]=self.iAUC180
        df.loc[0,"iauc240"]=self.iAUC240
        df.loc[0,"max_inc_time"]=self.max_inc[0]
        df.loc[0,"max_inc_val"]=self.max_inc[1]
        df.loc[0,"max_dec_time"]=self.max_dec[0]
        df.loc[0,"max_dec_val"]=self.max_dec[1]
        df.loc[0,"abs_max_time"]=self.abs_max_time
        df.loc[0,"abs_max_val"]=self.abs_max
        df.loc[0,"inf_pt_time"]=self.ip_after_max[0]  
        df.loc[0,"inf_pt_val"]=self.ip_after_max[1]
        df.loc[0,"tot_crit_pts"]=self.total_critical_pts
        df.loc[0,"local_max"]=self.local_max
        df.loc[0,"local_min"]=self.local_min
        self.df = df
        return df

def get_file():
    root = tk.Tk()
    root.wm_attributes('-topmost',1)
    root.withdraw()
    
    filename = askopenfilename(title="Please select .csv file",
                               filetypes = (("csv files","*.csv"),
                                            ("text files","*.txt")))
    return filename

def merge_data(col):
    filename = get_file()
    data = pd.read_csv(filename)
    file = filename.split("/")[-1].split(".")[0]
    #folder="/".join(filename.split("/")[:-1])+"/"
    new_file = file+"_calc_stats.csv"
    
    
    c1 = Curve(data.iloc[0,:],col)
    df = c1.calc_stats()
    n = len(data)

    for i in range(1,n):

        nums = [round((i+1)/n,3),round(1-(i+1)/n,3)]
        print(nums)

        c1 = Curve(data.iloc[i,:],col)

        
        df1 = c1.calc_stats()
        df = pd.concat([df,df1],axis=0)

    df.to_csv(new_file,index=False)
    return data,df


if __name__ == "__main__":
    pass
    #get_file()
        
        
    

