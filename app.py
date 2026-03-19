# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 09:24:32 2023

@author: grover.laporte
"""
### Section 1 #### 
### packages ####
import pandas as pd
import numpy as np
import plotly.express as px
import streamlit as st

st.write("""
          # Kondonculator 2.0
          Input your csv file with times (top row needs to be the times that the data was collected)
          as columns and observations as rows. Ensure you have the correct number of non-numerical 
          columns.
          """)

col = st.number_input(":red[Number of non-numerical columns:]",
                      min_value=0,max_value=20,value=2,step=1,
                help="The number of columns used to identify the observation")
st.divider()
## Section 2 ###
## get data ###

uploaded_file = st.file_uploader("Choose a file")
if uploaded_file is not None:
    df = pd.read_csv(uploaded_file)
# df = pd.read_csv("./../NPH_Research//20230622_Codonkulator/data/data_dev.csv")

### Section 3 ###
### functions / classes ###

def find_zeros(fp,f,t):
    """
    Function that uses the derivative to find the local max and mins of a
    single variable function
    Parameters
    ----------
    fp : np.array
        the derivative of f at each time, t -> f'(t)
    f : np.array
        the function value at each time, t -> f(t)
    t : np.array
        the independent variable, time

    Returns
    -------
    zeros : 2D list
        A list of critical points of the function f.
        [[t0,f(t0)],[t1,f(t1)],...]

    """
    fp = fp.reshape(len(fp),)
    f = f.reshape(len(f),)
    x0 = fp[0]
    i=0
    zeros = [[0,0]]
    for x1 in fp[1:]:
        if x0*x1<0:
            zeros.append([int((t[i]+t[i+1])/2),round((f[i]+f[i+1])/2,3)])
        i+=1  
        x0=x1
    zeros.pop(0)
    return (zeros)

def cubic_spline(x,a):
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
    if x<X[0] or x>X[len(X)-1]:
        print("Not a good approximation for this value")
        return(None)
    a = np.array(a);b = np.array(b);c = np.array(c);d = np.array(d)
    X = np.array(X)
    idx = np.where(X[x<=X][0]==X)
    if idx[0][0] != 0:
        idx[0][0] -= 1
    ans = a[idx]+b[idx]*(x-X[idx])+c[idx]*(x-X[idx])**2+d[idx]*(x-X[idx])**3
    return(ans)

@st.cache_data
def merge_df(data,col):
    c1 = Curve(data.iloc[0,:],col)
    df = c1.calc_stats()
    prog_bar = st.progress(value=1/len(data),text="Kodonculating")
    for i in range(1,len(data)):
        prog_bar.progress(value=(i+1)/len(data),text="Kodonculating")
        c1 = Curve(data.iloc[i,:],col)
        df1 = c1.calc_stats()
        df = pd.concat([df,df1],axis=0)
    prog_bar.progress(value=0,text="Complete")
    return df.to_csv().encode('utf-8')

class Curve (object):
    """
    A curve is the splines values for a set of discrete data. 
    Input:  row - the row of data (in a pandas Series)
            col - the number of non-numerical columns in the data.
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
        
    def display(self):
        """
        display is used to show the user the calculus characteristics of each
        curve which coincide with the graph.
        """
        self.calc_stats()
        labels = self.df.columns
        disp = {}
        for lbl in labels[1:15]:
            disp[lbl]=self.df[lbl]
        dspdf = pd.DataFrame(disp)
        dspdf.index = ["stats"]
        st.write(dspdf)

        disp = {}
        disp["local max"]={}
        for i in range(len(self.local_max)):
            disp["local max"][i]={}
            disp["local max"][i][0]=self.local_max[i][0]
            disp["local max"][i][1]=self.local_max[i][1]
        disp["local min"] = {}
        for i in range(len(self.local_min)):
            disp["local min"][i]={}
            disp["local min"][i][0]=self.local_min[i][0]
            disp["local min"][i][1]=self.local_min[i][1]
        st.write(pd.DataFrame(disp))
        
    
    def spline_coeff(self):
        a,b,c,d = cubic_spline(self.num_labels,self.num_data)
        self.a = a
        self.b = b
        self.c = c
        self.d = d
    
    def splines(self,X):
        try:
            #this allows for a list of numbers or just one
            res = []
            for x in X:
                res.append(S(x,self.num_labels,self.a,self.b,self.c,self.d))
            return np.array(res)
        except:
            return float(S(X,self.num_labels,self.a,self.b,self.c,self.d))

    def fprime(self,X):
        h = 1E-2       
        return 1/(12*h)*(self.splines(X-2*h)-8*self.splines(X-h)+
                          8*self.splines(X+h)-self.splines(X+2*h))
    
    def fpp (self,X):
        h = 1E-2
        return 1/(12*h)*(self.fprime(X-2*h)-8*self.fprime(X-h)+
                        8*self.fprime(X+h)-self.fprime(X+2*h))

    def plot_f(self):
        fig = px.line(x=self.t,y=self.f,color_discrete_sequence=["blue"],
                      title = "Continuous Curve").update_layout(
                      xaxis_title="time", yaxis_title="value")
        fig.add_scatter(x=self.num_labels,y=self.num_data,name = "data",
                        mode = "markers",marker_color="salmon",
                        marker_size=8)
        return fig
    
    def auc(self,time):
        tot = 0
        idx = np.where(self.t<=time)[0]
        for i in idx[:-1]:
            tot += (self.f[i]+self.f[i+1])/2*self.deltat
        AUC = float(tot)
        iAUC = AUC - time*self.baseline
        return AUC,iAUC
    
    def local_opt(self):
        self.local_max = []
        self.local_min = []
        t0 = np.where(self.fp==self.fp.max())[0][0]
        t1 = np.where(self.fp==self.fp.min())[0][0]
        self.max_inc = [self.t[t0],self.fp[t0]]
        self.max_dec = [self.t[t1],self.fp[t1]]
        self.critical_points = find_zeros(self.fp,self.f,self.t)
        self.total_critical_pts = len(self.critical_points)
        self.AUC180, self.iAUC180 = self.auc(180)
        self.AUC240, self.iAUC240 = self.auc(240)

        for point in self.critical_points:
            if point[0] < self.t[0]:
                point[0]=self.t[0]
            if self.fpp(point[0])<=0:
                self.local_max.append(point)
            else:
                self.local_min.append(point)

        self.local_max = np.array([np.array(e) for i,e in enumerate(self.local_max)])
        self.local_min = np.array([np.array(e) for i,e in enumerate(self.local_min)])

        
    def inflection_points(self):
        """
        finds inflection points, absolute max, and inflection pt after 
            absolute max
        """
        self.ips = find_zeros(self.fdp,self.f,self.t)
        self.ips = np.array([np.array(e) for i, e in enumerate(self.ips)])
        ex = self.local_max
        
        if len(ex)>0 and len(self.ips)>0:
            ex = np.concatenate([ex,[[self.t[0],self.f[0]]]])
            ex = np.concatenate([ex,[[self.t[-1],self.f[-1]]]])
            max_time=ex[np.where(ex[:,1]==ex[:,1].max())[0],:][0]
            bool=self.ips[:,0]>max_time[0]
            self.abs_max = max_time[1]
            self.abs_max_time = max_time[0]
            
        else:
            self.abs_max_time = 0 if self.f[0]>self.f[-1] else self.t[-1]
            self.abs_max = self.f[0] if self.f[0]>self.f[-1] else self.f[-1]
        try:
            self.ip_after_max = self.ips[bool,:][0]
        except:
            self.ip_after_max = [None,None]
        return self.ip_after_max
    
    def calc_stats(self):
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

try:
    st.download_button(label="Download csv",
                       data = merge_df(df,col),
                       file_name='rename.csv')
    
    st.divider()
    
    row = st.number_input("Select an observation to plot:",
                          min_value=0,max_value=len(df)-1,value=0,step=1,
                    help="Change the observation to see the curve and data.")
    
    c1 = Curve(df.iloc[row,:],col)
    stats = c1.calc_stats()
    st.plotly_chart(c1.plot_f())
    c1.display()
except:
    try:
        df
    except:
        st.write("""### Choose a :blue[**file**] and the :red[**number of non-numerical columns**] in that file.
              """)


        






