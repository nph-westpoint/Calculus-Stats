import numpy as np
import pandas as pd
from .splines import cubic_spline,S
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid
from scipy.signal import find_peaks
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

#from IPython.display import display


def cs_app(x):
    """
    helper function used to `apply` a row of dataframe 
        data to the cubic splines function
    """
    a = np.array(x.index)
    b = x.values
    return cubic_spline(a,b)

def consecutive_counts(x):
    """
    helper function to calculate consecutive counts 
        using `apply`.
    """
    is_null = x.isna()
    groups = is_null.ne(is_null.shift()).cumsum()
    consecutive_counts = is_null[is_null].groupby(groups).size()
    return consecutive_counts.max()


def interpolate_data(y):
    """
    interpolates a row of data using `apply` by first
        extrapolating endpoints, then removing the 
        negative values and replacing with the last 
        non-negative values.
    """
    t = np.array(y.index)
    idx = ~np.isnan(y)
    ## use linear interpolation with extrapolation
    f_linear = interp1d(t[idx],y[idx],kind = 'linear',fill_value = 'extrapolate')
    y1 = f_linear(t)
    ## if any interpolated values are not greater than zero, then interpolate using endpoints
    idx = np.where(y1>0)[0]
    y2 = np.interp(t,t[idx],y1[idx])
    return y2

class Curve (object):
    """
    A curve is the splines values for a set of discrete data. 
    
    Input:
        -------------------------------------------------------
        data - the dataframe with times in the columns and 
            individual meal test observations in the rows.
        col - the number of non-numerical columns in the data.
            these columns are used as the index / description 
 
    Output:
        -------------------------------------------------------
        1) A dataframe with the calculus statistics for each curve 
            (self.stats)
        2) A dataframe with the cubic splines values for 0.1 minutes 
            (self.df)
        3) A dataframe with linear interpolated imputed values (self.data) 
            a dataframe explaining the interpolation (self.missing_data)
        4) An easy way to plot the points / cubic splines curve
        5) An easy way to calculate the derivatives / second derivative of the
            curve.
    """
    def __init__(self,data,col):
        self.dt = 0.1
        cols = list(data.columns[:col])
        self.cols = cols
        self.data = data.copy().set_index(cols)
        self.data.columns =  [int(col) for col in self.data.columns]
        self.time_breaks = np.array(self.data.columns)
        self.t = np.arange(self.time_breaks[0],self.time_breaks[-1],self.dt)
        self.linear_interpolate()
        self.splines = self.data.apply(cs_app,axis = 1)
        self.current_idx = 0
        self.df = self.create_dataframe()
        self.stats = self.calc_stats()

    def linear_interpolate(self):
        """
        linear interpolates missing data for the time series ensuring that
            there are no zeros and records any anomalies in `missing_data`.
        
        """
        df = self.data.copy()
        self.raw_data = self.data.copy()
        cols = self.data.columns

        drop_data_total = df.isna().sum(axis=1)
        drop_data_beginning = df[cols[:2]].isna().sum(axis=1)
        drop_data_end = df[cols[-2:]].isna().sum(axis=1)
        drop_data_consecutive = df.apply(consecutive_counts,axis=1)
        self.missing_data = pd.concat([drop_data_total,
                                       drop_data_consecutive,
                                       drop_data_beginning,
                                       drop_data_end,
                                       ],
                                       axis=1,
                                      )
        self.missing_data.columns = ['total','max_inarow','first_two','last_two']
        df = self.data.apply(interpolate_data,axis=1)
        df = np.array([np.array(df.iloc[i]) for i in range(len(df))])
        df = np.round(df,decimals=1)
        self.data = pd.DataFrame(df,index=self.data.index,columns=self.data.columns)
        return self.data
    
    def __call__(self,idx,t):
        a,b,c,d = self.splines.iloc[idx]
        X = self.time_breaks
        return S(t,X,a,b,c,d)
        
    def fprime(self,idx,X = None,h=1E-3):
        """
        vectorized version of 5 point derivative for the cubic spline
            in the idx(th) row of the data.

        Inputs:
            idx - the row of the data that created the cubic spline
            X - (optional, default is self.t)either a single value 
                or iterable where you want the derivative. 
            h - the step size for the derivative algorithm
        Output:
            single or iterable derivatives of the cubic spline 
        """
        if X is None:
            X = self.t
        X = np.array(X)
        x1,x2,x4,x5 = self(idx,[X-2*h,X-h,X+h,X+2*h])
        return (x1-8*x2+8*x4-x5)/(12*h)
    
    def fdprime(self,idx,X=None,h=1E-3):
        if X is None:
            X = self.t
        X = np.array(X)
        x1,x2,x4,x5 = self.fprime(idx,[X-2*h,X-h,X+h,X+2*h])
        return (x1-8*x2+8*x4-x5)
    
    def find_closest_index(self,vec,val):
        """
        finds the index of the vector (vec) that is the closest
            to the value (val)
        """
        vec = np.array(vec)
        return np.abs(vec-val).argmin()

    def auc(self,x,time):
        """
        helper function to calculate auc using `apply`
        """
        idx = x.name
        idx_t = self.find_closest_index(x.index,time)
        base = self(idx,0)
        auc_ = trapezoid(x.iloc[:idx_t],x.index[:idx_t])
        iauc_ = auc_-(time*base)
        return np.array([auc_,iauc_])
    
    def f(self,x):
        """
        helper function to calculate the cubic spline using `apply`
        """
        idx = x.name
        return self(idx,self.t)

    def create_dataframe(self):
        """
        creates a dataframe for the cubic spline using f
        """
        d1 = (self.data
                   .reset_index()
                   .drop(self.cols,axis=1)
                   .apply(self.f, axis=1)
        )
        d1 = np.array([np.array(d1.iloc[i]) for i in range(len(d1))])
        d1 = np.round(d1,decimals=1)
        self.df = pd.DataFrame(d1,index=np.arange(len(d1)),columns=self.t)
        return self.df
    
    def max_increase_decrease(self,x):
        idx = x.name
        fp = self.fprime(idx,x.index)
        max_i = fp.argmax()
        min_i = fp.argmin()
        return np.array([self.t[min_i],fp[min_i],self.t[max_i],fp[max_i]])
    
    def abs_max(self,x):
        idx = int(x.argmax())
        return np.array([self.t[idx],x.iloc[idx]])
    
    def inflection_after_max(self,x):
        idx = x.name
        max_time = self.stats.loc[idx,'abs_max_time']
        idx_t = self.find_closest_index(x.index,max_time)
        fpp = self.fdprime(idx,x.index)
        fz = np.where(fpp[1:]*fpp[:-1]<0)[0]
        if (fz>idx_t).any():
            try:
                while fz[0]<idx_t:
                    fz = np.delete(fz,0)
            except:
                fz = np.where(fpp[1:]*fpp[:-1]<0)[0]
            return np.array([x.index[fz[0]],self(idx,x.index[fz[0]])])
        else:
            return (np.nan,np.nan)
    
    def critical_points(self,x):
        idx = x.name
        idx_max,_ = find_peaks(x)
        idx_min,_ = find_peaks(-x)
        t_max = np.round(self.t[idx_max],1)
        t_min = np.round(self.t[idx_min],1)
        val_max = np.round(self(idx,t_max),1)
        val_min = np.round(self(idx,t_min),1)
        total_crit_points = len(idx_max)+len(idx_min)
        local_max = [[t,v] for t,v in zip(t_max,val_max)]
        local_min = [[t,v] for t,v in zip(t_min,val_min)]
        return total_crit_points,local_max,local_min
    
    def plot_curve(self,idx,include):
        #fig = go.Figure()
        fig = make_subplots(specs=[[{"secondary_y":True}]])
        if len(include)>0:
            if 'f' in include:
                fig.add_trace(go.Scatter(
                    x = self.t,
                    y = self.df.iloc[idx],
                    mode = 'lines',
                    line = dict(color='blue'),
                    name = 'cubic splines',
                    
                ),secondary_y = False)
            if 'd' in include:
                fig.add_trace(go.Scatter(
                    x = self.time_breaks,
                    y = self.data.iloc[idx],
                    mode = 'markers',
                    marker = dict(color='salmon',size=8),
                    name = 'data',
                    
                ),secondary_y = False)

            if 'fp' in include:
                fig.add_trace(go.Scatter(
                    x = self.t,
                    y = self.fprime(idx,self.t),
                    mode = 'lines',
                    line = dict(color='green'),
                    opacity=0.2,
                    name = '1st derivative',
                    
                ),secondary_y = True)

            if 'fpp' in include:
                fig.add_trace(go.Scatter(
                    x = self.t[1:],
                    y = self.fdprime(idx,self.t[1:]),
                    mode = 'lines',
                    line = dict(color='orange'),
                    opacity=0.2,
                    name = '1st derivative',
                    
                ),secondary_y = True)

            
            
            fig.update_layout(
                title = 'Visualization',
                xaxis_title = "time",
                yaxis_title = "value",
            )
            fig['data'][0]['showlegend']=True
        return fig
                      

    def calc_stats(self):
        self.stats = pd.DataFrame()
        self.stats[['base']] = self.df[[0]]
        if self.time_breaks[-1]>=180:
            self.stats[['auc_180','iauc_180']] = self.df.apply(
                                                self.auc,
                                                axis=1,
                                                time=180,
                                                result_type = 'expand'
                                            )
        if self.time_breaks[-1]>=240:
            self.stats[['auc_240','iauc_240']] = self.df.apply(
                                                self.auc,
                                                axis=1,
                                                time=240,
                                                result_type = 'expand'
                                            )
        if self.time_breaks[-1]<180:
            max_time = self.time_breaks[-1]
            name = 'auc_'+str(self.time_breaks[-1])
            self.stats[[name,'i'+name]]= self.df.apply(
                                                self.auc,
                                                axis=1,
                                                time=max_time,
                                                result_type = 'expand'
                                            )

        
        self.stats[['max_dec_time','max_dec_val','max_inc_time','max_inc_val']] = \
            self.df.apply(
                self.max_increase_decrease,
                axis=1,
                result_type='expand'
            )
        
        self.stats[['abs_max_time','abs_max_val']] = self.df.apply(
                self.abs_max,
                axis=1,
                result_type='expand'
            )
        
        self.stats[['inf_pt_time','inf_pt_val']] = self.df.apply(
            self.inflection_after_max,
            axis=1,
            result_type='expand'
        )

        self.stats[['total_critical_points','local_max','local_min']]=(
            self.df.apply(
                self.critical_points,
                axis=1,
                result_type='expand'
            )
        )
        return self.stats
    def display(self,idx):
        """
        display is used to show the user the calculus characteristics of each
        curve which coincide with the graph.
        """
        cols = list(self.stats.columns)
        min_max_cols = cols[-2:]
        cols = cols[:-2]
        d1 = self.stats[cols].iloc[idx].to_frame()
        d2 = self.stats[min_max_cols].iloc[idx].to_frame()
        return d1,d2








        
        
        
        