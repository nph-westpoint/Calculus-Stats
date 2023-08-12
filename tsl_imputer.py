# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 07:47:28 2023

@author: grover.laporte
"""
"""
Created on Thu Jul 13 14:45:08 2023
    
    Linear interpolation imputer for time-series data
    This is an extension of the work that Skyler Chauff and I have been working
        on last semester (AY23-2).
        
    The next step is to create a ARIMA imputer in order to compare the results.
    
@author: grover.laporte
"""
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter.filedialog import askopenfilename

def get_file():
    root = tk.Tk()
    root.wm_attributes('-topmost',1)
    root.withdraw()
    
    filename = askopenfilename(title="Please select .csv file",
                               filetypes = (("csv files","*.csv"),
                                            ("text files","*.txt")))
    return filename

def linear_interpolation(admin_cols):
    filename = get_file()
    df = pd.read_csv(filename)
    imp = Imputer(df,admin_cols)
    if imp.impute_bool:
        df = imp.impute_data()
    return df

class Imputer(object):
    def __init__(self,data,admin_cols,times=None):
        """
        Imputer Object for the task of imputing the data.
        
        Input
        -------------------------------------------------------------------
        data - pandas DataFrame
        admin_cols - the number of cols in the data that will not be imputed
        times - use the list times if given, otherwise get it from the cols
        
        Output / Stored data from constructor
        -------------------------------------------------------------------
        self.times - a numerical list of times used for each observation
        self.admin_cols - the number of admin cols for each observation
        self.data - the original dataframe with numerical column names
        self.X - orginial dataframe without any admin cols 
        self.shape - the shape of self.X (rows and cols)
        self.impute_bool - is data missing (need to impute): True or False
        self.deleted_rows - blank list
        self.imputed_data = blank list
        
        self.rows - the rows in self.X that have missing data
        self.cols - list of np arrays that corresponds to self.rows (example)
        
        Example: self.rows = [ 2,  3,  5,  7,  9, 10, 14, 18, 19, 20]
                 self.cols = [array([6], dtype=int64),
                              array([0, 1], dtype=int64),
                              array([2, 6], dtype=int64),
                              array([1, 2, 3], dtype=int64),
                              array([0, 3], dtype=int64),
                              array([4], dtype=int64),
                              array([3, 4], dtype=int64),
                              array([0], dtype=int64),
                              array([1, 4], dtype=int64),
                              array([5, 6], dtype=int64)]
        
        row 2 is missing data in column 6
        row 3 is missing data in columns 0 and 1
        row 5 is missing data in columns 2 and 6
        ... and so on
        
        """
        ## Get the times from the list that is given or from the column headers
        ac = admin_cols
        if times != None:
            self.times = times
        else:
            self.times = [int(t) for t in data.columns[ac:]]
        times = self.times
        
        ## Record the number of admin columns
        self.admin_cols = admin_cols
        
        #self.data - original data (including admin) with numerical columns
        data.columns = list(data.columns)[:admin_cols]+times
        self.data = data
        
        #self.X - original data (no admin) with numerical columns
        self.X = pd.DataFrame(data.iloc[:,admin_cols:].values,
                              columns=self.times)
        
        # imputer shape is the shape of the data not including admin
        self.shape = self.X.values.shape

        # rows in the dataframe that have a null value
        rows = self.X.isnull().any(axis=1) #creates a boolean mask
        # does the data need imputing: if there are rows, then True else False
        self.impute_bool = True if rows.sum()>0 else False
        rows = np.where(rows==True)[0] #gets the indices for rows missing data
        # store these rows that have missing values 
        self.rows = rows
        
        #for each row that has a null value, find the cols that are null
        cols = []
        for row in rows:
            col_ = np.where(self.X.iloc[row].isnull()==True)[0]
            cols.append(col_)
        self.cols = cols
        self.deleted_rows = []
        self.imputed_data = []
        
    def consecutive_numbers(self,cols):
        return sorted(cols) == list(range(min(cols),max(cols)+1))
    
    
    def linear(self,df,times):
        """
        method linear >>> used to calculate a missing value given the two
                           nearest points.
        Input: df      >>> the data frame for the two nearest points where 
                            the x-values are in the index and y-values in the
                            values.
               times   >>> list of times with missing values
               
        Output: numpy array for all of the missing times.
        """
        x = df.index
        y = df.values
        pt1 = [x[0],y[0]]; pt2 = [x[1],y[1]]
        m = (pt2[1]-pt1[1])/(pt2[0]-pt1[0])
        b = pt1[1]-m*pt1[0]
        f = lambda t: m*t+b
        ans = []
        for t in times:
            ans.append(int(f(t)))
        return np.array(ans)
    
    def impute_or_not(self,cols):
        """
        use this method to determine whether or not to impute the row
            based on the number of columns that are missing data. More
            is needed to determine how many missing values are too much.
            
            I use a simple formula until that research is done. If there
            are more than three values missing return False, impute otherwise.
        """
        return len(cols)<3
    
    def impute_data(self):
        """
        method is used to impute the data that we want to impute while
            deleting the rows that are missing more than 2 values.
        """
        for i in range(len(self.rows)-1,-1,-1):
            if self.impute_or_not(self.cols[i]):
                row = self.rows[i]
                cols = self.cols[i]
                times,vals = self.values_to_impute(row,cols)
                for t,v in zip(times,vals):
                    self.X.loc[row,t] = v
                if np.any(vals<0):
                    self.deleted_rows.append(row)
                    self.X = self.X.drop(row)
                    self.data = self.data.drop(row)
            else:
                self.deleted_rows.append(self.rows[i])
                self.X = self.X.drop(self.rows[i])
                self.data = self.data.drop(self.rows[i])
                
        data = pd.concat([self.data.iloc[:,:self.admin_cols],
                          self.X],axis=1)
        self.imputed_data = data
        return data
    
    def values_to_impute(self,row,cols):
        """
        this method returns the data points and times of the missing 
            values to pass to the linear imputer.
        Currently, we are passing on data that have more than 2 data points
            missing, so we only check a few possibilities.
            group1 - missing 1 at the start.
            group2 - missing 1 in the middle.
            group3 - missing 1 at the end.
            group4 - missing 2 at the start.
            group5 - missing 2 in the middle.
            group6 - missing 2 at the end.
            group7 - missing 1 at the start and 1 in the middle.
            group8 - missing 1 in the middle and 1 at the end.
            group9 - missing 2 non-consecutive in the middle.
        """
        last_col = self.shape[1]-1
        N = len(cols)
        pts = self.X.iloc[row,:]
        times = np.array(self.times)
        pts_ = 0; t = 0; val = 0; col_times = 0
        if N == 1:
            # only groups 1,2,3
            if 0 in cols:
                #Group1
                t = [times[0]]
                idx = np.array([1,2])
                pts_ = pts.iloc[idx]
                val = self.linear(pts_,t)
                col_times = t
                
            elif last_col in cols:
                #Group3
                t = [times[last_col]]
                idx = np.array([last_col-2,last_col-1])
                pts_ = pts.iloc[idx]
                val = self.linear(pts_,t)
                col_times = t
                
            else:
                #Group2
                col = cols[0]
                t = [times[col]]
                idx = np.array([col-1,col+1])
                pts_ = pts.iloc[idx]
                val = self.linear(pts_,t)
                col_times = t
            
        elif N == 2:
            # can only be in group 4,5,6,7,8,9
            if self.consecutive_numbers(cols):
                # only groups 4,5,6
                if 0 in cols:
                    #group4
                    t = times[cols]
                    idx = np.array([2,3])
                    pts_ = pts.iloc[idx]
                    val = self.linear(pts_,t)
                    col_times = np.array(t)
                elif last_col in cols:
                    #group6
                    t = times[cols]
                    idx = np.array([last_col-3,last_col-2])
                    pts_ = pts.iloc[idx]
                    val = self.linear(pts_,t)
                    col_times = np.array(t)
                else:
                    #group5
                    t = times[cols]
                    idx = np.array([cols[0]-1,cols[1]+1])
                    pts_ = pts.iloc[idx]
                    val = self.linear(pts_,t)
                    col_times = np.array(t)
                    
            elif 0 in cols:
                #group7
                col = cols[1]
                val = []
                col_times=[]
                ########################################
                t = [times[0]]
                idx = np.array([1,2])
                pts_ = pts.iloc[idx]
                val.append(self.linear(pts_,t)[0])
                col_times.append(t[0])

                ########################################
                t = [times[col]]
                idx = np.array([col-1,col+1])
                pts_ = pts.iloc[idx]
                val.append(self.linear(pts_,t)[0])
                col_times.append(t[0])
                ########################################
                val = np.array(val)
                col_times = np.array(col_times)
                
            elif last_col in cols:
                #group8
                col=cols[0]
                val = []
                col_times=[]
                ########################################
                t = [times[col]]
                idx = np.array([col-1,col+1])
                pts_ = pts.iloc[idx]
                val.append(self.linear(pts_,t)[0])
                col_times.append(t[0])

                ########################################
                t = [times[last_col]]
                idx = np.array([last_col-2,last_col-1])
                pts_ = pts.iloc[idx]
                val.append(self.linear(pts_,t)[0])
                col_times.append(t[0])
                val = np.array(val)
                col_times = np.array(col_times)
            else:
                #group9
                val = []
                col_times=[]
                for col in cols:
                    t = [times[col]]
                    idx = np.array([col-1,col+1])
                    pts_ = pts.iloc[idx]
                    val.append(self.linear(pts_,t)[0])
                    col_times.append(t[0])
                col_times = np.array(col_times)
                
        return col_times,np.array(val)
    

        
if __name__ == "__main__":
    # filename = get_file()
    # df = pd.read_csv(filename)
    # imp = Imputer(df,2)

    pass            

