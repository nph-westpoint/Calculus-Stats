# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 09:24:32 2023

Vectorized on 24 March 2026

@author: grover.laporte
"""

import pandas as pd
import numpy as np

import streamlit as st
import glucose_splines as gs

if 'curve' not in st.session_state:
    st.session_state['curve'] = None

if 'current_file' not in st.session_state:
    st.session_state['current_file'] = None

if 'idx' not in st.session_state:
    st.session_state['idx'] = 0


current_file = st.session_state['current_file']
curve = st.session_state['curve']
idx = st.session_state['idx']

if current_file is None:
    st.write("""
          # Calculus metrics from Mixed Meal Tolerance Test
          ## Kondonculator 3.0 
          Input your csv file with times (top row needs to be the times that the data was collected)
          as columns and observations as rows. Ensure you have the correct number of non-numerical 
          columns.
          """)
    col = st.sidebar.number_input(":red[Number of non-numerical columns:]",
                      min_value=0,max_value=20,value=2,step=1,
                help="The number of columns used to identify the observation")
    current_file = st.file_uploader("Choose a file")
    st.session_state['current_file'] = current_file
    

if current_file is not None and curve is None:
    df = pd.read_csv(current_file)
    try:
        curve = gs.Curve(df,col)
        st.session_state['curve'] = curve

    except:
        res = 'Ensure your file is arranged with observations in the rows, '
        res+='times along the columns, the number of non-numerical columns '
        res+='identified, and the column names of the numerical columns as '
        res+='numbers.'
        st.write(res)
    

if curve is not None:
    curve_options = list(curve.data.index)
    idx_name = st.sidebar.selectbox(
        "Select observation:",
        options = curve_options,
        index=0
    )
    idx = curve_options.index(idx_name)
    st.header(idx_name)
    options = [
        "Load Data",
        "Visual - Calculus Metrics",
        "Visual - Derivatives",
        "Download Stats",
    ]
    select = st.sidebar.radio(
        label = "Select a View:",
        options = options,
        index = 0,
        key = 'view_key'
    )
    if select == options[0]:
        st.write("The file has been loaded.")
        st.write(curve.stats)

    if select == options[1]:
        st.plotly_chart(curve.plot_curve(idx,['f','d']))
        d1,d2 = curve.display(idx)
        d1.columns = ['stats']
        d2.columns = ['critical points']
        cols = st.columns([1,2])
        with cols[0]:
            st.write(d1.style.format(precision=1))
        with cols[1]:
            st.write(d2.style.format(precision=1))

    if select == options[2]:
        check_opts = [
            'cubic spline',
            '1st Derivative',
            'Second Derivative',
            'Data',
        ]
        cols = st.columns(4)
        charts = []
        with cols[0]:
            fig1 = st.checkbox("Cubic Spline:",
                        key = 'chk_cs',
                        )
            if fig1:
                charts.append('f')
        with cols[1]:
            fig2 =st.checkbox("Data",
                        key = 'chk_data',
                        )
            if fig2:
                charts.append('d')
        with cols[2]:
            fig3 = st.checkbox("Derivative",
                        key = 'chk_fp',
                        )
            if fig3:
                charts.append('fp')
        with cols[3]:
            fig4 = st.checkbox("Second Derivative",
                        key = 'chk_fpp',
                        )
            if fig4:
                charts.append('fpp')
        st.plotly_chart(curve.plot_curve(idx,charts))

    if select == options[3]:
        st.download_button("Download CSV",
                           data = curve.stats.to_csv().encode("utf-8"),
                           file_name="stats.csv",
                           mime='text/csv',
                           icon = ":material/download:"
                           )
        st.write(curve.stats)

        


    








        






