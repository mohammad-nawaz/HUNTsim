#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 15:13:35 2022

@author: nawaz
"""

import numpy as np
import matplotlib.pyplot as plt
import dataprep as dp
import pandas as pd
from huntsim import ReadConfig
# from displayvals import FollowDotCursor

class Plotter:
    """
    """
    
    def plotDuo(self, sname, var1, var2, midday_data = None):
        var1 = var1.upper(); var2 = var2.upper()
        ddc = dp.DayDataCheck()
        df = ddc.dayDataCheck(sname, midday_data =  midday_data)

        x = np.array(df['Date'].tolist())
        y1 = np.array(df[var1].tolist())
        y2 = np.array(df[var2].tolist())
        
        # Create figure and axis objects 
        fig, ax = plt.subplots()
        fig.autofmt_xdate(rotation=90)  # rotate the x-labels
        
        ax.plot(x, y1, color='red', marker='o' )
        
        # set x- and y-axis labels
        ax.set_xlabel('Dates')
        ax.set_ylabel(var1)
        
        # title
        ax.set_title(sname)
        
        ax2=ax.twinx()
        ax2.plot(x, y2, color='blue',marker='o')
        #set y-axis label
        ax2.set_ylabel(var2)
        
        # cursor1 = FollowDotCursor(ax, x, y1)
        # cursor2 = FollowDotCursor(ax2, x, y2)
        
        plt.show()
        
    def plotDuo_n(self, sname, var1=None, var2=None, midday_data = None):
        # var1 = var1.upper(); var2 = var2.upper()
        ddc = dp.DayDataCheck()
        df = ddc.dayDataCheck(sname, midday_data =  midday_data)
        
        # Read configuration file
        rconf = ReadConfig()
        rconf.readConfig()

        x = df['Date'].tolist()
        # y1 = df['COTP'].tolist()
        y1 = df['CP'].tolist()
        y2 = df['HB/HS'].tolist()
        mbs_ratio = df['MB/MS'].tolist()
        rbv = df['RBV'].tolist()
        rbvd = df['RBVD'].tolist()
        rbv_diff = rbvd[-1]
        
        # Create figure and axis objects 
        fig, ax = plt.subplots(2)
        fig.autofmt_xdate(rotation=90)  # rotate the x-labels
        
        ax[0].plot(x, y1, color='black', marker='o' )
        
        # set x- and y-axis labels
        ax[0].set_xlabel('Dates')
        ax[0].set_ylabel('cotp')
        
        # title
        ax[0].set_title(sname)
        
        ax2=ax[0].twinx()
        ax2.plot(x, y2, color='blue',marker='o')
        #set y-axis label
        ax2.set_ylabel('hb/hs')
        
        # Display buy flag if the condition is met 
        if rconf.Buy_params_count == 1:
            if y2[-1]>=rconf.b_hbs_ratio:
                ax[0].text(x[-1],y1[-1], 'B', color = 'green', fontsize=16)
            
        if rconf.Buy_params_count == 2:
            if y2[-1]>=rconf.b_hbs_ratio and mbs_ratio[-1]>=rconf.b_mbs_ratio:
                ax[0].text(x[-1],y1[-1], 'B', color = 'green', fontsize=16)            
            
        if rconf.Buy_params_count == 3: 
            if y2[-1]>=rconf.b_hbs_ratio and mbs_ratio[-1]>=rconf.b_mbs_ratio and rbv_diff>rconf.b_rbv_diff:
                ax[0].text(x[-1],y1[-1], 'B', color = 'green', fontsize=16) 
                
        print ('Latest price:', y1[-1])
        print('Latest cotp:', df['COTP'].tolist()[-1])
        print('Latest hb/hs:', y2[-1]) 
        print('Latest mb/ms:', mbs_ratio[-1])
        
        # Sell 
        # set x- and y-axis labels
        y2 = df['RBV'].tolist()
        ax[1].set_xlabel('Dates')
        ax[1].set_ylabel('cotp')
        
        ax[1].plot(x, y1, color='black', marker='o' )
        
        # set x- and y-axis labels
        ax[1].set_xlabel('Dates')
        ax[1].set_ylabel('cotp')
        
        # title
        # ax[1].set_title(sname)
        
        ax2=ax[1].twinx()
        ax2.plot(x, y2, color='blue',marker='o')
        #set y-axis label
        ax2.set_ylabel('rbv')
 
        # Display sell flag if the condition is met 
        if (rbv_diff<=rconf.s_rbv_diff):
            ax[1].text(x[-1],y1[-1], 'S', color = 'r', fontsize=16)
        
        print('Latest rbv:', y2[-1]) 
        print('Latest rbv_diff:', rbv_diff)   
        
        plt.show()   
        
        
    def plotBarMvp(self, sname, date):
        # Read data from file: 
        sname = sname.upper()+'.csv'
        dfmin = pd.read_csv(sname)
            
        dfmin = dfmin[dfmin['Date']==date]
        
        t = dfmin["Time"].tolist()
        value = dfmin["Value"].tolist()
        
        fig = plt.figure(figsize = (10,7))
        
        plt.bar(t, value)
        plt.xticks(rotation=90, ha='right')

        plt.show()
        

            
        

