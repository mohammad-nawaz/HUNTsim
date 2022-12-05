#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 22:28:19 2022

@author: nawaz
"""

import pandas as pd

def interpolation(df, dfn, idx):
    # create a temporary dataframe for the neutral data
    dftemp = pd.DataFrame(columns=df.columns)
    pstep0 = df["Price"][idx-1]-df["Price"][idx-2]
    vstep0 = df["Vstep"][idx-1]
    dftemp = pd.DataFrame(columns=df.columns)
    # create a temporary dataframe to store a copy of dftemp
    dfs = pd.DataFrame(columns=df.columns)
    while idx<len(df):
        pstep = df["Price"][idx]-df["Price"][idx-1]
        vstep1 = df["Vstep"][idx]
        # Find the block of neutral data
        if pstep==0:
            dftemp_1 = dftemp.append(df[idx:idx+1])
            dftemp = dftemp_1
            dftemp = dftemp.reset_index(drop=True)
            idx += 1
        # We change the neutral data with the interpolated data as follows:
        #       If Vsteps before and after the neutral data is the same, say
        #       Buy, the neutral data type changes to Buy
        #       If, Vstep_before and Vstep_after are not the same, we calculate 
        #       weights of buy and sell Volumes and Values from the Vsteps. 
        #       Here we use a scale of 1 for Vstep_before-Vstep_after.
        #       Example: Vstep_before = 2; Vstep_after = -1 splits a data as
        #       follows:
        #           scale: abs(2) + abs(-1) --> 1
        #           Vol1 = 2/(abs(2)+abs(-1) * Vol
        #           Vol2 = 1/(abs(2)+abs(-1) * Vol

        if pstep!=0:
            pstep1 = pstep
            bfrac = abs(max(vstep0,vstep1)/(abs(vstep0)+abs(vstep1)))
            sfrac = abs(min(vstep0,vstep1)/(abs(vstep0)+abs(vstep1)))
            for row in dftemp.index:
                if pstep0>0 and pstep1>0:
                    dftemp.at[row, "Type"]="Buy"
                elif pstep0<0 and pstep1<0:
                    dftemp.at[row, "Type"]="Sell"
                # elif pstep0>0 and pstep1<0: 
                else:
                    df2s = pd.DataFrame(dftemp.loc[[row]])  
                    
                    dftemp.at[row,"Type"]="Buy"
                    dftemp.at[row,"Volume"]=int(bfrac*dftemp.at[row,"Volume"])
                    dftemp.at[row,"Value"]=bfrac*dftemp.at[row,"Value"]
                    
                    df2s.at[row,"Type"]="Sell"
                    df2s.at[row,"Volume"]=int(sfrac*df2s.at[row,"Volume"])
                    df2s.at[row,"Value"]=sfrac*df2s.at[row,"Value"]
                    
                    dfs=pd.concat([dfs,df2s])

            dfn = pd.concat([dfn,dftemp, dfs])
            dfn = dfn.reset_index(drop=True)
            break
     
    return dfn, idx

def interpolatedData(df):
    # create an empty dataframe for the interpolated data
    dfn = pd.DataFrame(columns=df.columns)
    
    # Find and concatenate the morning data
    idx = -1
    while idx<len(df):
        idx+=1
        if df["Vstep"][idx]!=0:
            break
    mdf = pd.DataFrame(df.loc[0:idx])
    dfn = pd.concat([dfn,mdf])   
    
    # find and concatenate the day neutral data
    while idx<len(df):
        price0 = df["Price"][idx-1]
        price1 = df["Price"][idx]
        pstep = price1-price0
        if pstep!=0:
            dfn_ = dfn.append(df[idx:idx+1])
            dfn=dfn_
            dfn=dfn.reset_index(drop=True)
            idx += 1
        if pstep==0:
            dfn, idx = interpolation(df,dfn, idx)
            lidx = idx
            
    # afternoon data
    last_type = dfn['Type'].iat[-1]
    idx = len(df)-1
    print(idx)
    while idx<(len(df)):
        idx = idx - 1
        print (df["Vstep"][idx])
        if df["Vstep"][idx]!=0:
            break
        
    edf = pd.DataFrame(df.loc[idx+1:len(df)])
    for row in edf.index:
        edf.at[row, "Type"]=last_type
    dfn = pd.concat([dfn,edf])
    dfn = dfn.reset_index(drop=True)
    
    return dfn

# read datafile
df = pd.read_csv('test.csv')
dfn = interpolatedData(df)
        
dfn.to_csv('interp.csv')        



