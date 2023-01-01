#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 11:35:17 2022

@author: nawaz
"""

import matplotlib.pyplot as plt
import dataprep as dp
import pandas as pd

with open('paths.txt') as f:
    paths = [line.rstrip() for line in f]
    
spath = paths[2]

sname = 'kohinoor'

df = pd.read_csv(spath+sname+'interpolated.csv')

dt = 'dec12'
dfd = df[df['Date']==dt]

dfdb = dfd[dfd["Type"]=="Buy"]
dfds = dfd[dfd["Type"]=="Sell"]

fig = plt.figure(1, figsize=(8,10))

lmargin = rmargin = bmargin = tmargin = 0.1; pad = 0.01

left0 = lmargin
bottom0 = bmargin
width0 = 0.35
height0 = 0.8

left1 = lmargin+width0+pad
bottom1 = bottom0
width1 = width0
height1 = height0

ax0 = plt.axes([left0, bottom0, width0, height0])
ax1 = plt.axes([left1, bottom1, width1, height1])

# Calculate COTP
price_list = dfd["Price"].tolist()
value_list = dfd["Value"].tolist()
pval = [i*j for i, j in zip(price_list,value_list)]
cotp = round(sum(pval)/sum(value_list),1)
cotp = str(cotp)
                
# Price list
price = dfd["Price"].tolist()
pb = dfdb["Price"].tolist()
tb = dfdb["Type"].tolist()
ps = dfds["Price"].tolist()
ts = dfds["Type"].tolist()


# Value list
vb = dfdb["Value"].tolist()
vs = dfds["Value"].tolist()

lwmax = 10
maxv = 50
for i,j in zip(pb,vb):
    lw = lwmax*(j-min(vb))/(maxv-min(vb))
    alpha = (j-min(vb))/(maxv-min(vb))
    ax0.axhline(i, linewidth=lw, alpha=alpha, color='g')
ax0.axhline(cotp, linewidth=3, color='k')
ax0.set_title(sname +'  '+ dt)
    
# Set y-lim 
ymin = min(price)#-0.1*min(price)
ymax = max(price)#+0.1*max(price)
ax0.set_ylim(ymin, ymax)

for i,j in zip(ps,vs):
    lw = lwmax*(j-min(vs))/(maxv-min(vs))
    alpha = (j-min(vs))/(maxv-min(vs))
    ax1.axhline(i, linewidth=lw, alpha=alpha, color='r')
ax1.axhline(cotp, linewidth=3, color='k')
    
# Set y-lim 
ymin = min(price)#-0.1*min(price)
ymax = max(price)#+0.1*max(price)
ax1.set_ylim(ymin, ymax)