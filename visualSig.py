# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 03:33:22 2022

@author: Reza
"""

import os
import ast
import pandas as pd
import dataprep
import seaborn as sns
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib import gridspec

name = 'GENEXIL'

with open('paths.txt') as f:
    paths = [line.rstrip() for line in f]

parent_dir = paths[3]
road = parent_dir + name
if not os.path.exists(road):
    os.mkdir(road)
else:
    if os.path.exists(road):
        os.rmdir(road)
    os.mkdir(road)
road += '/'

ddc = dataprep.DayDataCheck()
df = ddc.dayDataCheck(name)
high = int(1.1*max(df['HP'])); low = int(0.9*min(df['LP']))
tdates = df['Date'].tolist()
floorP = pd.read_csv('floorPrice.csv')
fprice = floorP[floorP['Name']==name]['Floor'].item()
COT = []; Val = []; y = []; buy_sig =0; sell_sig = 0;
tot_pair = []; len_y = []; i =0
for td in tdates:
    df_td = df[df['Date']==td]
    
    cotp = df_td['COTP'].item()
    cp = df_td['CP'].item()
    if np.isnan(df_td['CP'].item()):
        continue
    
    ##Calculating sig###
    val = df_td['TVAL'].item()
    # strVal = str(val)
    # if not strVal == 'nan':
    COT.append(cotp)
    Val.append(val)
    tot_pair.append([cotp,val])
    if i>1:
        y.append(tot_pair[-3])
        y.sort()
        while val>0 and len(y)>0:
            if val> y[0][1]:
                val -= y[0][1]
                y.pop(0)
            elif val<= y[0][1]:
                y[0][1] -= val
                val = 0
    len_y.append(len(y))
    pRatio = cotp/fprice;
    if i>1:
        if len_y[-1]-len_y[-2] <= -3:
                if pRatio >=1.1:
                    buy_sig = 1
                    sell_sig = 0
                    
        else:
            buy_sig = 0
        if len_y[-1] - len_y[-2] >= 1:
            sell_sig = 1
    
    # if td == 'nov22':
    #     break
    # print('y:', y)
    # print('\nLength of y:', len_y)
    # if i>1:
    #     print('diff:', len_y[-1]-len_y[-2])
    # if buy_sig ==1:
    #     print('Buy_Signal')
    # if sell_sig == 1:
    #     print('Sell_Signal')
    # print('Cot:',cotp)
    # print('---------------')
    
    fig = plt.figure(facecolor='cyan')

    fig.set_figheight(5)
    fig.set_figwidth(5)
    fig.suptitle(f'Date: {td} \n {name}', fontsize=10)
    gs1 = gridspec.GridSpec(1, 2)
    gs1.update(wspace=0.01, hspace=0.1)
    
    ###################--------------ax1-------------###############
    ax1 = plt.subplot(gs1[0])
    ax1.margins(0.05)          
    ax1.set_ylim([low,high])
    ax1.get_xaxis().set_visible(False)
    k = 1
    if i>1:
        ax1.annotate(f'current:{len_y[-1]}, prev:{len_y[-2]}', (-0.05, high*0.97))
    for j in range(len(y)):
        k = -1*k
        ax1.scatter(0, y[j][0], c='yellow',marker ="o", alpha=1, edgecolor = 'black', s= y[j][1]//4)
        if k>0:
            ax1.annotate(f'{int(y[j][1])}', (-0.02,y[j][0]+0.05))
        if k<0:
            ax1.annotate(f'{int(y[j][1])}', (0.005,y[j][0]+0.05))
    ax1.axhline(y= cotp, xmin=0, xmax=1, c="black", linewidth=2, zorder=0)
    ax1.annotate(f'{cotp}', (-0.05,cotp*1.01))
    ###########################--------------ax2-----------###############################
    ax2 = plt.subplot(gs1[1])
    ax2.set_ylim([low,high])  
    ax2.get_yaxis().set_visible(False)
    ax2.get_xaxis().set_visible(False)
    ax2.axhline(y= cotp, xmin=0, xmax=1, c="black", linewidth=2, zorder=0)
    if buy_sig == 1:
        ax2.annotate(f'BUY_SIGNAL', (0.05,high*0.9))
    
    fignum = i
    if fignum//10==0:
        outputname = road + f'00{fignum}.png'   
    else:
        outputname = road + f'0{fignum}.png'
    
    fig.savefig(outputname,bbox_inches='tight')
    plt.close()    
    i += 1
# # plt.show()
# # if i> last-100:
# #     plt.show()
# # else:
#     # plt.close()


    
# fig.savefig(outputname,bbox_inches='tight')
# plt.close()
# i += span
# fignum += 1


    
    
        
    # if i//10==0:
    #     outputname = store + f'00{i}.png'   
    # else:
    #     outputname = store + f'0{i}.png'
        
    # fig.savefig(outputname,bbox_inches='tight')

















# y = [m/n if n else 0 for m,n in zip(df['Volume'], df['TimeDiff'])]
# fig, ax1 = plt.subplots(figsize=(100,5))
# sns.lineplot(data = df['Price'],color = 'green', marker='o', sort = False, ax=ax1)
# ax2 = ax1.twinx()

# ax2.stem(y, linefmt = 'blue')

# filename = store_path + name+ '.csv'
# df = pd.read_csv(filename)
# store_path = 'E:/Excel_Stock/Archive_mindata/'
# name = 'Ehl'
# name = name.upper()
# filename = store_path + name + '.csv'
def recordDates(name): 
    with open('recordDates.txt') as f:
        data = f.read()
    rD = ast.literal_eval(data)
    f.close()
    if name in rD.keys():
        recordDate = rD[name]
    else: recordDate = []
    return recordDate
    

