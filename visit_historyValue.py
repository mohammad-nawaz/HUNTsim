# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 22:01:19 2022

@author: Reza
"""
import os
import ast
import pandas as pd
import tcalendar
import numpy as np
from readconfig import ReadConfig
import dataprep
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import shutil

class Visual:
    '''
    '''
    def dataVis(self, name, start_date, end_date, storepng= None):
        name = name.upper()
        storepng = 'on'
        if storepng == 'on':
            with open('paths.txt') as f:
                paths = [line.rstrip() for line in f]
    
            parent_dir = paths[5]
            history_store = paths[6]
            destination = parent_dir + name + '_history'
            if not os.path.exists(destination):
                os.mkdir(destination)
            else:
                if os.path.exists(destination):
                    shutil.rmtree(destination)
                os.mkdir(destination)
            destination += '/'
            
            df = pd.read_csv(history_store+name+'.csv')
            df = df.dropna(how= 'all')
        else:
            df = pd.read_csv(history_store+name+'.csv')
            df = df.dropna(how= 'all')
        sizedata = len(df)
        
        
        high = int(1.1*max(df['High'])); low = int(0.9*min(df['Low']))
        tdates = df['Date'].tolist()
        idx1 = tdates.index(start_date)
        idx2 = tdates.index(end_date)
        COT = []; Val = []; cotVal = []; #buy_sig =0; sell_sig = 0; bigTr = []
        tot_pair = []; len_cotVal = []; i =0; count = -1; #flag = True; highflag = True; buysell = 0
        for i in range(idx1-1, idx2):
            count += 1
            rawval = df['Value(mn)'][i]
            if type(rawval) == str:
                rawval = rawval.replace(",","")
            val = float(rawval)*10
            cotp = round(df['Avg (VWAP)'][i],2)
            cp = round(df['Close'][i],1)
            COT.append(cotp)
            Val.append(val)
            tot_pair.append([cotp,val])
            if i>idx1+1:
                cotVal.append(tot_pair[-3])
                cotVal.sort()
                while val>0 and len(cotVal)>0:
                    if val> cotVal[0][1]:
                        val -= cotVal[0][1]
                        cotVal.pop(0)
                    elif val<= cotVal[0][1]:
                        cotVal[0][1] -= val
                        val = 0
            len_cotVal.append(len(cotVal))
            i += 1
        
            if count> 1:
                fig = plt.figure(facecolor='cyan')
                # print('val:', Val)
                # print('count', count)
                fig.set_figheight(5)
                fig.set_figwidth(5)
                fig.suptitle(f'Date: {tdates[i]} \n {name}', fontsize=10)
                gs1 = gridspec.GridSpec(1, 1)
                gs1.update(wspace=0.01, hspace=0.1)
                
                ###################--------------ax1-------------###############
                ax1 = plt.subplot(gs1[0])
                ax1.margins(0.05)          
                ax1.set_ylim([low,high])
                ax1.get_xaxis().set_visible(False)
                
                # ax1.annotate(f'mat_val:{int(Val[count-2])}, cur_val: {int(Val[-1])}', (-0.05, high*0.97))
                k = 1
                for j in range(len(cotVal)):
                    k = -1*k
                    # ax1.scatter(0, cotVal[j][0], c='yellow',marker ="o", alpha=1, edgecolor = 'black', s= cotVal[j][1]//1000)
                    if k>0:
                        ax1.annotate(f'{int(cotVal[j][1])}', (-0.01-j*0.001,cotVal[j][0]+0.05))
                    if k<0:
                        ax1.annotate(f'{int(cotVal[j][1])}', (0.001+j*0.003,cotVal[j][0]+0.05))
                ax1.axhline(y= cotp, xmin=0, xmax=1, c="black", linewidth=2, zorder=0)
                ax1.annotate(f'{int(Val[count-1])}', (-0.05,COT[count-1]+0.05), color = 'blue')
                ax1.annotate(f'cot:{cotp}', (0.03,cotp*1.01))
                ax1.annotate(f'(mat_val:{int(Val[count-2])})', (0.025,COT[count-2]), color = 'green')
                ax1.annotate(f'(cur_val: {int(Val[-1])})', (-0.05,cotp*0.94), color = 'blue')
                
                if storepng == 'on':
                    fignum = count
                    if fignum//10==0:
                        outputname = destination + f'00{fignum}.png'
                    else:
                        outputname = destination + f'0{fignum}.png'
                    
                    fig.savefig(outputname,bbox_inches='tight')
                    plt.close()
vs = Visual()
name = 'BXPHARMA'
start_date = '25-Jun-20'
end_date = '7-Nov-21'
vs.dataVis(name, start_date, end_date)