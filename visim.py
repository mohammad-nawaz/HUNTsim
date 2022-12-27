# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 12:33:36 2022

@author: Reza
"""

import os
import ast
import pandas as pd
import tcalendar
import configparser
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
    def dataVis(self, name, storepng= None, start_date= None):
        name = name.upper()
        if storepng == 'on':
            with open('paths.txt') as f:
                paths = [line.rstrip() for line in f]
    
            parent_dir = paths[5]
            destination = parent_dir + name
            if not os.path.exists(destination):
                os.mkdir(destination)
            else:
                if os.path.exists(destination):
                    shutil.rmtree(destination)
                os.mkdir(destination)
            destination += '/'
            
            ddc = dataprep.DayDataCheck()
            df = ddc.dayDataCheck(name)
        else:
            ddc = dataprep.DayDataCheck()
            df = ddc.dayDataCheck(name, midday_data = 'on')
        sizedata = len(df)
        
        
        high = int(1.1*max(df['HP'])); low = int(0.9*min(df['LP']))
        tdates = df['Date'].tolist()
        floorP = pd.read_csv('floorPrice.csv')
        fprice = floorP[floorP['Name']==name]['Floor'].item()
        COT = []; Val = []; cotVal = []; buy_sig =0; sell_sig = 0; bigTr = []
        tot_pair = []; len_cotVal = []; i =0; count = -1; flag = True; highflag = True; buysell = 0
        for td in tdates:
            df_td = df[df['Date']==td]
            count += 1
            cotp = df_td['COTP'].item()
            cp = df_td['CP'].item()
            if np.isnan(df_td['CP'].item()):
                continue
            
            ##Calculating sig###
            val = df_td['TVAL'].item()
            COT.append(cotp)
            Val.append(val)
            tot_pair.append([cotp,val])
            if i>1:
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
            pRatio = cotp/fprice;
            
            if i>1:
                if len_cotVal[-1]-len_cotVal[-2] <= -4:
                        if pRatio >=1.1:
                            buy_sig = 1
                            flag = True
                            cat = 1
                            sell_sig = 0
                            bigTr.append(cotp)
                            if not len(bigTr) == 0:
                                if bigTr[-1] >= 1.4*bigTr[0]:
                                    buy_sig =0
                                    bigTr = []
                                    highflag = False
                            
                elif i>2 and len_cotVal[-1]-len_cotVal[-3]<=-6:
                    # if pRatio >=1.1:
                    buy_sig = 1
                    flag = True
                    cat = 2
                    sell_sig = 0
                elif i>3 and len_cotVal[-1]-len_cotVal[-4]<=-7:
                    # if pRatio >=1.1:
                    buy_sig = 1
                    flag = True
                    cat =3
                    sell_sig = 0
                elif i>1 and len_cotVal[-1]==0 and len_cotVal[-2]==0:
                    if 1.8 >= pRatio >= 1.2:
                        buy_sig = 1
                        cat = 0
                        sell_sig = 0
                        flag = True
                    #and len_cotVal[-3]==0
                else:
                    buy_sig = 0
                # if len_cotVal[-1] - len_cotVal[-2] >= 1:
                #     sell_sig = 1
                #incorporating sell signal into buy signal
                #to avoid creating buy signal on the sell signal
                rconf = ReadConfig()
                rconf.readConfig()
                if rconf.Sell_params_count == 1: 
                    if df_td['RBVD'].item()>=rconf.s_rbv_diff:
                        sell_sig = 0
                    if df_td['RBVD'].item()<rconf.s_rbv_diff:
                        buy_sig = 0
                        sell_sig = 1
                        flag = False
            i += 1
        
            if storepng == 'on':
                idx = 1
                if not start_date == None:
                    idx = tdates.index(start_date)
            else:
                idx = sizedata - 4
            if count>= idx:
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
                
                ax1.annotate(f'current:{len_cotVal[-1]}, prev:{len_cotVal[-2]}', (-0.05, high*0.97))
                k = 1
                for j in range(len(cotVal)):
                    k = -1*k
                    # ax1.scatter(0, cotVal[j][0], c='yellow',marker ="o", alpha=1, edgecolor = 'black', s= cotVal[j][1]//1000)
                    if k>0:
                        ax1.annotate(f'{int(cotVal[j][1])}', (-0.01-j*0.001,cotVal[j][0]+0.05))
                    if k<0:
                        ax1.annotate(f'{int(cotVal[j][1])}', (0.001+j*0.003,cotVal[j][0]+0.05))
                ax1.axhline(y= cotp, xmin=0, xmax=1, c="black", linewidth=2, zorder=0)
                ax1.annotate(f'{int(Val[i-2])}', (-0.05,COT[i-2]+0.05), color = 'blue')
                # ax1.annotate(f'{cotp}', (-0.05,cotp*1.01))
                ###########################--------------ax2-----------###############################
                ax2 = plt.subplot(gs1[1])
                ax2.set_ylim([low,high])  
                ax2.get_yaxis().set_visible(False)
                ax2.get_xaxis().set_visible(False)
                ax2.axhline(y= cotp, xmin=0, xmax=1, c="black", linewidth=2, zorder=0)
                ax2.annotate(f'{cotp}', (0.04,cotp*1.01))
                ax2.annotate(f'(mat_val:{int(Val[i-3])})', (0.5,COT[i-3]+0.05), color = 'green')
                ax2.annotate(f'(cur_val: {int(Val[-1])})', (0.005,cotp*0.983), color = 'blue')
                
                if buy_sig == 1:
                    if cat == 0 and flag:
                        ax2.annotate('BUY_SIGNAL: Continuous high\ntransactions without gap', (0.005,high*0.9))
                    if cat == 1 and flag:
                        ax2.annotate('BUY_SIGNAL: High transaction\ngobbled at least\n4 of the\nprevious levels of stocks', (0.005,high*0.88))
                    if cat == 2 and flag:
                        ax2.annotate('BUY_SIGNAL: High transaction\ncontinued for two days', (0.005,high*0.9))
                    if cat == 3 and flag:
                        ax2.annotate('BUY_SIGNAL: High transaction\ncontinued for three days', (0.005,high*0.9))
                    if sell_sig == 1 and flag == False:
                        ax2.annotate('BUY_SIGNAL couldnot generate\nbecause RBVD is\nbelow threshold', (0.005,high*0.9))
                        buysell =1
                        flag = True
                    if highflag == False:
                        ax2.annotate('BUY_SIGNAL couldnot generate\nbecause price is\nsignificantly high', (0.005,high*0.9))
                        highflag = True
                if sell_sig == 1 and buysell == 0:
                    ax2.annotate('Sell Signal', (0.005,high*0.95))
                
                buysell = 0
                    # ax2.annotate('BUY_SIGNAL', (0.005,high*0.9))
                if storepng == 'on':
                    fignum = count
                    if fignum//10==0:
                        outputname = destination + f'00{fignum}.png'
                    else:
                        outputname = destination + f'0{fignum}.png'
                    
                    fig.savefig(outputname,bbox_inches='tight')
                    plt.close()
        
        print('date:',td)
        print('RBVD:',df_td['RBVD'].item())
vs = Visual()
name = 'SEAPEARL'
vs.dataVis(name)