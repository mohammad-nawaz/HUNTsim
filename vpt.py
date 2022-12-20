#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 23:21:19 2022

@author: nawaz
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import dataprep as dp
import os
import math
import matplotlib
import tcalendar as tc


# plt.ioff()
plt.ion()

class MomentumSignal:

    def distFunction(self, lis):
        pass
    
    def sumListItems(self, num_entry, lis, idx):
        ''' Return sum of items of a list '''
        if idx<=num_entry:
            return 0
        else: 
            return sum(lis[idx-num_entry:idx])
        
    def buyFrac(self, num_entry, lis, typ, cut_val=None):
        lsval = lis.copy()
        lstyp = typ.copy()
        newval=[]
        for i,j in zip(lsval, lstyp):
            if j=='Buy':
                newval.append(i)
            else:
                newval.append(-i)
        
        buyfrac=[]
        buyfract=[]
        count =0
        while count<len(lsval):
            if count<num_entry:
                buyfrac.append(0)
                buyfract.append(0)
            if count>=num_entry:
                cutval=cut_val
                cur_window = newval[count-num_entry:count]
                bval = sum(x for x in cur_window if x>cutval)
                sval = abs(sum(x for x in cur_window if x<-cutval)) 
                count_high = sum(i >= cutval for i in cur_window)
                sval_t = abs(sum(x for x in cur_window if x<0)) 
                bval_t = bval
                if count_high==1:
                    sval = sval_t
                if (bval+sval)==0:
                    bval=0.5; sval=0.5 # so that the denominator becomes 1 with equal bval & sval
                if (bval_t+sval_t)==0:
                    bval_t = 0.5
                    sval_t = 0.5
                bfrac = bval/(bval+sval)
                bfract = bval_t/(bval_t+sval_t)
                buyfrac.append(bfrac)
                buyfract.append(bfract)
            count+=1
                
        return buyfrac, buyfract
    
    def getListIndices(self, ls):
        '''Return two list of lists zs=[[]] and nzs=[[]]
        zs contains list of indices for islands of zeros 
        nzs contains list of indices for islands of nonzeros
        '''
        # make a copy of the list
        lis = ls.copy()
        
        # Define lists and lists of lists    
        z=[]; zs=[]
        nz=[]; nzs=[]
        # Define a flag to identify zeros and nonzeros in a loop
        if lis[0]==0:
            pflag=0
        else: 
            pflag=1
        
        # Define a counter
        count = 0
        while len(lis)!=0:
            if lis[0]==0:
                flag = 0
                z.append(count)
            if lis[0]!=0:
                flag = 1
                nz.append(count)
            if (pflag-flag)==-1:
                zs.append(z)
                z=[]
            if (pflag-flag)==1:
                nzs.append(nz)
                nz=[]
            pflag = flag
            # pop the first item of the list before reentering into the loop
            lis.pop(0)
            count +=1 
        
        # add the last unappended list if any
        if len(z)!=0:
            zs.append(z)
        if len(nz)!=0:
            nzs.append(nz)
            
        return zs, nzs
    
    def modelList(self, lsprice, lstype):
        ''' Model the null psteps with the neighbouring psteps'''
        zs, nzs = self.getListIndices(lsprice)
        # modelled price step
        nls = lsprice.copy()
        if len(nzs)==0: 
            pass
        for i in zs: 
            if i.count(0)>0 and len(nzs)!=0:
                idx = i[-1]+1
                for j in i:
                    nls[j]=lsprice[idx]                               
            if i.count(len(lsprice)-1)>0:
                idx = i[0]-1
                for j in i: 
                    nls[j]=lsprice[idx]
            if i.count(0)==0 and i.count(len(lsprice)-1)==0: 
                idx1=i[0]-1
                idx2=i[-1]+1
                bpstep = max(lsprice[idx1], lsprice[idx2])
                spstep = min(lsprice[idx1], lsprice[idx2])
                if lsprice[idx1]*lsprice[idx2]>0:
                    bpstep = spstep = (lsprice[idx1]+lsprice[idx2])/2
                for j in i: 
                    if lstype[j]=='Buy':
                        nls[j]=bpstep
                    if lstype[j]=='Sell':
                        nls[j]=spstep
        return nls
    
    def momentumSignal(self, sname=None, sdate=None, edate=None, midday_data='on'):    
        from datetime import date, time
        
        with open('paths.txt') as f:
            paths = [line.rstrip() for line in f]
            
        spath = paths[2]    # stock path
        mspath = paths[3]   # path for saving figures
        stpath = paths[4]    # path for the file StockList.txt
        
        # Reading name of stocks 
        with open(stpath+'StockList.txt') as f: 
            dfile = [line.rstrip() for line in f]
        
        # dfile = dfile[0:1].copy()
        
        if midday_data == 'on':
            dfile = [sname]
        
        font = {'family' : 'normal',
                'weight' : 'bold',
                'size'   : 12}
        
        markersize = 100
        
        num_entry = 50
        
        midday_data = 'on'
        
        
        
        for i in dfile: 
            sname = i
            sname = sname.upper()
            
            df = pd.read_csv(spath+sname+'interpolated.csv')
            
            # midday data prep
            md = dp.MinData()
            df_mid = md.minMiddayDataPrep(name = sname)
            df = pd.concat([df, df_mid])
            df = df.reset_index(drop=True)
            
            months=['jan','feb','mar','apr','may','jun',
                    'jul','aug','sep','oct','nov','dec']
            
            tod = date.today()
            m = months[tod.month-1]
            d = tod.day
            if d//10==0:
                d = '0'+str(d)
            else: 
                d = str(d)
                
            if edate==None: 
                end_date = m+d
            else: 
                end_date = edate
                
            if sdate==None: 
                start_date = end_date
            else: 
                start_date = sdate
            
            # Make list of tdates 
            tdates = tc.TradingDates()
            td = tdates.tradingDates(i, start_date, end_date)
            # td = []
            # [td.append(x) for x in df['Date'].tolist() if x not in td]
            
            # Plot Val per trades
            # destination of figures
            figpath = os.path.join(mspath,sname)
            if not os.path.isdir(figpath):
                os.mkdir(figpath)
            figpath = figpath+'/'
            
            
            count = 0
            for date in td:
                print(sname, date)
                # dataframe of a day
                dfn=df[df['Date']==date]
                
                # Model pstep and vstep
                psteps = dfn["Pstep"].tolist()
                types = dfn["Type"].tolist()
                values=dfn["Value"].tolist()
                Npsteps=self.modelList(psteps,types)
                Nvsteps= [abs(i*j) for i,j in zip(Npsteps,values)] 
                dfn["Npsteps"]=Npsteps
                dfn["Nvsteps"]=Nvsteps
                
                # buy dataframe of a day
                dfnb = dfn[dfn["Type"]=='Buy']
                # sell dataframe of a day
                dfns = dfn[dfn["Type"]=='Sell']
                
                # Time axis
                timeb = dfnb["Time"].tolist()
                times = dfns["Time"].tolist()
                first_time = datetime.strptime(dfn["Time"].tolist()[0], '%H:%M:%S')
                timeb_ins = [datetime.strptime(i,'%H:%M:%S') for i in timeb]
                times_ins = [datetime.strptime(i,'%H:%M:%S') for i in times]
                tb_ins = [(i-first_time).total_seconds() for i in timeb_ins]
                ts_ins = [(i-first_time).total_seconds() for i in times_ins]
                
                time = dfn["Time"].tolist()
                time_ins = [datetime.strptime(i,'%H:%M:%S') for i in time]
                t_ins = [(i-first_time).total_seconds() for i in time_ins]
                
                # val/trade of a day
                vpt = dfn["ValperTrade"].tolist()
                # buy_val/trade of a day
                vptb = dfnb["ValperTrade"].tolist()
                # sell_val/trade of a day
                vpts = dfns["ValperTrade"].tolist()
                        
                # buy_val of a day
                valb = dfnb["Value"].tolist()
                # sell_val of a day
                vals = dfns["Value"].tolist()
                
                # NVsteps of a day
                nvstb = dfnb["Nvsteps"].tolist()
                # sell_val of a day
                nvsts = dfns["Nvsteps"].tolist()
               
                # Calculate COTP
                price_list = dfn["Price"].tolist()
                value_list = dfn["Value"].tolist()
                pval = [i*j for i, j in zip(price_list,value_list)]
                cotp = round(sum(pval)/sum(value_list),1)
                cotp = str(cotp)
            
                # calculate sum_val for a number of entries
                vptb_sum=[]; valb_sum = []; pvstb_sum=[]
                vpts_sum=[]; vals_sum = []; pvsts_sum=[]
                valb_frac = []
                idx = 1
                while idx<len(valb)+1:
                    vptbsum = self.sumListItems(num_entry, vptb, idx)
                    vptb_sum.append(vptbsum)
                    vbsum = self.sumListItems(num_entry, valb, idx)
                    valb_sum.append(vbsum)
                    pvstbsum=self.sumListItems(num_entry, nvstb, idx)
                    pvstb_sum.append(pvstbsum)
                    idx += 1
                idx = 1
                while idx<len(vals)+1:
                    vptssum = self.sumListItems(num_entry, vpts, idx)
                    vpts_sum.append(vptssum)
                    vssum = self.sumListItems(num_entry, vals, idx)
                    vals_sum.append(vssum)
                    pvstssum=self.sumListItems(num_entry, nvsts, idx)
                    pvsts_sum.append(pvstssum)
                    idx += 1
                    
                # calculate buy frac values for a number of entries 
                nentry = num_entry
                # buyfrac_vpt, bft_vpt = buyFrac(nentry, vpt, types, cut_val=1)
                buyfrac_val, bft_val = self.buyFrac(nentry, values , types, cut_val=30)
                # buyfrac_pvsteps, bft_pvsteps = buyFrac(nentry, Nvsteps, types,cut_val=5)
        
                # Set figure
                fig, ax = plt.subplots(4, sharex=True)
                # fig.set_size_inches(30, 20, forward=True)
                # ax.axis([0, 5000, 0,10])
            
                # ## Val/Trade
                # ax[0].scatter(ts_ins,vpts, marker='x', color='b')
                # ax[0].scatter(tb_ins,vptb, marker='x', color='k')
                # ax[0].set_ylabel('vpt')
                
                # # ax[1].scatter(ts_ins,vpts_sum, marker='x', color='b')
                # # ax[1].scatter(tb_ins,vptb_sum, marker='x', color='k')
                # # ax[1].set_ylabel('vpt_sum')
                
                # ax[1].scatter(t_ins,buyfrac_vpt, marker='.', color='k',s=markersize)
                # ax[1].set_ylabel('bfr_vpt')
                # ax[1].axhline(y=0.8, linewidth=2, color='g')
                # ax[1].axhline(y=0.5, color='k')
                # ax[1].axhline(y=0.2, color='r')
                # ax[1].set_ylim(0,1)
                
                # ## Pstep*Value
                # ax[2].scatter(ts_ins,nvsts, marker='x', color='b')
                # ax[2].scatter(tb_ins,nvstb, marker='x', color='k')
                # ax[2].set_ylabel('pvst')
                     
                # # ax[3].scatter(ts_ins,pvsts_sum, marker='x', color='b')
                # # ax[3].scatter(tb_ins,pvstb_sum, marker='x', color='k')
                # # ax[3].set_ylabel('PVst_sum')
                
                # # buy fraction of pvst
                # ax[3].scatter(t_ins,buyfrac_pvsteps, marker='.', color='k', s=markersize)
                # ax[3].set_ylabel('bfr_pvst')
                # ax[3].axhline(y=0.8, linewidth=2, color='g')
                # ax[3].axhline(y=0.5, color='k')
                # ax[3].axhline(y=0.2, color='r')
                # ax[3].set_ylim(0,1)
                
                ax[0].scatter(ts_ins,vals, marker='x', color='b')
                ax[0].scatter(tb_ins,valb, marker='x', color='k')
                ax[0].set_ylabel('value')
                
                ax[1].scatter(t_ins,buyfrac_val, marker='.', color='k', s=markersize)
                ax[1].set_ylabel('bfr_val')
                ax[1].axhline(y=0.8, linewidth=2, color='g')
                ax[1].axhline(y=0.5, color='k')
                ax[1].axhline(y=0.2, color='r')
                ax[1].set_ylim(0,1)
                
                ax[2].scatter(t_ins,bft_val, marker='.', color='k', s=markersize)
                ax[2].set_ylabel('bfrt_val')
                ax[2].axhline(y=0.7, linewidth=2, color='g')
                ax[2].axhline(y=0.5, color='k')
                ax[2].axhline(y=0.3, color='r')
                ax[2].set_ylim(0,1)
                
                # ax[5].scatter(ts_ins,vals_sum, marker='x', color='b')
                # ax[5].scatter(tb_ins,valb_sum, marker='x', color='k')
                # ax[5].set_ylabel('val_sum')
                
                ax[3].scatter(t_ins, price_list, marker = 'x', color='r')
                ax[3].set_ylabel('price')
                 
                # # Histograms of val/trade
                # ax.hist(vptb,50, ec='black', fc='none', lw=2, histtype='step',label='buy', 
                #           cumulative='True')
                # ax.hist(vpts,50, ec='blue', fc='none', lw=2, histtype='step',label='sell', 
                #           cumulative='True')
                
                # # Add price to the figure
                cp = dfn['Price'].tolist()[-1]
                cp = str(cp)
                # ax.text(6,450, cp)
                
                sname = dfn['Name'].tolist()[0]
                ax[0].set_title(sname+' '+date+' '+cp+' '+ cotp)
            
                
                
                # # Histograms of values
                # ax.hist(valb,50, ec='black', fc='none', lw=2, histtype='step',label='buy', 
                #          cumulative='True')
                # ax.hist(vals,50, ec='blue', fc='none', lw=2, histtype='step',label='sell', 
                #          cumulative='True')
                
            
        
                matplotlib.rc('font', **font)
                matplotlib.rcParams['axes.linewidth']=1
                
                fname = date
                if not midday_data:
                    plt.savefig(figpath+str(count).zfill(4))
                    plt.close()
                if midday_data:
                    plt.show()

                count += 1
              
                
            # Create movie
            # import os
            # import moviepy.video.io.ImageSequenceClip
            # image_folder=path
            # fps=1
            
            # image_files = [os.path.join(image_folder,img)
            #                for img in os.listdir(image_folder)
            #                if img.endswith(".png")]
            # clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
            # clip.write_videofile('my_video.mp4')
            
            
                    
