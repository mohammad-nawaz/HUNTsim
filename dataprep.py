# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 08:11:51 2022

@authors: ali.reza.buet@gmail.com
          ali.nawaz.md@gmail.com
"""
import os
import os.path
import pandas as pd
import numpy as np
import math
from numpy import nan
from datetime import date, time
from datetime import datetime
import tcalendar
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from readconfig import ReadConfig


class MinData:
    
    def minDataPrep(self, name, date):
        name = name.upper()
        td = tcalendar.TradingDates()
        date = td.dateStyle(date)
        calendarD = td.calendarDates('sep01', 'dec31')
        idx1 = calendarD.index('oct31')
        idx2 = calendarD.index(date)
        with open('paths.txt') as f:
            paths = [line.rstrip() for line in f]
        rawSN_store = paths[0]
        rawAS_store = paths[1]
        if idx2< idx1:
            df = self.minDataPrep1(name, date, rawSN_store)
            # if len(df) == 0:
            #     df = self.minDataPrep2(name, date, rawAS_store)
        else:
            df = self.minDataPrep2(name, date)
        return df
            
    def storeData(self, name, start_date, end_date):
        name = name.upper()
        td = tcalendar.TradingDates()
        start_date = td.dateStyle(start_date)
        end_date = td.dateStyle(end_date)
        
        calendarD = td.calendarDates('sep01', 'dec31')
        idx = calendarD.index('oct31')
        startidx = calendarD.index(start_date)
        endidx = calendarD.index(end_date)
        with open('paths.txt') as f:
            paths = [line.rstrip() for line in f]
        storage1 = paths[0]; storage2 = paths[1]
        finstore = paths[2];
        if startidx <= idx < endidx:
            
            df1 = self.minDataGroup(name, start_date, 'oct30', storage1)
            df1.reset_index(drop=True)
            df1.to_csv(finstore+name+'old.csv', index = False)
            
            df2 = self.minDataGroup(name, 'oct31', end_date, storage2)
            df2.reset_index(drop=True)
            df2.to_csv(finstore+name+'new.csv', index = False)
            df2 = self.minDataGroup(name, 'oct31', end_date, storage2, interpol = 'on')
            df2.reset_index(drop=True)
            df2.to_csv(finstore+name+'interpolated.csv', index = False)
            # elif endidx< idx:
            #     df = self.minDataGroup(name, start_date, end_date, storage1)
            #     df.to_csv(finstore+name+'1.csv')
        elif startidx <= idx and endidx<idx:
            df = self.minDataGroup(name, start_date, end_date, storage1)
            df.reset_index(drop=True)
            df.to_csv(finstore+name+'old.csv',index=False)
            
        elif startidx>=idx and endidx >= idx:
           df = self.minDataGroup(name,start_date, end_date, storage2);
           df.reset_index(drop=True)
           df.to_csv(finstore+name+'new.csv', index = False)
           df = self.minDataGroup(name,start_date, end_date, storage2, interpol = 'on');
           df.reset_index(drop=True)
           df.to_csv(finstore+name+'interpolated.csv', index = False)
        
            
    
    
    def minDataPrep2(self, name, date, interpol= None):
        name = name.upper()
        td = tcalendar.TradingDates()
        date = td.dateStyle(date)
        with open('paths.txt') as f:
            paths = [line.rstrip() for line in f]
        store_path = paths[1]
        head_name = store_path + date+ '_dvp.csv'
        file_exists = os.path.exists(head_name)
        if not file_exists:
            return print('file does not exist')
        df = pd.read_csv(head_name,low_memory=False)
        df = df.dropna(axis=1,how='all')
        cols = df.columns.tolist()
        find_name = '    "Scrip": "' + name +  '",'
        my_stock = []
        for col_name in cols:
            if df[col_name][2] == find_name:
                my_stock = df[col_name].dropna()
                break
        if len(my_stock) == 0:
            return pd.DataFrame()
        if interpol == 'on':
            dfn = self.minuteData(my_stock, name, date, interpol = 'on')
        
        else:
            dfn = self.minuteData(my_stock, name, date)
        # dfn = self.priceStep(dfn)
        # dfn = self.prevFormat(dfn)
        return dfn

    
   
    def minuteData(self, my_stock, name, date, interpol = None):
        '''convert rawdata_as to rawdataframe_as'''
        
        Name = []; Date = []; Time = []; Type =[]; Price = []; Volume = []; Trade = []; cumTrade= []; 
        Value = [];  ValperTrade = []; Color = []
		##
        step = []; p_step = []; v_step = []
		##
        Pdel =[]; Pdelsq =[]; Valdel =[]; Valdelsq =[]; Tmdiff = [];
        Cvpt = []; Bcvpt = []; Scvpt = []; Ncvpt = []
		
        gap = 16
        start = 6; end = len(my_stock)-11
        num_steps = 1+ (end-start)//gap 
        k = start
        flag = True
        curPrice = float(my_stock[end].replace('    "Close": ',"").replace(',', ""))
        initPrice = float(my_stock[start].replace('    "Close": ',"").replace(',', ""))
        tb,ts, tn =0,0,0
        bvpt,svpt,nvpt = 0,0,0
        for i in range(num_steps):
            
            ln = len(Price)
			
            #price, volume, trade, time
            p = float(my_stock[k].replace('    "Close": ',"").replace(',', ""))
            v = int(my_stock[k+3].replace('    "Volume": ',"").replace(',', ""))
            t = int(my_stock[k+5].replace('    "Trade": ',"").replace(',', ""))
            cumTrade.append(t)
            val = p*v/100000;
            if not i == 0:
                t = t-cumTrade[i-1]
            tm = my_stock[k+7][-10:].replace('",',"")

            
            if flag and p == initPrice:
                perc = (curPrice - p)/p
                
                bv = max(int((v//2) + perc*(v//2)),1)
                sv = max(v - bv,1)
                
                bt = max(int((t//2) + (perc*t//2)),1)
                st = max(t -bt,1)
                
                bval = round(p*bv/100000,10); sval = round(p*sv/100000,10);
                #######     buy      ########
                Name.append(name); Date.append(date); Time.append(tm)
                Price.append(p); Volume.append(bv);  Type.append('Buy'); Color.append('yellow');
                Trade.append(bt); Value.append(bval); ValperTrade.append(round(bval/max(bt,1), 5));
                # bvpt += bval; tb += bt;
                # Bcvpt.append(round(bvpt/max(tb,1), 4))
                # Scvpt.append(0)
                # Ncvpt.append(0)
                # Pdel.append(0); Pdelsq.append(0);
                ##########   sell     ##########\
                Name.append(name); Date.append(date); Time.append(tm)
                Price.append(p); Volume.append(sv); Type.append('Sell'); Color.append('red');
                Trade.append(st); Value.append(sval); ValperTrade.append(round(sval/max(st,1), 4));
                # Pdel.append(0); Pdelsq.append(0);
                # svpt += sval; ts += st;
                # Scvpt.append(round(svpt/max(ts,1), 4))
                # Bcvpt.append(0)
                # Ncvpt.append(0)
                
            else:
                flag = False
                Name.append(name); Date.append(date); Time.append(tm); 
                Price.append(p); Volume.append(v); Trade.append(t); 
                val = round(p*v/100000,4); Value.append(val); ValperTrade.append(round(val/max(t,1), 4))
                if p == Price[ln-1]:
                    Type.append('Neutral')
                    Color.append('Blue')
                    
                elif p>Price[ln-1]:
                    Type.append('Buy')
                    Color.append('yellow')
                else: 
                    Type.append('Sell')
                    Color.append('red')
				###
                # diffP = p - Price[ln-1]
                # Pdel.append(diffP)
                # Pdelsq.append(diffP**2)
                
                ##----------CVPT-------#########
                # if Type[ln] == 'Buy':
                #     bvpt += val
                #     tb += t
                #     Bcvpt.append(round(bvpt/max(tb,1),3))
                #     Scvpt.append(0)
                #     Ncvpt.append(0)
                # if Type[ln] == 'Neutral':
                #     nvpt += val
                #     tn += t
                #     Ncvpt.append(round(nvpt/max(tn,1),3))
                #     Bcvpt.append(0)
                #     Scvpt.append(0)
                # if Type[ln] == 'Sell':
                #     svpt += val
                #     ts += t
                #     Scvpt.append(round(svpt/max(ts,1),3))
                #     Bcvpt.append(0)
                #     Ncvpt.append(0)
                    

            k += gap
            #end of loop
        
        dfn = pd.DataFrame()
        dfn['Name'] = Name; dfn['Date'] = Date; 
        dfn['Time'] = Time; dfn['Type'] = Type;
        dfn['Price'] = Price; dfn['Volume'] = Volume; 
        dfn['Trade']= Trade; dfn['Value'] = Value; 
        dfn['ValperTrade'] = ValperTrade; 
        dfn['Color'] = Color;
		
        # dfn['CVPT'] = Cvpt
        # dfn['BCVPT'] = Bcvpt
        # dfn['SCVPT'] = Scvpt
        # dfn['NCVPT'] = Ncvpt
        dfn = self.priceStep(dfn)
        if interpol == 'on':
            dfn = self.interpolatedData(dfn)
            with open('paths.txt') as f:
                paths = [line.rstrip() for line in f]
            store_path = paths[1]
            dfn.to_csv(store_path+'temporaryfiletodelete.csv', index = False)
            dfn = pd.read_csv(store_path+'temporaryfiletodelete.csv')
            os.remove(store_path+'temporaryfiletodelete.csv')
            # dfn = dfn.drop(['Step','Pstep','Vstep','Tstep', 'Prate'], axis = 1)
            dfn = self.priceStep(dfn)
            dfn = self.addnewCols(dfn)
        return dfn
    
    def interpolatedData(self,df):
        df2 = pd.DataFrame()
        df3 = pd.DataFrame()
        idx1 = 0
        rows = []; R = []
        i = 0
        
        while i< len(df)-1:
            rows.append(i)
            i+=1
            if df["Type"][i] == 'Neutral':
                R.append(rows)
                rows =[]    
                while df['Type'][i] == 'Neutral':
                    if i == len(df)-1:
                        break
                    i += 1
        if len(R) == 0:
            return df
        elif len(R) == 1:
            df2= pd.DataFrame()
            df3 = df.iloc[R[0],:]
            idx1 = R[0][-1]
            idx2 = len(df)
            if df['Type'][idx1] == 'Buy':
                df2['Name'] = df['Name'][idx1+1:idx2]
                df2['Date'] = df['Date'][idx1+1:idx2]
                df2['Time'] = df['Time'][idx1+1:idx2]
                df2['Type'] = ['Buy' for x in range(idx1+1, idx2)]
                df2['Price'] = df['Price'][idx1+1:idx2]
                df2['Volume'] = df['Volume'][idx1+1:idx2]
                df2['Trade'] = df['Trade'][idx1+1:idx2]
                df2['Value'] = df['Value'][idx1+1:idx2]
                df2['ValperTrade'] = df['ValperTrade'][idx1+1:idx2]
                df2['Color'] =  ['yellow' for x in range(idx1+1, idx2)]
                df3 = pd.concat([df3,df2])
                df2 = pd.DataFrame()
                return df3
            else:
                df2['Name'] = df['Name'][idx1+1:idx2]
                df2['Date'] = df['Date'][idx1+1:idx2]
                df2['Time'] = df['Time'][idx1+1:idx2]
                df2['Type'] = ['Sell' for x in range(idx1+1, idx2)]
                df2['Price'] = df['Price'][idx1+1:idx2]
                df2['Volume'] = df['Volume'][idx1+1:idx2]
                df2['Trade'] = df['Trade'][idx1+1:idx2]
                df2['Value'] = df['Value'][idx1+1:idx2]
                df2['ValperTrade'] = df['ValperTrade'][idx1+1:idx2]
                df2['Color'] =  ['red' for x in range(idx1+1, idx2)]
                df3 = pd.concat([df3,df2])
                df2 = pd.DataFrame()
                return df3
        vs = df['Vstep'].tolist()
        vol = df['Volume'].tolist()
        for i in range(len(R)-1):            
            df2 = df.iloc[R[i],:]               
            df3 = pd.concat([df3,df2])
            df2 = pd.DataFrame()
            idx1 = R[i][-1];
            idx2 = R[i+1][0];
            if vs[idx1]> 0 and vs[idx2] > 0:
                df2['Name'] = df['Name'][idx1+1:idx2]
                df2['Date'] = df['Date'][idx1+1:idx2]
                df2['Time'] = df['Time'][idx1+1:idx2]
                df2['Type'] = ['Buy' for x in range(idx1+1, idx2)]
                df2['Price'] = df['Price'][idx1+1:idx2]
                df2['Volume'] = df['Volume'][idx1+1:idx2]
                df2['Trade'] = df['Trade'][idx1+1:idx2]
                df2['Value'] = df['Value'][idx1+1:idx2]
                df2['ValperTrade'] = df['ValperTrade'][idx1+1:idx2]
                df2['Color'] =  ['yellow' for x in range(idx1+1, idx2)]
                df3 = pd.concat([df3,df2])
                df2 = pd.DataFrame()
                # ro = []
                # for j in range(idx1+1, idx2):
                #    df['Type'][j] = 'Buy';
                #    df['Color'][j] = 'yellow'
                #    ro.append(j)
                # df2 = df.iloc[ro,:]
                # df3 = pd.concat([df3,df2])
                # df2 = pd.DataFrame()
            elif vs[idx1]< 0 and vs[idx2]< 0:
                df2['Name'] = df['Name'][idx1+1:idx2]
                df2['Date'] = df['Date'][idx1+1:idx2]
                df2['Time'] = df['Time'][idx1+1:idx2]
                df2['Type'] = ['Sell' for x in range(idx1+1, idx2)]
                df2['Price'] = df['Price'][idx1+1:idx2]
                df2['Volume'] = df['Volume'][idx1+1:idx2]
                df2['Trade'] = df['Trade'][idx1+1:idx2]
                df2['Value'] = df['Value'][idx1+1:idx2]
                df2['ValperTrade'] = df['ValperTrade'][idx1+1:idx2]
                df2['Color'] =  ['red' for x in range(idx1+1, idx2)]
                df3 = pd.concat([df3,df2])
                df2 = pd.DataFrame()
                # ro = []
                # for j in range(idx1+1, idx2):
                #    df['Type'][j] = 'Sell'
                #    df['Color'][j] = 'red'
                #    ro.append(j)
                # df2 = df.iloc[ro,:]
                # df3 = pd.concat([df3,df2])
                # df2 = pd.DataFrame() 
            
            elif vs[idx1]*vs[idx2] < 0:
                perc1 =  vs[idx1]/(abs(vs[idx1])+abs(vs[idx2]))
                perc2 =  vs[idx2]/(abs(vs[idx1])+abs(vs[idx2]))
                if perc1>0:
                    perc3 = perc1
                    perc1 = perc2
                    perc2 = perc3
                gap = idx2-idx1 -1
                dfnew = pd.DataFrame()
                Name = []; Date =[];Time =[]; Type=[];
                Price =[]; Volume =[]; Trade=[];Trade=[];
                Value=[]; ValperTrade= []; Color= []; Step =[]; Pstep = [];
                Vstep =[]; Tstep=[]; Prate=[]
                for j in range(idx1+1,idx2):
                    vol1 = max(int(abs(perc1)*vol[j]),1); vol2 = max(int(abs(perc2)*vol[j]),1)
                    val1 = round(df['Price'][j]*vol1/100000,5)
                    val2 = round(df['Price'][j]*vol2/100000,5)
                    tr1 = max(int(abs(perc1)*df['Trade'][j]),1);
                    tr2 = max(df['Trade'][j] -tr1,1)
                    Name.append(df['Name'][j])
                    Name.append(df['Name'][j])
                    Date.append(df['Date'][j])
                    Date.append(df['Date'][j])
                    Time.append(df['Time'][j])
                    Time.append(df['Time'][j])
                    Type.append('Sell')
                    Type.append('Buy')
                    Price.append(df['Price'][j])
                    Price.append(df['Price'][j])
                    Volume.append(vol1)
                    Volume.append(vol2)
                    Trade.append(tr1)
                    Trade.append(tr2)
                    Value.append(val1)
                    Value.append(val2)
                    ValperTrade.append(round(val1/tr1,5))
                    ValperTrade.append(round(val2/tr2,5))
                    Color.append('red')
                    Color.append('yellow')
                    Step.append(0); Step.append(0)
                    Vstep.append(0); Vstep.append(0)
                    Pstep.append(0); Pstep.append(0)
                    Tstep.append(0); Tstep.append(0)
                    Prate.append(0); Prate.append(0)
                
                dfnew['Name'] = Name; dfnew['Date'] = Date;
                dfnew['Time'] = Time; dfnew['Type'] = Type;
                dfnew['Price'] = Price; dfnew['Volume'] = Volume; dfnew['Trade'] = Trade
                dfnew['Value'] = Value; dfnew['ValperTrade']= ValperTrade;
                dfnew['Color'] = Color
                dfnew['Step'] = Step; dfnew['Pstep'] = Pstep;
                dfnew['Vstep'] = Vstep; dfnew['Tstep'] = Tstep;
                dfnew['Prate'] = Prate
                df3 = pd.concat([df3,dfnew])
        
        idx1 = R[-1][-1]
        idx2 = len(df)
        if df['Type'][idx1] == 'Buy':
            df2['Name'] = df['Name'][idx1+1:idx2]
            df2['Date'] = df['Date'][idx1+1:idx2]
            df2['Time'] = df['Time'][idx1+1:idx2]
            df2['Type'] = ['Buy' for x in range(idx1+1, idx2)]
            df2['Price'] = df['Price'][idx1+1:idx2]
            df2['Volume'] = df['Volume'][idx1+1:idx2]
            df2['Trade'] = df['Trade'][idx1+1:idx2]
            df2['Value'] = df['Value'][idx1+1:idx2]
            df2['ValperTrade'] = df['ValperTrade'][idx1+1:idx2]
            df2['Color'] =  ['yellow' for x in range(idx1+1, idx2)]
            df3 = pd.concat([df3,df2])
            df2 = pd.DataFrame()
        else:
            df2['Name'] = df['Name'][idx1+1:idx2]
            df2['Date'] = df['Date'][idx1+1:idx2]
            df2['Time'] = df['Time'][idx1+1:idx2]
            df2['Type'] = ['Sell' for x in range(idx1+1, idx2)]
            df2['Price'] = df['Price'][idx1+1:idx2]
            df2['Volume'] = df['Volume'][idx1+1:idx2]
            df2['Trade'] = df['Trade'][idx1+1:idx2]
            df2['Value'] = df['Value'][idx1+1:idx2]
            df2['ValperTrade'] = df['ValperTrade'][idx1+1:idx2]
            df2['Color'] =  ['red' for x in range(idx1+1, idx2)]
            df3 = pd.concat([df3,df2])
            df2 = pd.DataFrame()
        # df2 = df.iloc[ro,:]
        # df3 = pd.concat([df3,df2])
        # df2 = pd.DataFrame()
        return df3
    
    def addnewCols(self,df):
        Pbuy = [0 for x in range(len(df))]
        Psell = [0 for x in range(len(df))]
        # Pneut = [0 for x in range(len(df))]
        Vbuy, Vsell = [0 for x in range(len(df))], [0 for x in range(len(df))]
        Valb, Vals = [0 for x in range(len(df))], [0 for x in range(len(df))]
        Tbuy, Tsell = [0 for x in range(len(df))], [0 for x in range(len(df))]
        BCVPT, SCVPT = [0 for x in range(len(df))], [0 for x in range(len(df))]
        bt, st = 0, 0
        b,s = 0,0
        for i in range(len(df)):
            if df['Type'][i] == "Buy":
                Pbuy[i] = df['Price'][i]
                Vbuy[i] = df['Volume'][i]
                Valb[i] = df['Value'][i]
                Tbuy[i] = df['Trade'][i]
                b += Valb[i]
                bt += Tbuy[i]
                BCVPT[i] = (round(b/bt,4))
            elif df['Type'][i] == "Sell":
                Psell[i] = df['Price'][i]
                Vsell[i] = df['Volume'][i]
                Vals[i] = df['Value'][i]
                Tsell[i] = df['Trade'][i]
                s += Vals[i]
                st += Tsell[i]
                SCVPT[i] = (round(s/st,4))
            # else:
            #     Pneut[i] = df['Price'][i]
            #     Vneut[i] = df['Volume'][i]
            #     Valn[i] = df['Value'][i]
            #     Tneut[i] = df['Trade'][i]
                
        df['PriceB'] = Pbuy; df['VolumeB'] = Vbuy; df['ValueB'] = Valb; df['TradeB'] = Tbuy;
        df['PriceS'] = Psell; df['VolumeS'] = Vsell; df['ValueS'] = Vals; df['TradeS'] = Tsell;
        # df['PriceN'] = Pneut; df['VolumeN'] = Vneut; df['ValueN'] = Valn; df['TradeN'] = Tneut;
        df['BCVPT'] = BCVPT
        df['SCVPT'] = SCVPT
        return df
    
    
    
    def minDataGroup(self, name, start_date, end_date, storage, interpol = None):
        """
        Create minute data for multiple dates by concatenation
        """
        name = name.upper()
        td = tcalendar.TradingDates()
        tdates = td.tradingDates(name,start_date, end_date)
        result = pd.DataFrame()
        calendarD = td.calendarDates('sep01', 'dec31')
        with open('paths.txt') as f:
            paths = [line.rstrip() for line in f]
        storage1 = paths[0]; storage2 = paths[1]
        
        idx = calendarD.index('oct31')
        
        for dt in tdates:
            idxDt = calendarD.index(dt)
            if idxDt< idx:
                df = self.minDataPrep1(name,dt,storage1)
                print (df)
                result = pd.concat([result, df]) 
                result = result.reset_index(drop=True)
            else:
                if interpol == 'on':
                    df = self.minDataPrep2(name,dt,interpol='on')
                else:
                    df = self.minDataPrep2(name,dt) 
                print (df)
                result= pd.concat([result, df])  
                result = result.reset_index(drop=True) 
             
            
        return result
    
    def minDataUpdate(self, sname=None, enddate=None):
        """
        Concatenate new data to the existing minute data file
        """
        sname = sname.upper()
        with open('paths.txt') as f:
            paths = [line.rstrip() for line in f]
        finstore = paths[2]
        df = pd.read_csv(finstore+sname+"interpolated.csv")
        sname = df["Name"].tolist()[0]
        months=['jan','feb','mar','apr','may','jun',
                'jul','aug','sep','oct','nov','dec']
        
        today = date.today()
        m = months[today.month-1]
        d = today.day
        if d//10==0:
            d = '0'+str(d)
        else: 
            d = str(d)
            
        last_tdate = df["Date"].tolist()[-1]
        
        td = tcalendar.TradingDates()
        start_date = td.nextTradingDate(sname,last_tdate)

        if enddate==None: 
            end_date = m+d
        else: 
            end_date = enddate
        tdates = td.tradingDates(sname,start_date, end_date)
        dfu = pd.DataFrame()        
        dfu = pd.concat([dfu,df])
        
        for dt in tdates: 
            dfn = self.minDataPrep2(sname, dt, interpol = 'on')
            dfu = pd.concat([dfu,dfn])
        
        dfu = dfu.reset_index(drop=True)
        
        dfu = dfu.drop_duplicates(keep='first')
        
        dfu.to_csv(finstore+sname+"interpolated.csv",index=False)
            
        return dfu
    
    def minMiddayDataPrep(self,name):    
        """
        temporary addition of current data to analyse the latest market
        """
        name = name.upper()
        file_name = "temp.csv"
        file_exists = os.path.exists(file_name)
        if not file_exists:
            return print('file does not exist')
        
        df = pd.read_csv('temp.csv',low_memory=False)
        df = df.dropna(axis=1,how='all')
        cols = df.columns.tolist()
        cols = cols[::-1]
        find_name = '    "Scrip": "' + name +  '",'
        my_stock = []
        for col_name in cols:
            if df[col_name][2] == find_name:
                my_stock = df[col_name].dropna()
                break
        if len(my_stock) == 0:
            return pd.DataFrame()
        
        # date
        months=['jan','feb','mar','apr','may','jun',
                'jul','aug','sep','oct','nov','dec']
        today = date.today()
        m = months[today.month-1]
        d = today.day
        if d//10==0:
            d = '0'+str(d)
        else: 
            d = str(d)
        midday_date = m+d
        
        dfn = self.minuteData(my_stock, name, midday_date, interpol='on')
        return dfn
    
    def priceStep(self, df):
        """
        Calculate steps and update dataframe
        """
        step = []
        p_step = []
        v_step = []
        t_step = []
        for i in range(len(df)):
            if i==0: 
                step.append(0); p_step.append(0); v_step.append(0)
                t_step.append(0)
            else: 
                price_step = df['Price'][i]-df['Price'][i-1]
                if price_step>0:
                    step.append(1)

                elif price_step<0:
                    step.append(-1)
                elif price_step==0:
                    step.append(0)
                if not type(df['Time'][i-1]) == str:
                    x = df['Time'][i-1].strftime("%H:%M:%S") ; y = df['Time'][i].strftime("%H:%M:%S")   
                else:
                    x = df['Time'][i-1];  y = df['Time'][i]
                if x[0:5] == y[0:5]:
                    seconds = int(y[-2:])- int(x[-2:])
                elif y[0:2] == x[0:2]:
                    s1 = int(y[3:5]); s2 = int(x[3:5])
                    seconds = 60*(s1-s2) + int(y[-2:]) - int(x[-2:])
                elif not y[0:2] == x[0:2]:
                    seconds =3600*(int(y[0:2])-int(x[0:2])-1) +60*(60- int(x[3:5])) + 60*(int(y[3:5])) + int(y[-2:]) - int(x[-2:])
                
                t_step.append(seconds)

                p_step.append(price_step)
                v_step.append(price_step*df['Value'][i])
        df['Step'] = step 
        df['Pstep'] = p_step
        df['Vstep'] = v_step
        df['Tstep'] = t_step
        df['Prate'] = [round(m/max(1,n),3) for m,n in zip(p_step,t_step)]
        return df 
    
    ##############################################################
    def minDataPrep1(self,name, date, store_path=None):    
        """
        Create raw minute dataframe from raw data of stocknow
        """
        with open('paths.txt') as f:
            paths = [line.rstrip() for line in f]
        if store_path == None:
            store_path = paths[0]
        name = name.upper()
        name_list = name
        name += " - Latest Trades"
        tcal = tcalendar.TradingDates()
        head_name = tcal.dateStyle(date)
        head_name = store_path+ head_name + "_dvp.xlsx"
        file_exists = os.path.exists(head_name)
        if not file_exists: 
            print('file does not exist')
        if file_exists:
            finder = ["invested", "gainer", "loser"]
            for s_name in finder:
                df = pd.read_excel(head_name, sheet_name = s_name)
                data_raw = []
                if name in df.head(0):
                    data_raw = df[name]
                    break
            if len(data_raw) == 0:
                # return pd.DataFrame(), pd.DataFrame()
                return pd.DataFrame()
            data=data_raw.dropna(how='all')
        else: return pd.DataFrame() #, pd.DataFrame()
        
        data = self.readFirstData(data)

        Time = []; Type = [] ; Price = []; Volume = []; Trade = []; Value = []; 
        ValperTrade = []; Color = []; Date = []; Name = []
        k = 7
        for i in range((len(data)-7)//6):    
            Time.append(data[k]); Type.append(data[k+1]); Price.append(data[k+2])
            Volume.append(data[k+3]); Trade.append(data[k+4]); Date.append(date); 
            Name.append(name_list)
            Value.append(data[k+2]*data[k+3]/100000)
            ValperTrade.append(data[k+2]*data[k+3]/(max(1,100000*data[k+4])))
            if data[k+1] == "Buy":
                Color.append("yellow")
            else: Color.append("red")
            k += 6
        #df frame
        df = pd.DataFrame()
        df['Name'] = Name
        df['Date'] = Date
        df['Time'] = Time[::-1]
        df['Type'] = Type[::-1]
        df['Price'] = Price[::-1]
        df['Volume'] = Volume[::-1]
        df['Trade'] = Trade[::-1]
        df['Value'] = Value[::-1]
        df['ValperTrade'] = ValperTrade[::-1]
        df['Color'] = Color[::-1]
        
        df = self.priceStep(df)
        
        return df
     
    def readFirstData(self, data):
        """
        read first data by analysing data data (stocknow)
        """
        count = 5
        dsize = len(data)
        
        while count<dsize:
            if isinstance(data[dsize-count], time):
                # print(data[dsize-count])
                count += 5
            else: 
                # print(data[dsize-count])
                count += -5
                # print (count)
                break
        drest = data.drop(data.index[dsize-count:])
        dmorning = data.drop(data.index[:dsize-count]).reset_index(drop=True)
        
        # If no morning data is available return the original data
        if len(dmorning)==0:
            return data
        
        # Price difference (in percentage) between the latest price and 
        # the morning price
        if isinstance(data[8],str):
            price_diff = round((drest[9]-dmorning[1])/dmorning[1]*100)
        else:
            price_diff = 0
        
        # Number of input in the morning data
        md_count = round(len(dmorning)/5)
        
        # Update morning data with derived 'buy' and 'sell' flags
        # dmorn_new = pd.DataFrame()
        dmorn_new = []
        for start in range(0, dmorning.shape[0], 5):
            # buy chunk
            dmorning_chunk = dmorning.iloc[start: start+5].reset_index(drop=True)
            price = dmorning_chunk[1]
            bvol =round(dmorning_chunk[2]/2+ price_diff/100.*dmorning_chunk[2]/2)
            btrade = max(1,round(dmorning_chunk[3]/2))
            bval =round(price*bvol/100000,2)
            dmorn_new.append(dmorning[0])
            dmorn_new.append('buy')
            dmorn_new.append(price)
            dmorn_new.append(bvol)
            dmorn_new.append(btrade)
            dmorn_new.append(bval)
            
            # sell chunk
            svol = dmorning_chunk[2] - bvol
            strade = max(1, dmorning_chunk[3]-btrade)
            sval = round(price*svol/100000, 2)
            dmorn_new.append(dmorning[0])
            dmorn_new.append('sell')
            dmorn_new.append(price)
            dmorn_new.append(svol)
            dmorn_new.append(strade)
            dmorn_new.append(sval)
            
        # Concatenate morning data to the rest 
        dmnew = pd.Series(dmorn_new)
        data_edited = pd.concat([drest,dmnew]).reset_index(drop=True)
   
        return data_edited        
    ###################################################################

class DayData:
    """
    Data Preparation
    ...
    Attributes:
    -----------
    dataframe:  minute dataframe
    buycut:     Set lower cut for high buy values; default: 100
    sellcut:    Set lower cut for high sell values; default: 80
    ...
    Methods:
    --------
    dayData:
    weightedAverage:
    """
    
    def __init__(self, dataframe, buycut=None, sellcut=None):
        # Read configuration file
        rconf = ReadConfig()
        rconf.readConfig()
        
        self.dataframe = dataframe 
        self.buycut = rconf.hb_thresh
        self.sellcut = rconf.hs_thresh
        
        self.mbuycut = rconf.mb_thresh
        self.msellcut = rconf.ms_thresh
        
    def dayDataPrep(self, start_date=None, end_date=None, midday_data = None):
        if midday_data == 'on':
            # print(self.dataframe)
            # print (self.dataframe['Name'].tolist()[0])
            sname = self.dataframe['Name'].tolist()[0]
            md = MinData()
            df_mid = md.minMiddayDataPrep(name = sname)
            df = pd.concat([self.dataframe, df_mid])
            df = df.reset_index(drop=True)
            # print (self.dataframe)
            self.dataframe = df
        dat = self.dataframe["Date"].tolist()
        dates = [];[dates.append(x) for x in dat if x not in dates]
        
        # trading dates: 
        if start_date == None: 
            start_date = dates[0]
        if end_date == None: 
            end_date = dates[-1]
        td = tcalendar.TradingDates()
        tdates = td.tradingDates(sname,start_date, end_date)        
        
        
        day_count = []
        for i in range(len(tdates)):
            day_count.append(i+1)
        
        # cot : center of trade
        # tvalue : total value
        # tval_b : total buy value
        # tval_s : total sell value
        # hcotbp : high center of trade price for buy
        # hcotbv : high center of trade value for buy
        # hcotsp : high center of trade price for sell
        # hcotsv : high center of trade value for sell
        # bsp_rat: buy sell pressure ratio
        name = []
        op=[]; cp = []; hp=[]; lp = []
        cot =[]; tvalue = []; tval_b=[]; tval_s=[]
        hcotbp=[]; hcotbv=[]
        hcotsp=[]; hcotsv=[]
        bs_rat=[]; hbhs_rat=[]; mbms_rat = []
        bsp_rat=[]
        step_sum = []
        pstep_sum = []
        vstep_sum=[]
        for i in tdates:
            if i not in dates:
                cot.append(np.nan); tvalue.append(np.nan); tval_b.append(np.nan)
                op.append(np.nan); cp.append(np.nan); hp.append(np.nan); lp.append(np.nan)
                tval_s.append(np.nan); hcotbp.append(np.nan); hcotbv.append(np.nan)
                hcotsp.append(np.nan); hcotsv.append(np.nan); bs_rat.append(np.nan)
                hbhs_rat.append(np.nan); mbms_rat.append(np.nan)
                bsp_rat.append(np.nan); 
                step_sum.append(np.nan); pstep_sum.append(np.nan);
                vstep_sum.append(np.nan)
                name.append(self.dataframe["Name"].tolist()[0])
            else: 
                dftemp = self.dataframe[self.dataframe.Date==i]
                # print(dftemp)
                # Center of Trade
                cot.append(round(self.weightedAverage(dftemp, "Price", "Value"),1))
                # Total Buy value
                dftemp_b = dftemp[dftemp["Type"]=="Buy"]
                tval_b.append(int(dftemp_b["Value"].sum()))
                # Total Sell value
                dftemp_s = dftemp[dftemp["Type"]=="Sell"]
                tval_s.append(int(dftemp_s["Value"].sum()))          
                # Total value
                tvalue.append(int(dftemp["Value"].sum()))
                # High-center-of-trade- price and value for buy
                dftemp_hb = dftemp[(dftemp["Type"]=="Buy") & (dftemp["Value"]>=self.buycut)]
                # dftemp_mb = dftemp[(dftemp["Type"]=="Buy") & (dftemp["Value"]>=self.mbuycut)]
                hcotbprice = self.weightedAverage(dftemp_hb, "Price", "Value")
                hcotbp.append(int(hcotbprice))
                hcotbv.append(int(dftemp_hb["Value"].sum()))
                # High-center-of-trade- price and value for sell
                dftemp_hs = dftemp[(dftemp["Type"]=="Sell") & (dftemp["Value"]>=self.sellcut)]
                # dftemp_ms = dftemp[(dftemp["Type"]=="Sell") & (dftemp["Value"]>=self.msellcut)]
                hcotsprice = self.weightedAverage(dftemp_hs, "Price", "Value")
                hcotsp.append(int(hcotsprice))
                hcotsv.append(int(dftemp_hs["Value"].sum()))
                # name column
                name.append(self.dataframe["Name"].tolist()[0])
                # ratio between buy and sell values 
                bs_rat.append(round(dftemp_b["Value"].sum()/max(dftemp_s["Value"].sum(),50),1))
                # ratio between high buy and high sell values 
                hbhs_rat.append(round(dftemp_hb["Value"].sum()/max(dftemp_hs["Value"].sum(),50),1))
                
                # medium buy sell ratio
                
                dftemp_mb = dftemp[(dftemp["Type"]=="Buy") 
                                    & (dftemp["Value"]>=self.mbuycut)
                                    & (dftemp["Value"]<self.buycut)]
                dftemp_ms = dftemp[(dftemp["Type"]=="Sell") 
                                    & (dftemp["Value"]>=self.msellcut)
                                    & (dftemp["Value"]<self.buycut)]
                # ratio between medium buy and medium sell values 
                mbms_rat.append(round(dftemp_mb["Value"].sum()/max(dftemp_ms["Value"].sum(),50),1))
                
                # Calculate ratio betweeen buy and sell pressure points
                # Points: 100+ 10; 80-100: 7; 50-80: 5; 30-50: 3; 10-30: 1
                bsp_rat.append(self.buySellPratio(dftemp))
                
                
                step_sum.append(round(dftemp["Step"].sum()))
                pstep_sum.append(round(dftemp["Pstep"].sum()))
                vstep_sum.append(round(dftemp["Vstep"].sum()))
                
                
                # Opening, closing, high, and low prices
                op.append(round(dftemp["Price"].tolist()[0],1))
                cp.append(round(dftemp["Price"].tolist()[-1],1))
                hp.append(round(dftemp["Price"].max()))
                lp.append(round(dftemp["Price"].min()))
                
                                
        
        # Dataframe for day data
        dfday = pd.DataFrame({"Name":name,
                              "Date":tdates, "Day":day_count, 
                              "OP": op, "CP":cp, "HP":hp, "LP":lp,
                              "COTP":cot, "TVAL":tvalue, 
                              "TBV":tval_b, "TSV":tval_s, 
                              "HBP":hcotbp, "HBV":hcotbv,
                              "HSP":hcotsp, 
                              "HSV":hcotsv, 
                              "B/S":bs_rat, "HB/HS":hbhs_rat, "MB/MS":mbms_rat,
                              "BP/SP":bsp_rat,
                              "STEPS":step_sum,
                              "PSS":pstep_sum,
                              "VSS":vstep_sum})
        
        # Calculate remaining buy value (RBV) and 
        #           remaining high center of trade value (RHBV)
        curHCOTBV = dfday["HBV"].tolist()
        curHCOTSV = dfday["HSV"].tolist()
        curTVALB = dfday["TBV"].tolist()
        curTVALS = dfday["TSV"].tolist()
        
        # remHCOTBV : remaining high center of trade value for buy
        # remBV     : remaining buy value
        
        remHCOTBV=[]
        remBV = []
        remhcbv=0; rembv = 0 
            
        for i in range(len(curHCOTBV)):
            if np.isnan(curHCOTBV[i]):
                remhcbv = 0; rembv=0
                chbv = 0; chsv = 0
                ctbv = 0; ctsv = 0
            elif not np.isnan(curHCOTBV[i]):
                chbv = curHCOTBV[i]; chsv = curHCOTSV[i]
                ctbv = curTVALB[i]; ctsv = curTVALS[i]
               
            remhcbv = remhcbv + chbv - chsv
            rembv = rembv + ctbv - ctsv
            
            if rembv<=0: # Condition that all recent buys are expired
                remhcbv=0; rembv=0
                
            remHCOTBV.append(remhcbv)  
            remBV.append(rembv)
        
        dfday=dfday.assign(RHBV=remHCOTBV)
        dfday=dfday.assign(RBV=remBV)
        
        # Calculate: remainign high center of trade buy value difference
        #            with the previous day's value (RHBVD)
        #            remaining buy value difference with the previous day (RBVD)
        remHCOTBVD=[]
        remBVD = []
        RHBVdiff = 0
        RBVdiff = 0
        for i in range(len(remHCOTBV)):
            if i==0:
                remHCOTBVD.append(round(0,2)); remBVD.append(round(0,2))
            else:
                if np.isnan(dfday['COTP'].tolist()[i]):
                    remHCOTBVD.append(round(RHBVdiff,2))
                    remBVD.append(round(RBVdiff,2))  
                else:
                    remHCOTBVD.append(round(remHCOTBV[i]-remHCOTBV[i-1],2))
                    remBVD.append(round(remBV[i]-remBV[i-1],2))               
                RHBVdiff = remHCOTBV[i]-remHCOTBV[i-1]
                RBVdiff = remBV[i]-remBV[i-1]

        dfday=dfday.assign(RHBVD=remHCOTBVD)
        dfday=dfday.assign(RBVD=remBVD)         
   
        return dfday
    
    def buySellPratio(self, dframe = None, date = None):
        """
        """
        df = dframe 
        
        pres_point = []
        for i in df.Value.tolist():
            if i>=100: 
                pres_point.append(i*10)
            elif 100>i>=80:
                pres_point.append(i*5)
            elif 80>i>=50:
                pres_point.append(i*4)
            elif 50>i>=30:
                pres_point.append(i*3)
            elif 30>i>=10:
                pres_point.append(i*2)
            else: 
                pres_point.append(i*1)                
        
        df=df.assign(PPOINTS=pres_point)
        bp_point = df[df["Type"]=='Buy']["PPOINTS"].sum()
        sp_point = df[df["Type"]=='Sell']["PPOINTS"].sum()
        
        
        
        return round(bp_point/max(1,sp_point),1)
        
                
        
            
        
        
    def weightedAverage(self, dataframe, value, weight):
        val = dataframe[value]
        wt = dataframe[weight]
        if wt.sum()==0:
            wa = 0
        else:
            wa = (val*wt).sum()/wt.sum()
        return wa
    
    
class DayDataCheck:
    """
    """
    def dayDataCheck(self, sname, midday_data = None):
        sname = sname.upper()
        with open('paths.txt') as f:
            paths = [line.rstrip() for line in f]
        store = paths[2]
        sname = store + sname +'interpolated.csv'
        dfmin = pd.read_csv(sname)
        
        # Read configuration file
        rconf = ReadConfig()
        rconf.readConfig()
        
        # dprep = DayData(dfmin,buycut=rconf.hb_thresh,sellcut=rconf.hs_thresh)
        dprep = DayData(dfmin)
        dfd = dprep.dayDataPrep(midday_data = midday_data)
        
        return dfd 
        
        

           

# name = 'KOHINOOR'
# date = 'nov10'
# md = MinData()
# dfn = md.minDataPrep(name,date)
# dfn = md.minMiddayDataPrep(name)
