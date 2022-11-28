#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 07:51:50 2022

@author: nawaz
"""

import dataprep
import configparser
import numpy as np
from readconfig import ReadConfig
import tcalendar
import pandas as pd

class HuntSim:
    """
    """
    
    def __init__(self, sname):
        self.sname = sname
        
        # Set initial configuration from file: 'simtrade_config.ini'
        # self.readConfig()
        rconf = ReadConfig()
        rconf.readConfig()
        
        # Set running trade initial variables
        self.curledg = rconf.i_inv
        self.invday1 = 0
        self.invday2 = 0
        self.invday3 = 0
        self.invmat = 0
        self.startdate = None
        self.enddate = None
        self.scount = 0
        
        # To calculate total g/l
        self.cinv = {}
        self.compinv = {}
        self.mscount = 0
        self.tscount = 0
        self.totinv = 0
        self.totsell = 0
        self.totgain = 0
        self.totgrowth = 0
        self.success = 0
        self.fail = 0
        self.successrate = 0
        
        # To calculate avg and bep
        self.avg = 0
        self.bep = 0
        
        # Updating the sell signals
        self.sellsig = 0
        
        # Trade Summary
        self.tradesummary = []

    
    def simTrade(self, curledg=None, matinv=None, immatinv=None, lastinvdate=None,
                 startdate=None, enddate=None, pout = False):
        ddc = dataprep.DayDataCheck()
        self.df = ddc.dayDataCheck(self.sname)
        self.tdates = self.df['Date'].tolist()
        
        rconf = ReadConfig()
        rconf.readConfig()
        
        # print('Date   #stock_count   COTP CBEP  RBVD Ledger   Total_inv   G/L')        
        for td in self.tdates:
            df_td = self.df[self.df['Date']==td]
            self.cotp = df_td['COTP'].item()
            
            self.cp = df_td['CP'].item()
                      
            # if np.isnan(df_td['COTP'].item()):
            if np.isnan(df_td['CP'].item()):
                # self.matureStock(td, call='buy')
                # self.matureStock(td, call='sell')
                continue
            self.rbvd = df_td['RBVD'].item()
            self.hbsrat = df_td['HB/HS'].item()
            self.mbsrat = df_td['MB/MS'].item()
            
            if pout:
                print ('********************* %s ******************'%td)
            self.countStock(td)
            self.calcGrowth(td)
            
            if pout:
                print(self.sname, td, self.tscount, self.mscount, self.cotp, 
                      round(self.bep,1), 
                      '(', self.hbsrat, self.mbsrat, self.rbvd, ')',
                      round(self.totinv/100000,2), round(self.totsell/100000,2), 
                      round(self.totgain), round(self.totgrowth,1))
                print('--------- BUY call ----------')            
        
            self.buyCall(td)
            
            self.countStock(td) 
            self.calcGrowth(td)

            if pout:
                print(self.sname, td, self.tscount, self.mscount, self.cotp, 
                      round(self.bep,1), 
                      '(', self.hbsrat, self.mbsrat, self.rbvd,  ')',
                      round(self.totinv/100000,2), round(self.totsell/100000,2), 
                      round(self.totgain), round(self.totgrowth,1))
            
            if pout:
                print('-------- SELL call ---------')
            self.sellCall(td)
            self.countStock(td)  
            self.calcGrowth(td)
            if pout: 
                print(self.sname, td, self.tscount, self.mscount, self.cotp, 
                      round(self.bep,1), 
                      '(', self.hbsrat, self.mbsrat, self.rbvd, ')',
                      round(self.totinv/100000,2), round(self.totsell/100000,2), 
                      round(self.totgain), round(self.totgrowth,1))
        
        
        # Summary of the trade for the given trading dates
        date = td
        total_inv = round(self.totinv/100000,2)
        total_sell = round(self.totsell/100000,2)
        total_gain = round(self.totgain)
        total_growth = round(self.totgrowth,1)
        success = self.success
        fail = self.fail
        self.trade_summary = [self.sname, self.tdates[0], self.tdates[-1],
                              total_inv, total_sell, total_gain,
                              total_growth, success, fail]
        
        # print(self.trade_summary)
            
        self.resetVal()
        
    def countStock(self, date):
        """
        Count total and mature stocks for a given date.

        Parameters
        ----------
        date : string

        Returns
        -------
        None.

        """
        totst_count = 0
        matst_count = 0
        for d in self.cinv.keys():
            if (self.tdates.index(date)-self.tdates.index(d))>=2:
                matst_count += self.cinv[d][0]
            totst_count += self.cinv[d][0]
        self.mscount = matst_count
        self.tscount = totst_count  
        
    def calcGrowth(self, date):
        tot_inv = 0
        tot_sell = 0
        tot_success = 0
        tot_fail = 0
        for d in self.compinv.keys():
            tot_inv += self.compinv[d][1]
            tot_sell += self.compinv[d][2]
            tot_success += self.compinv[d][5]
            tot_fail += self.compinv[d][6]
        tot_gain = tot_sell-tot_inv
        tot_growth = tot_gain/max(1, tot_inv)*100
        self.totinv = tot_inv
        self.totsell = tot_sell
        self.totgain = tot_gain
        self.totgrowth = tot_growth
        self.success = tot_success
        self.fail = tot_fail
        success_rate = round(tot_success/max(1,(tot_success+tot_fail))*100,2)
        self.successrate = success_rate
 
             
    
    def buyCall(self, date):
        df = self.df
        df_td = df[df['Date']==date]
        # td = tcalendar.TradingDates()
        # prev_date = td.prevTradingDate(name,date)
        # dates = df['Date'].tolist()
        
        # if date==dates[0]:
        #     prev_date = dates[0]
        # df_prev = df[df['Date']==prev_date]

        # rbv = df_td['RBV'].item()
        # rbv_prev = df_prev['RBV'].item()
        
        # rbv_diff = rbv-rbv_prev
        rbv = df_td['RBV'].item()
        rbv_diff = df_td['RBVD'].item()
        
        
        step_sum = df_td['STEPS'].item()
        pss = df_td['PSS'].item()
        vss = df_td['VSS'].item()
        
        
        rconf = ReadConfig()
        rconf.readConfig()
        
        # stock_count =  int(rconf.ind_inv/df_td['COTP']/1.004)
        stock_count =  int(rconf.ind_inv/df_td['CP']/1.004)
        
        prev_scount = self.tscount
        prev_avg = self.avg
        prev_bep = self.bep
        # buy_price = df_td['COTP'].item()*1.004
        buy_price = df_td['CP'].item()*1.004        
        cur_inv = stock_count*buy_price 
        
        # print(round(self.avg), round(self.bep), round(cur_inv))
        
        if rconf.Buy_params_count == 1: 
            if df_td['HB/HS'].item()>=rconf.b_hbs_ratio:                
                self.cinv[date] = [stock_count, buy_price, cur_inv]
                self.avg = (prev_avg*prev_scount+cur_inv)/(prev_scount+stock_count)
                self.bep = (prev_bep*prev_scount+cur_inv)/(prev_scount+stock_count)


        if rconf.Buy_params_count == 2: 
            if df_td['HB/HS'].item()>=rconf.b_hbs_ratio and df_td['MB/MS'].item()>=rconf.b_mbs_ratio:
                self.cinv[date] = [stock_count, buy_price, cur_inv]   
                self.avg = (prev_avg*prev_scount+cur_inv)/(prev_scount+stock_count)
                self.bep = (prev_bep*prev_scount+cur_inv)/(prev_scount+stock_count)

        if rconf.Buy_params_count == 3: 
            if df_td['HB/HS'].item()>=rconf.b_hbs_ratio and df_td['MB/MS'].item()>=rconf.b_mbs_ratio and rbv_diff>=rconf.b_rbv_diff:
                self.cinv[date] = [stock_count, buy_price, cur_inv]   
                self.avg = (prev_avg*prev_scount+cur_inv)/(prev_scount+stock_count)
                self.bep = (prev_bep*prev_scount+cur_inv)/(prev_scount+stock_count)
                
        if rconf.Buy_params_count == 4: 
            if df_td['HB/HS'].item()>=rconf.b_hbs_ratio and df_td['MB/MS'].item()>=rconf.b_mbs_ratio and rbv_diff>=rconf.b_rbv_diff and rbv<=rconf.b_rbv:
                self.cinv[date] = [stock_count, buy_price, cur_inv]   
                self.avg = (prev_avg*prev_scount+cur_inv)/(prev_scount+stock_count)
                self.bep = (prev_bep*prev_scount+cur_inv)/(prev_scount+stock_count)
                
        if rconf.Buy_params_count == 5: 
            if vss>=rconf.b_vss:
                self.cinv[date] = [stock_count, buy_price, cur_inv]   
                self.avg = (prev_avg*prev_scount+cur_inv)/(prev_scount+stock_count)
                self.bep = (prev_bep*prev_scount+cur_inv)/(prev_scount+stock_count)        
    
    def sellCall(self, date):
        df = self.df
        df_td = df[df['Date']==date]
        
        # Find previous trading dates sell signal
        td = tcalendar.TradingDates()
        prev_date = td.prevTradingDate(name,date)
        dates = df['Date'].tolist()
        
        if date==dates[0]:
            prev_date = dates[0]
        df_prev = df[df['Date']==prev_date]

        rbvd_prev = df_prev['RBVD'].item()
        

        
        rconf = ReadConfig()
        rconf.readConfig()
        
        self.countStock(date)
        prev_scount = self.tscount
        sell_count = self.mscount

        minv_avg = sell_count*self.avg   # mature investment with average value
        # sell_val = sell_count*df_td['COTP'].item()*0.995
        sell_val = sell_count*df_td['CP'].item()*0.995
        gain = sell_val - minv_avg
        growth = (gain/max(1,minv_avg))*100
        
        
        success = fail = 0
        if gain>0: 
            success = 1
        if gain<=0:
            fail = 1
        
        
        prev_bep = self.bep
        
        # Exit sellCall() if there is nothing to sell
        if sell_count ==0: 
            return 
        
        if rconf.Sell_params_count == 1: 
            if df_td['RBVD'].item()<rconf.s_rbv_diff:
                    
                # Sell mature stocks
                for d in list(self.cinv.keys()):
                    if (self.tdates.index(date)-self.tdates.index(d))>=2:
                        self.compinv[date] = [sell_count, minv_avg, sell_val, 
                                              gain, growth, success, fail]
                        del self.cinv[d]
                    
        if rconf.Sell_params_count == 2: 
            if df_td['RBVD'].item()<rconf.s_rbv_diff or rbvd_prev<rconf.s_rbv_diff:
                    
                # Sell mature stocks
                for d in list(self.cinv.keys()):
                    if (self.tdates.index(date)-self.tdates.index(d))>=2:
                        self.compinv[date] = [sell_count, minv_avg, sell_val, 
                                              gain, growth, success, fail]
                        del self.cinv[d]        
                    
            self.countStock(date)
            
            self.bep = (prev_bep*prev_scount-sell_val)/max(1,(prev_scount+self.tscount))
        
    # def matureStock(self, date, call):
    #     if call=='buy':
    #         if 2>(self.tdates.index(date)-self.tdates.index(self.invdate_latest))>=1:
    #             self.invday2 += self.invday1
    #             self.invday1 = 0
    #             self.invdate_previous = self.invdate_latest
                
    #         elif 2>(self.tdates.index(date)-self.tdates.index(self.invdate_previous))>=1:
    #             self.invday3 += self.invday2
    #             self.invday2 = 0
    #             self.invdate_previous2 = self.invdate_previous          
                
    #     if call=='sell': 
    #         if (self.tdates.index(date)-self.tdates.index(self.invdate_latest))>=2:
    #             self.invmat += self.invday1
    #             self.invday1 = 0 
            
    #         elif (self.tdates.index(date)-self.tdates.index(self.invdate_previous))>=2:
    #             self.invmat += self.invday2
    #             self.invday2 = 0

    #         elif self.invdate_previous2 and (self.tdates.index(date)-self.tdates.index(self.invdate_previous2))>=2:
    #             self.invmat += self.invday3
    #             self.invday3 = 0
                
    def resetVal(self):
        rconf = ReadConfig()
        rconf.readConfig()
        # Set running trade initial variables
        self.curledg = rconf.i_inv
        self.invday1 = 0
        self.invday2 = 0
        self.invmat = 0
        self.invdate_latest = None
        self.invdate_previous = None
        self.startdate = None
        self.enddate = None
        self.scount = 0
        
        # To calculate total g/l
        self.total_inv = 0    
        self.total_sell = 0
        self.total_gain = 0
        self.invest = 0
        self.cinv = {}
        self.compinv = {}
        self.mscount = 0
        self.tscount = 0
        self.totinv = 0
        self.totsell = 0
        self.totgain = 0
        self.totgrowth = 0
        self.success = 0
        self.fail = 0
        self.successrate = 0
        
        # To calculate avg and bep
        self.avg = 0
        self.bep = 0
        
        self.sellsig = 0
        
        # self.trade_summary = []
        
class SimDataOut:
    
    def simDataOut(self, stname):
        st = HuntSim(stname)
        st.simTrade()
        
        return st.trade_summary
    
    def simOutSummary(self, dfile):
        simout_list = []
        
        for i in dfile: 
            output = self.simDataOut(i)
            simout_list.append(output)
            print(output)
 

        simout = pd.DataFrame(simout_list, columns=['sname', 'start', 'end', 'tinv', 'tsell',
                                            'tgain', 'tgrowth', 'success', 'fail'])
        
        print (simout)
        
        numofst = len(simout)        
        total_inv = simout.tinv.sum()
        total_sell = simout.tsell.sum()
        total_gain = simout.tgain.sum()
        total_growth = total_gain/max(1,total_inv*100000)*100
        total_success = simout.success.sum()
        total_fail = simout.fail.sum()
        success_rate = total_success/(total_success+total_fail)*100
        
        simout_summary = pd.DataFrame([[numofst, total_inv, total_sell, 
                                        total_gain, total_growth, success_rate]]
                                      , 
                                      columns = ['#st', 'tinv', 'tsell',
                                                  'tgain','tgrowth', 'srate'])
        
        self.simout = simout
        self.simoutsum = simout_summary
        
        return self.simoutsum
        
        
        

# datafile = ["AAMRATECH", "ACMELAB","ACIFORMULA","ADNTEL","ANWARGALV", "BBS",
#             "BBSCABLES",
#             "BDCOM","BEACONPHAR","BEXIMCO","BPML","BSC","BSCCL","BXPHARMA",
#             "COPPERTECH","DELTALIFE","EHL","FAREASTLIF","GEMINISEA","GENEXIL",
#             "IBP", "INDEXAGRO",
#             "INTRACO","JHRML","JMISMDL","KDSALTD", "KEYACOSMET","KOHINOOR",
#             "LHBL","LRBDL",
#             "MAKSONSPIN","MALEKSPIN","METROSPIN","MIRAKHTER","MONOSPOOL", 
#             "NAHEEACP","NPOLYMER","OLYMPIC","ORIONINFU","ORIONPHARM", "PAPERPROC",
#             "PHARMAID",
#             "PENINSULA","SEAPEARL","SONALIPAPR","SONALIANSH","SPCERAMICS",
#             "SPCL","UNIQUEHRL"]

# Output simulation data
# sdo = SimDataOut()
# simout_list = []
# for i in datafile: 
#     output = sdo.simDataOut(i)
#     simout_list.append(output)

# simout = pd.DataFrame(simout_list, columns=['sname', 'start', 'end', 'tinv', 'tsell',
#                                     'tgain', 'tgrowth'])

# Simulation output summary
# sdo = SimDataOut()
# sdo.simOutSummary(datafile)
