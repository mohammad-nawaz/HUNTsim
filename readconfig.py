#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 11:40:09 2022

@author: nawaz
"""

import configparser

class ReadConfig:
    """
    """
            
    def readConfig(self):       
        config = configparser.ConfigParser()
        config.read('simtrade_config.ini')
        
        self.i_inv = float(config['Initial Investment']['Initial_inv'])
        self.ind_inv = float(config['Initial Investment']['Individual_inv'])
        self.Buy_params_count = float(config['Buy Condition']['Buy_params_count'])
        self.b_hbs_ratio = float(config['Buy Condition']['b_hbs_ratio'])
        self.b_mbs_ratio = float(config['Buy Condition']['b_mbs_ratio'])
        self.b_rbv_diff = float(config['Buy Condition']['b_rbv_diff']) 
        self.b_rbv = float(config['Buy Condition']['b_rbv']) 
        self.b_rbv_ratio = float(config['Buy Condition']['b_rbv_ratio'])
        self.b_rhbv_ratio = float(config['Buy Condition']['b_rhbv_ratio'])
        self.b_bs_ratio = float(config['Buy Condition']['b_bs_ratio'])
        self.b_hbv = float(config['Buy Condition']['b_hbv'])
        self.b_ss = float(config['Buy Condition']['b_ss'])
        self.b_pss = float(config['Buy Condition']['b_pss'])
        self.b_vss = float(config['Buy Condition']['b_vss'])
        
        self.Sell_params_count = float(config['Sell Condition']['Sell_params_count'])
        self.s_rbv_ratio = float(config['Sell Condition']['s_rbv_ratio'])
        self.s_rhbv_ratio = float(config['Sell Condition']['s_rhbv_ratio'])
        self.s_rbv_diff = float(config['Sell Condition']['s_rbv_diff'])
        
        self.hb_thresh = float(config['Threshold Values']['hb_thresh'])
        self.hs_thresh = float(config['Threshold Values']['hs_thresh'])
        self.mb_thresh = float(config['Threshold Values']['mb_thresh'])
        self.ms_thresh = float(config['Threshold Values']['ms_thresh'])
