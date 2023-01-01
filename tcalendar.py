#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 18:03:44 2022

@author: nawaz
"""
import ast
import re

class TradingDates:
    """
    Produce trading dates from start and end dates.
    ...
    Attributes:
    -----------
    start_date: starting date 
    end_date: end date
    ...
    Methods:
    --------
    Weekends(year)
        Return a list of all weekends for a given year
    calenderDates(start_date,end_date)
        Return all dates between given start and end dates
    tradingDates(start_date,end_date)
        Return all trading days (weekdays) between given start and end dates
    ----------
    """
    
    
    def dateStyle(self, date):
        date_split = re.split(r'(\d+)', date)
        m = date_split[0]
        d = int(date_split[1])
        if d//10==0:
            d = '0' + str(d)
        else: 
            d = date_split[1]
        
        return m+d
    
    def Weekends(self, year=2022):
        # all days of the year
        cdates = self.calendarDates('jan01','dec31')
        # 1st friday of 2022: jan07
        # 1st saturday of 2022: jan01
        fri1 = 'jan07'
        sat1 = 'jan01'
        # list of fridays and saturdays
        fridays = cdates[cdates.index(fri1)::7]
        saturdays = cdates[cdates.index(sat1)::7]
        weekends = fridays + saturdays
        
        return weekends
    
    def Weekends2023(self, year=2023):
        cdates = self.calendarDates('jan01','sep30')
        # 1st friday of 2022: jan07
        # 1st saturday of 2022: jan01
        fri1 = 'jan06'
        sat1 = 'jan07'
        # list of fridays and saturdays
        fridays = cdates[cdates.index(fri1)::7]
        saturdays = cdates[cdates.index(sat1)::7]
        weekends = fridays + saturdays
        
        return weekends
    
    def checkyear(self, start_date, end_date):
        newmonths = ['jan','feb','mar','apr','may','jun',
                'jul','aug']
        first_month = re.split(r'(\d+)',start_date)[0]
        last_month = re.split(r'(\d+)',end_date)[0]
        if last_month in newmonths:
            if first_month in newmonths:
                return 2
            return 1
        
        return 2022
        
    
    def recordDates(self, name): 
        with open('recordDates.txt') as f:
            data = f.read()
        rD = ast.literal_eval(data)
        f.close()
        if name in rD.keys():
            recordDate = rD[name]
        else: recordDate = []
        return recordDate

    def calendarDates(self, start_date, end_date, year=None):
        months=['jan','feb','mar','apr','may','jun',
                'jul','aug','sep','oct','nov','dec']
        newmonths = ['jan','feb','mar','apr','may','jun',
                'jul','aug']
        month_size = {'jan':31,
                      'feb':28,
                      'mar':31,
                      'apr':30,
                      'may':31,
                      'jun':30,
                      'jul':31,
                      'aug':31,
                      'sep':30,
                      'oct':31,
                      'nov':30,
                      'dec':31}
        
            
        
        first_month = re.split(r'(\d+)',start_date)[0]
        last_month = re.split(r'(\d+)',end_date)[0]
        tmonths = [months[x] for x in 
                   range(months.index(first_month),months.index(last_month)+1)]
        if last_month in newmonths:
            # newend_date = end_date
            # end_date = 'dec31'
            if first_month not in newmonths:
                newlast_month = last_month
                last_month = 'dec'
                tmonths = [months[x] for x in 
                        range(months.index(first_month),months.index(last_month)+1)] + [newmonths[x] for x in range(newmonths.index(newlast_month)+1)]
            
            
        tdays = []
        for m in tmonths:
            if tmonths.index(m)==0:
                dstart = int(re.split(r'(\d+)',start_date)[1])
            else: 
                dstart = 1
            if tmonths.index(m)==len(tmonths)-1:
                dend = int(re.split(r'(\d+)',end_date)[1])+1
            else: 
                dend = month_size[m]+1
            
            for d in range(dstart,dend):
                tdate = self.dateStyle(m+str(d))
                tdays.append(tdate)
                
        return tdays
        
    def tradingDates(self,name, start_date,end_date):
        weekends = self.Weekends()
        check = self.checkyear(start_date, end_date)
        if not check==2022:
            weekends2023 = self.Weekends2023()
            if check == 1:
                cdates = self.calendarDates('jan01','dec31')
                fri1 = 'jan07'
                sat1 = 'jan01'
                fridays = cdates[cdates.index(fri1)::7]
                saturdays = cdates[cdates.index(sat1)::7]
                cd = self.calendarDates(start_date, 'dec31')
                idx = []
                for date in cd:
                    if date in fridays:
                        idx.append(cdates.index(date))
                    elif date in saturdays:
                        idx.append(cdates.index(date))
                    if len(idx)==2:
                        break
                if len(idx)==0:
                    fridays = []; saturdays = []
                elif len(idx) == 1:
                    fridays = cdates[idx[0]::7]
                    saturdays = []
                elif len(idx) == 2:
                    fridays = cdates[idx[0]::7]
                    saturdays = cdates[idx[1]::7]
                
                weekends = fridays + saturdays
                weekends += weekends2023
            else:
                weekends = weekends2023
            
        with open('holidays.txt') as f:
            lines = [line.rstrip() for line in f]
            f.close()
        lines.pop(0)
        holidays = lines
        record = self.recordDates(name) 
        
        tdates = [i for i in self.calendarDates(start_date,end_date)
                  if i not in weekends]
        
        tdates = [i for i in tdates
                  if i not in holidays]
        tdates = [i for i in tdates if i not in record]
        
        return tdates
    
    
    def nextTradingDate(self,name, date):
        
        date = self.dateStyle(date)
        if date in ['dec29','dec30','dec31']:
            record = self.recordDates(name)
            # print(record)
            if date == 'jan01':
                if 'jan01' in record:
                    return 'jan02'
            return 'jan01'
        
        
        tdates = self.tradingDates(name, date, 'dec31')
        cdates = self.calendarDates(date, 'dec31')
            
        # weekends = self.Weekends()
        # with open('holidays.txt') as f:
        #     lines = [line.rstrip() for line in f]
        #     f.close()
        # lines.pop(0)
        # holidays = lines
        
        # record = self.recordDates(name)
        
        # cdates = [i for i in self.calendarDates('jan01','dec31')]
        # tdates = [i for i in self.calendarDates('jan01','dec31')
        #           if i not in weekends
        #           if i not in holidays
        #           if i not in record]
        i = 0
        while i<10:
            i += 1
            if cdates[cdates.index(date)+i] in tdates:
                break    
        
        return cdates[cdates.index(date)+i]
    
    def prevTradingDate(self,name, date):
        
        date = self.dateStyle(date)
        record = self.recordDates(name)
        # print(record)
        if date == 'jan01':
            if 'dec29' in record:
                return 'dec28'
            return 'dec29'
        
        newmonths = ['jan','feb','mar','apr','may','jun',
                'jul','aug']
        month = date[0:3]
        if month not in newmonths:
            tdates = self.tradingDates(name, 'oct31', 'dec31')
            cdates = self.calendarDates('oct31',date)
        
        else:
            tdates = self.tradingDates(name, 'jan01', date)
            cdates = self.calendarDates('jan01',date)
        
        # weekends = self.Weekends()
        # with open('holidays.txt') as f:
        #     lines = [line.rstrip() for line in f]
        #     f.close()
        # lines.pop(0)
        # holidays = lines
        # record = self.recordDates(name)
        
        # cdates = [i for i in self.calendarDates('jan01','dec31')]
        # tdates = [i for i in self.calendarDates('jan01','dec31')
        #           if i not in weekends
        #           if i not in holidays
        #           if i not in record]
        
        i = 0
        while i<10:
            i += 1
            if cdates[cdates.index(date)-i] in tdates:
                break   
        
        return cdates[cdates.index(date)-i]
        


# td = TradingDates()

# tdates = td.tradingDates('sep1','oct13')
