

import pandas as pd
import matplotlib.pyplot as plt

cutoffs = {'h_b' : 100, 'h_s': 70, 'm_b': 50, 'm_s' : 30}
name = 'EHL'

def scatPlot(df, cutoffs):
    name = name.upper()
    df = pd.read_csv(name+'.csv')
    h_b = cutoffs['h_b'] ;
    h_s = cutoffs['h_s'] ;
    m_b = cutoffs['m_b'] ;
    m_s = cutoffs['m_s'] ;
    Date = df['Date'].tolist(); Price = df['Price'].tolist(); Value = df['Value'].tolist(); Type = df['Type'].tolist()
    
                
    price = []; date = []; priceB = []; priceS = []; colorS = []; colorB = [];
    dateB = [];  dateS = []; dateMS = [];
    priceHB = []; dateHB = []; priceMB = []; dateMB = []
    priceHS = []; dateHS = [];  priceMS= []; dateMS =[]
    colorHB = []; colorHS = []; colorMB = []; colorMS = []
    dateP = []; pClose = [];
    k = max(Price)//10
    
    for i in range(len(df)):
        
        if i>1 and not Date[i-1] == Date[i]:
            if len(dateHB) == 0 or not dateHB[-1] == Date[i-1]:
                priceHB.append(Price[i-1]-k); dateHB.append(Date[i-1]); colorHB.append('white')  
            if len(dateHS) == 0 or not dateHS[-1] == Date[i-1]:
                priceHS.append(Price[i-1]-k); dateHS.append(Date[i-1]); colorHS.append('white'); 
            dateP.append(Date[i-1]); pClose.append(Price[i-1]);
            
        if Type[i] == 'Buy':
            u = 1
            if Value[i] > h_b:
                priceHB.append(Price[i]); dateHB.append(Date[i]); colorHB.append('yellow')
                
            elif m_b<=Value[i]<=h_b:
                priceMB.append(Price[i]); dateMB.append(Date[i]); colorMB.append('blue')
        else:
            if Value[i] > h_s:
                priceHS.append(Price[i]); dateHS.append(Date[i])
                colorHS.append('black')

            elif m_s<=Value[i]<=h_s:
                priceMS.append(Price[i]); dateMS.append(Date[i]); colorMS.append('white')
    
    dateP.append(Date[i-1]); pClose.append(Price[i-1])
    
    fig = plt.figure() 
    fig.set_figheight(10)
    fig.set_figwidth(20)

    hB = plt.scatter(dateHB, priceHB, c=colorHB,marker ="o", alpha=1, s= 250)
#     mB = plt.scatter(dateMB, priceMB, c=colorMB,marker ="s", alpha=0.4, s= 90)

    hS = plt.scatter(dateHS, priceHS, c=colorHS,marker ="^", alpha=1, s= 50)
#     mS = plt.scatter(dateMS, priceMS, c =colorMS, edgecolor = 'black', alpha = 0.5,marker =".", s = 80)
    Prc = plt.scatter(dateP, pClose, c = 'black', marker = 'x', alpha = 1, s =250 )

    plt.xticks(rotation=90, ha='right')
    plt.xlabel("Date")
    plt.ylabel("Price")

#     plt.legend((hB, mB, hS, mS, Prc),
#                ('HighBuy', 'MidBuy', 'HighSell', 'MidSell', 'Closing'),
#                scatterpoints=2,
#                loc='lower left',
#                ncol=1,
#                fontsize=14)
    plt.legend((hB, hS, Prc),
               ('HighBuy', 'HighSell', 'Closing'),
               scatterpoints=2,
               loc='lower left',
               ncol=1,
               fontsize=14)
    plt.show()