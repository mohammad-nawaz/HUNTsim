import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from numpy import sign


##df = pd.read_csv('EHL.csv')
##storepath = 'E:/EHL/'


def save_png(df, storepath):
    name = df['Name'][0]
    ###############
    Date = df['Date'].tolist(); Price = df['Price'].tolist(); Volume = df['Volume'].tolist();
    Value = df['Value'].tolist(); Type = df['Type'].tolist()

    valB ,volB ,pB = [],[],[]
    valS ,volS ,pS = [],[],[]

    cotB, cotS,coT, TB, TS, TF= [],[],[],[], [],[]
    V = []
    date = []

    if Type[0] =='Buy':
        valB.append(Value[0])
        volB.append(Volume[0])
        pB.append(Price[0])
        valS.append(0)
        volS.append(0)
        pS.append(0)
    else:
        valS.append(Value[0])
        volS.append(Volume[0])
        pS.append(Price[0])
        valB.append(0)
        volB.append(0)
        pB.append(0)

    k = 0
    for i in range(1,len(df)):
        if Type[i] =='Buy':
            valB.append(Value[i])
            volB.append(Volume[i])
            pB.append(Price[i])
            valS.append(0)
            volS.append(0)
            pS.append(0)
        else:
            valS.append(Value[i])
            volS.append(Volume[i])
            pS.append(Price[i])
            valB.append(0)
            volB.append(0)
            pB.append(0)

        if not Date[i-1] == Date[i]:
            date.append(Date[i-1])
            cotB.append(round(1e5*sum(valB[k:i])/sum(volB[k:i]), 2))
            cotS.append(round(1e5*sum(valS[k:i])/sum(volS[k:i]), 2))
            TB.append(sum(valB[k:i]))
            TS.append(sum(valS[k:i]))
            TF.append(sum(Value[k:i]))
            coT.append(round(1e5*sum(Value[k:i])/sum(Volume[k:i]), 2))
            V.append(Volume[k:i])
            k = i
    date.append(Date[i-1])
    cotB.append(round(1e5*sum(valB[k:i])/sum(volB[k:i]), 2))
    cotS.append(round(1e5*sum(valS[k:i])/sum(volS[k:i]), 2))
    coT.append(round(1e5*sum(Value[k:i])/sum(Volume[k:i]), 2))
    TB.append(sum(valB[k:i]))
    TS.append(sum(valS[k:i]))
    TF.append(sum(Value[k:i]))
    V.append(Volume[k:i])
    
    cTB = np.cumsum(TB).tolist()
    cTS = [0,0]+ np.cumsum(TS[2:]).tolist()
    c_cotB, c_cotS = [],[]
    uplim = int(max(Price)*1.1)
    lowlim = int(min(Price)*0.9)
    ################
    av_S, av_P = [],[]
    chunks = 7
    gap = (uplim-lowlim)//chunks
    Tsum = [0 for k in range(chunks)]
    Totsum = [0 for k in range(chunks)]

    for i in range(len(date)):

        #-------rem buy&sell value calc-------##

        if i>=2:
            s = cTS[i]
            av_S = TB[0:i-1]
            av_P = cotB[0:i-1]
            while len(av_S)>0:
                if s> av_S[0]:
                    s -= av_S[0]
                    av_S.pop(0)
                    av_P.pop(0)
                else:
                    av_S[0] -= s
                    break

        #######---------------------------######



        fig = plt.figure(facecolor='cyan')
        fig.set_figheight(10)
        fig.set_figwidth(20)
        fig.suptitle(f'Date: {date[i]} \n {name}', fontsize=20)


        gs1 = gridspec.GridSpec(1, 4)
        gs1.update(wspace=0.01, hspace=0.1)

        ###################--------------ax1-------------###############
        ax1 = plt.subplot(gs1[0])
        ax1.margins(0.05)          
        ax1.set_ylim([lowlim,uplim])

        y1 = cotB[i]; y2 = cotS[i]
        ycot = coT[i]

        ax1.scatter(0, y1, c='yellow',marker ="o", alpha=1, edgecolor = 'black', s= TB[i]//4)
        ax1.scatter(0, y2, c='red',marker ="o", alpha=0.5, s= TS[i]//4)

        if (y1-y2) == 0: k = y1//30
        else: k = sign(y1-y2)*y1//30

        ax1.annotate(f'buy:{int(TB[i])}', (0.005,y1+k))
        ax1.annotate(f'sell:{int(TS[i])}', (-0.025,y2-k))
        ax1.annotate(f'total:{int(TF[i])}', (-0.025,y1+k))
        ax1.axhline(y= ycot, xmin=0, xmax=1, c="black", linewidth=2, zorder=0)


        ax1.get_xaxis().set_visible(False)
        ax1.set_title(f'Buy&Sell at {date[i]}')

        ###########################--------------ax2-----------###############################
        ax2 = plt.subplot(gs1[1])
        ax2.set_ylim([lowlim,uplim])

        if i<2:
            ax2.scatter(0, 80.1, c='white', marker='.', s =1)

        elif len(av_S) == 0:
            ax2.scatter(0, lowlim+.1, c='white', marker='.', s =1)
            ax2.scatter(0, cotS[i], c='red',marker ="o", alpha=1, edgecolor = 'black', s= s//2)
            ax2.annotate(f'{int(s)}', (0.002,cotS[i]+cotS[i]//40))
        else:
            for j in range(len(av_S)):

                y1 = av_P[j]

                ax2.scatter(0, av_P[j], c='yellow',marker ="o", alpha=1, edgecolor = 'black', s= y1//2)

                ax2.annotate(f'{int(av_S[j])}', (0.002,y1+y1//40))

        ax2.axhline(y= ycot, xmin=0, xmax=1, c="black", linewidth=2, zorder=0)


        ax2.get_xaxis().set_visible(False)
        ax2.set_title('Rem Available Sellers')

        ##############------------------ax3-----------------#################
        ax3 = plt.subplot(gs1[2])

        ax3.set_ylim([lowlim,uplim])

        for j in range(chunks):
            if not Tsum[j] ==0:
                ax3.scatter(0, cotval[j], c = 'red', marker ='o', alpha= 0.5,s =Tsum[j]//100)
                ax3.annotate(f'{int(Tsum[j])}', (0.005,cotval[j]+cotval[j]//100))
                ax3.axhline(y= cotval[j], xmin=0, xmax=1, c="black", linewidth=2, ls = '-.', zorder=0)
        for j in range(chunks):
            if lowlim+j*gap  < coT[i] <= lowlim +(j+1)*gap:
                Tsum[j] += TF[i]
                Totsum[j] += coT[i]*TF[i]
        ax3.scatter(0, coT[i], c = 'blue', marker ='o', alpha= 1,s =80)
        ax3.annotate(f'{int(TF[i])}', (-0.015,coT[i]+coT[i]//100))
        ax3.axhline(y= coT[i], xmin=0, xmax=1, c="black", linewidth=2, zorder=0)
        cotval = [round(m/n,0) if n else 0 for m, n in zip(Totsum,Tsum)]

        ax3.get_xaxis().set_visible(False)
        ax3.set_title('available seller')   

        ################--------------ax4-------------#################
        ax4 = plt.subplot(gs1[3])
        ax4.set_ylim([lowlim,uplim])
        ax4.axhline(y= ycot, xmin=0, xmax=1, c="black", linewidth=2, zorder=0)
        ax4.set_title('Pending')

        plt.close()
 #       plt.show()
        if i//10==0:
            outputname = storepath + f'00{i}.png'
        
        else:
            outputname = storepath+ f'0{i}.png'
            
        fig.savefig(outputname,bbox_inches='tight')
      