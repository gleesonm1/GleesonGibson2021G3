## code to look at and analyse Perple_X outputs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import ndimage as nd
import glob
import os
import csv

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def loadPerpleData(Name):

    Start=False

    with open('Data/'+Name,'r') as f:
        reader = csv.reader(f,delimiter=' ',skipinitialspace=True)
        for row in reader:
            A=row

            if Start==True:
                B=[float(i) for i in A]
                arr=np.array(B)
                if i==0:
                    Arr=arr
                else:
                    Arr=np.vstack((Arr,arr))
                i=i+1

            if A[0]=='T(K)' and A[1]=='P(bar)':
                AA=A
                Start=True
                i=0

    df=pd.DataFrame(data=Arr, columns=AA[:-1])

    return df

def loadAllPerpleData(Name):

    Start=False

    with open('Data/'+Name,'r') as f:
        reader = csv.reader(f,delimiter=' ',skipinitialspace=True)
        for row in reader:
            A=row

            if Start==True:
                #B=[float(i) for i in A]
                arr=np.array(A)
                if len(arr)>54:
                    arr=arr[:-1]

                if i==0:
                    Arr=arr
                else:
                    Arr=np.vstack((Arr,arr))
                i=i+1

            if A[0]=='Name' and A[1]=='Counter':
                AA=A[0:53]
                Start=True
                i=0

        # for Name in Names:
        #     B=np.zeros(len(arr),Grid**2)


    df=pd.DataFrame(data=Arr, columns=AA[:-1])

    return df

def SplitSolid(Min,Names,Factor,Divide):
    for i in range(0,len(Min[Names[0]]['wt,%'])):
        if Min[Names[0]][Factor].loc[i]==Min[Names[1]][Factor].loc[i]:
            if Min[Names[0]][Factor].loc[i]>Divide:
                Min[Names[1]].loc[i]=Min[Names[0]].loc[i]*0
                Min[Names[1]]['T(K)'].loc[i]=Min[Names[0]]['T(K)'].loc[i]
                Min[Names[1]]['P(bar)'].loc[i]=Min[Names[0]]['P(bar)'].loc[i]
            else:
                Min[Names[0]].loc[i]=Min[Names[0]].loc[i]*0
                Min[Names[0]]['T(K)'].loc[i]=Min[Names[1]]['T(K)'].loc[i]
                Min[Names[0]]['P(bar)'].loc[i]=Min[Names[1]]['P(bar)'].loc[i]
        else:
            A=Min[Names[0]].loc[i]
            B=Min[Names[1]].loc[i]
            if Min[Names[0]][Factor].loc[i]<Min[Names[1]][Factor].loc[i]:
                Min[Names[0]].loc[i]=B
                Min[Names[1]].loc[i]=A

    return Min

def PerpleMesh(df,Var):
    Index=np.where(df['P(bar)']>1)
    Length=Index[0][0]
    Length=Length
    i=0
    for j in range(0,Length):
        arr=np.array(df[Var].values[i:i+Length])#np.array(df[Var].loc[i:i+Length-1])
        arr=arr.reshape(1,Length)
        i=i+Length
        if j==0:
            Arr=arr
        else:
            Arr=np.vstack((Arr,arr))

    T=np.array(df['T(K)'].values[0:Length])-273.15
    n=np.linspace(0,(Length)**2,Length+1)
    P=np.array(df['P(bar)'].loc[n[0:-1]])/10
    X,Y=np.meshgrid(T,P)

    Arr[np.isnan(Arr)]=0
    Arr = nd.median_filter(Arr, size=3)


    return Arr, X, Y

def MinAbsent(wt,Names,Trange=None):
    A=np.zeros(np.shape(wt[Names[0]]))
    D=np.zeros(np.shape(wt[Names[0]]))
    for Name in Names:
        B=np.where(wt[Name]==0)
        A[B]=A[B]+1

    C=np.where(A==len(Names))
    D[C]=True
    E=np.where(A<len(Names))
    D[E]=False

    return D

def Traverse(Min,Pressure=None,Temperature=None,Target=None, Mineral=None):
    if Mineral is None:
        if Pressure is not None:
            P=find_nearest(Min['Melt']['P(bar)'],Pressure*10)
            A=pd.DataFrame()
            for key in Min.keys():
                A[key]=Min[key][Target][Min['Melt']['P(bar)']==P]
            A['T']=Min['Melt']['T(K)'][Min['Melt']['P(bar)']==P]-273.15

        if Temperature is not None:
            T=find_nearest(Min['Melt']['T(K)'],Temperature+273.15)
            A=pd.DataFrame()
            for key in Min.keys():
                A[key]=Min[key][Target][Min['Melt']['T(K)']==T]
            A['P']=Min['Melt']['P(bar)'][Min['Melt']['T(K)']==T]/10

    if Mineral is not None:
        Oxides={'SiO2': 60.04, 'TiO2': 79.866, 'Al2O3':101.96, 'Cr2O3': 151.99, 'FeO': 71.844, 'MgO': 40.3044, 'CaO':56.0774, 'Na2O': 61.9789, 'O2': 15.999}
        if Pressure is not None:
            P=find_nearest(Min[Mineral]['P(bar)'],Pressure*10)
            A=pd.DataFrame()
            for ox in Oxides:
                A[ox]=Min[Mineral][ox+',wt%'][Min[Mineral]['P(bar)']==P]
            A['T']=Min[Mineral]['T(K)'][Min[Mineral]['P(bar)']==P]-273.15

        if Temperature is not None:
            T=find_nearest(Min[Mineral]['P(bar)'],Temperature+273.15)
            A=pd.DataFrame()
            for ox in Oxides:
                A[ox]=Min[Mineral][ox+',wt%'][Min[Mineral]['T(K)']==T]
            A['P']=Min[Mineral]['P(bar)'][Min[Mineral]['T(K)']==T]/10



    return A

def ConvertWt(Min):
    Oxides={'SiO2': 60.04, 'TiO2': 79.866, 'Al2O3':101.96, 'Cr2O3': 151.99, 'FeO': 71.844, 'MgO': 40.3044, 'CaO':56.0774, 'Na2O': 61.9789, 'O2': 15.999}
    for key in Min.keys():
        B=np.zeros(len(Min['Melt']['T(K)']))
        for ox in Oxides:
            A=Min[key][ox+',wt%']*Oxides[ox]
            B=B+A

        for ox in Oxides:
            Min[key][ox+',wt%']=100*Min[key][ox+',wt%']*Oxides[ox]/B

    return Min