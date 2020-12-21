import os

import pandas as pd
import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import Messinger
from Messinger import fem_velocidades  
import Messinger.Procesado_aerodinamico
from Messinger.Procesado_aerodinamico import Procesado_aerodinamico
from Messinger.Procesado_termico import procesado_termico,calcular_fraccion_congelacion
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pywt
import tkinter as tk


import numpy as np
import pandas as pd

from scipy.signal import find_peaks, peak_prominences
import scipy
from tkinter import filedialog
from tkinter import *
#%%
V_inf = int(70)
os.chdir('Messinger')
alpha=10

(Elementos,Velocidades,Superficie,Puntos)=Procesado_aerodinamico(70,10)
print(os.getcwd())
Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie,Puntos)
Modelo_termico = fem_velocidades.analisis_termico()
Modelo_termico.punto_remanso(Modelo)
print(Modelo_termico.x_remanso)


# In[23]:

df = pd.DataFrame(data =[[0,4,9,18,28,38,48,58],['A8','A7','A6','A5','A4','A3','A2','A1']])
df=df.T
df.columns = ['s(mm)','sensor']
print(df)


# In[16]:


MVD=40
T_remanso=-5
LWC=0.9
V=int(70)
x_experimental=[0,0.01,0.02]
T_experimental=[0,0.0,0.0]
angulo_ataque =5
zona='intrados'
(Elementos,Velocidades,Superficie,Puntos)=Procesado_aerodinamico(70,angulo_ataque)
Superficie = Superficie[Superficie['y']>=0]
Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie,Puntos)

x_superficie=np.array(Superficie['x'])
y_superficie = np.array(Superficie['y'])
s_perfil = [0]
for i in range(1,len(Modelo.x_superficie)):
    
    delta_S= (Modelo.x_superficie[i]-Modelo.x_superficie[i-1])**2+(Modelo.y_superficie[i]-Modelo.y_superficie[i-1])**2
    s_perfil.append(np.sqrt(delta_S)+s_perfil[-1])
X=[] 
df['s(m)'] =df['s(mm)']*10**-3  
for s in df['s(m)']:    
    for i in range(1,len(s_perfil)):     
            if s>=s_perfil[i-1] and s<=s_perfil[i]:
                X.append(Modelo.x_superficie[i-1]+(Modelo.x_superficie[i]-Modelo.x_superficie[i-1])/(s_perfil[i]-s_perfil[i-1])*(s-s_perfil[i-1]))
                break


# In[17]:
df['x(m)']=X
print(df)


# In[28]:


X={}
for i in range(len(df['x(m)'])-1):
    n_sensor=int(df['sensor'].loc[i][-1])
    X.update({'FBG'+str(n_sensor):df['x(m)'].loc[i]})
    

print(X)
#%%

fig,ax = plt.subplots(nrows=3,ncols=3,figsize =(30,15))
fig.tight_layout(pad=10.0)

columna = 0
fila = 0
Datos_rime=[]
for mvD in ['4','5','6']:
    fila = 0
    for lwC in ['A','B','C']:
        nombre_archivo ='C://Users//gonzalezvm//programa hielo//ensayos_P35_0_grados_glaze//ensayo_'+mvD+lwC+'.txt'
        file = open(nombre_archivo,'r')
        lineas = file.readlines()
        file.close()
        LWC = lineas[1]
        for caracter in range(len(LWC)):
            if LWC[caracter]=='=':
                lwc = float(LWC[caracter+1:-1])
        MVD = lineas[2]
        for caracter in range(len(MVD)):
            if MVD[caracter]=='=':
                mvd = int(MVD[caracter+1:-1])
        df = pd.read_csv(nombre_archivo,header=3,sep='\t',index_col=False)
        df = df[df['Tiempo (s)']<140]
        df.index =df['Tiempo (s)']
        df = df.drop(columns=['Tiempo (s)'])
        zona = 'extrados'
        for column in df.columns:
            if column=='FBG9':zona = 'intrados'
            if zona == 'intrados':
                print(int(column[-1]))
            else:
                ax[fila,columna].plot(df.index,df[column],label=column)
                ax[fila,columna].set_title(lineas[1]+lineas[2])
                ax[fila,columna].set_ylabel('T(ºC)')
                
       
        ax[fila,columna].legend()
        fila+=1
        
    ax[fila-1,columna].set_xlabel('x(m)')
    
    columna = columna + 1
#plt.savefig('glaze_10_grados.eps')   
plt.show()



#%%

fig,ax = plt.subplots(nrows=3,ncols=3,figsize =(30,15))

fig.tight_layout(pad=10.0)
columna = 0
fila = 0
Datos_rime=[]
for mvD in ['4','5','6']:
    fila = 0
    for lwC in ['A','B','C']:
        
        nombre_archivo ='C://Users//gonzalezvm//programa hielo//ensayos_P35_0_grados_glaze//ensayo_'+mvD+lwC+'.txt'
        file = open(nombre_archivo,'r')
        lineas = file.readlines()
        file.close()
        LWC = lineas[1]
        for caracter in range(len(LWC)):
            if LWC[caracter]=='=':
                lwc = float(LWC[caracter+1:-1])
        MVD = lineas[2]
        for caracter in range(len(MVD)):
            if MVD[caracter]=='=':
                mvd = int(MVD[caracter+1:-1])
        df = pd.read_csv(nombre_archivo,header=3,sep='\t',index_col=False)
        if mvD+lwC=='6B':df['Tiempo (s)'] = df['Tiempo (s)']+20
        df = df[df['Tiempo (s)']<140]
        df.index =df['Tiempo (s)']
        x_experimental = [0,0.01,0.02,0.039]
        T_experimental = [-2.5, -3,-5.5,-5.5]
        LWC=lwc
        V=70
        ax[fila,columna].set_title(str(lineas[1])+ str(MVD))
        T_remanso=np.mean(df['FBG1'].loc[70:80])
        MVD = mvd
        alpha=0
        zona_estudio = 'extrados'
        (x_ext,T_sur_ext)=procesado_termico(alpha,LWC,T_remanso,MVD,zona_estudio,V)
        for i in X:
            x = X[i]
            for j in range(len(x_ext)-1):
                if x>= x_ext[j] and x<= x_ext[j+1] and x>=0 and T_sur_ext[j]>=-200:
                    # ax[fila,columna].plot(x_ext[j],T_sur_ext[j],'ob')
                    # ax[fila,columna].text(x_ext[j],T_sur_ext[j],i)
                    Datos_rime.append([x_ext[j],T_sur_ext[j],i,mvd,lwc])
                    break

        
        zona_estudio = 'intrados'
        (x_int,T_sur_int)=procesado_termico(alpha,LWC,T_remanso,MVD,zona_estudio,V)
        for iteracion,T in enumerate(T_sur_int):
            
            if T <-200:T_sur_int[iteracion]=np.nan
        for i in X:
            x = X[i]
            for j in range(len(x_int)-1):
                if x<= x_int[j] and x>= x_int[j+1] and x<=0 and T_sur_int[j]>=-200:
                    # ax[fila,columna].plot(x_int[j],T_sur_int[j],'ob')
                    # ax[fila,columna].text(x_int[j],T_sur_int[j],i)
                    Datos_rime.append([x_int[j],T_sur_int[j],i,mvd,lwc])
                    break
      
        
        
        ax[fila,columna].plot(x_ext,T_sur_ext,'+k',label='analitic')
        ax[fila,columna].plot(x_int,T_sur_int,'+k')
        
        
        ax[fila,columna].set_xlabel('x(m)')
        ax[fila,columna].set_ylabel('T_sur (ºC)')
        ax[fila,columna].grid()
        
        Posiciones = []
        temperaturas =[]
        for sensor in X:
            try:
                temperaturas.append(float(np.mean(df[sensor].loc[70:80])))
            except:
                temperaturas.append(np.nan)
            Posiciones.append(X[sensor])
            ax[fila,columna].plot([Posiciones[-1],Posiciones[-1]],[temperaturas[-1],temperaturas[-1]],'ob',label='experimental')
        ax[fila,columna].plot([min(x_int),max(x_ext)],[T_remanso,T_remanso],'r',label='experimental')
        fila = fila + 1
        
    ax[fila-1,columna].set_xlabel('x(m)')
    columna = columna + 1
plt.savefig('ensayos_0_grados_raze.svg')   
plt.show()

#%%fig,ax = plt.subplots(nrows=3,ncols=3,figsize =(30,15))
fig.tight_layout(pad=10.0)

columna = 0
fila = 0
Datos_rime=[]
for mvD in ['4','5','6']:
    fila = 0
    for lwC in ['A','B','C']:
        nombre_archivo ='C://Users//gonzalezvm//programa hielo//ensayos_P35_0_grados_glaze//ensayo_'+mvD+lwC+'.txt'
        file = open(nombre_archivo,'r')
        lineas = file.readlines()
        file.close()
        LWC = lineas[1]
        for caracter in range(len(LWC)):
            if LWC[caracter]=='=':
                lwc = float(LWC[caracter+1:-1])
        MVD = lineas[2]
        for caracter in range(len(MVD)):
            if MVD[caracter]=='=':
                mvd = int(MVD[caracter+1:-1])
        df = pd.read_csv(nombre_archivo,header=3,sep='\t',index_col=False)
        df = df[df['Tiempo (s)']<140]
        df.index =df['Tiempo (s)']
        df = df.drop(columns=['Tiempo (s)'])
        zona = 'extrados'
        for column in df.columns:
            if column=='FBG9':zona = 'intrados'
            if zona == 'intrados':
                print(int(column[-1]))
            else:
                ax[fila,columna].plot(df.index,df[column],label=column)
                ax[fila,columna].set_title(lineas[1]+lineas[2])
                ax[fila,columna].set_ylabel('T(ºC)')
                
       
        ax[fila,columna].legend()
        fila+=1
        
    ax[fila-1,columna].set_xlabel('x(m)')
    
    columna = columna + 1
#plt.savefig('glaze_10_grados.eps')   
plt.show()



#%%

fig,ax = plt.subplots(nrows=3,ncols=3,figsize =(30,15))

fig.tight_layout(pad=10.0)
columna = 0
fila = 0
Datos_rime=[]
for mvD in ['4','5','6']:
    fila = 0
    for lwC in ['A','B','C']:
        
        nombre_archivo ='C://Users//gonzalezvm//programa hielo//P35_0_grados_rime//ensayo_'+mvD+lwC+'.txt'
        file = open(nombre_archivo,'r')
        lineas = file.readlines()
        file.close()
        LWC = lineas[1]
        for caracter in range(len(LWC)):
            if LWC[caracter]=='=':
                lwc = float(LWC[caracter+1:-1])
        MVD = lineas[2]
        for caracter in range(len(MVD)):
            if MVD[caracter]=='=':
                mvd = int(MVD[caracter+1:-1])
        df = pd.read_csv(nombre_archivo,header=3,sep='\t',index_col=False)
        if mvD+lwC=='6B':df['Tiempo (s)'] = df['Tiempo (s)']+20
        df = df[df['Tiempo (s)']<140]
        df.index =df['Tiempo (s)']
        x_experimental = [0,0.01,0.02,0.039]
        T_experimental = [-2.5, -3,-5.5,-5.5]
        LWC=lwc
        V=70
        ax[fila,columna].set_title(str(lineas[1])+ str(MVD))
        T_remanso=np.mean(df['FBG1'].loc[70:80])
        MVD = mvd
        alpha=0
        zona_estudio = 'extrados'
        (x_ext,T_sur_ext)=procesado_termico(alpha,LWC,T_remanso,MVD,zona_estudio,V)
        for i in X:
            x = X[i]
            for j in range(len(x_ext)-1):
                if x>= x_ext[j] and x<= x_ext[j+1] and x>=0 and T_sur_ext[j]>=-200:
                    # ax[fila,columna].plot(x_ext[j],T_sur_ext[j],'ob')
                    # ax[fila,columna].text(x_ext[j],T_sur_ext[j],i)
                    Datos_rime.append([x_ext[j],T_sur_ext[j],i,mvd,lwc])
                    break

        
        zona_estudio = 'intrados'
        (x_int,T_sur_int)=procesado_termico(alpha,LWC,T_remanso,MVD,zona_estudio,V)
        for iteracion,T in enumerate(T_sur_int):
            
            if T <-200:T_sur_int[iteracion]=np.nan
        for i in X:
            x = X[i]
            for j in range(len(x_int)-1):
                if x<= x_int[j] and x>= x_int[j+1] and x<=0 and T_sur_int[j]>=-200:
                    # ax[fila,columna].plot(x_int[j],T_sur_int[j],'ob')
                    # ax[fila,columna].text(x_int[j],T_sur_int[j],i)
                    Datos_rime.append([x_int[j],T_sur_int[j],i,mvd,lwc])
                    break
      
        
        
        ax[fila,columna].plot(x_ext,T_sur_ext,'+k',label='analitic')
        ax[fila,columna].plot(x_int,T_sur_int,'+k')
        
        
        ax[fila,columna].set_xlabel('x(m)')
        ax[fila,columna].set_ylabel('T_sur (ºC)')
        ax[fila,columna].grid()
        
        Posiciones = []
        temperaturas =[]
        for sensor in X:
            try:
                temperaturas.append(float(np.mean(df[sensor].loc[70:80])))
            except:
                temperaturas.append(np.nan)
            Posiciones.append(X[sensor])
            ax[fila,columna].plot([Posiciones[-1],Posiciones[-1]],[temperaturas[-1],temperaturas[-1]],'ob',label='experimental')
        ax[fila,columna].plot([min(x_int),max(x_ext)],[T_remanso,T_remanso],'r',label='experimental')
        fila = fila + 1
        
    ax[fila-1,columna].set_xlabel('x(m)')
    columna = columna + 1
plt.savefig('ensayos_0_grados_rime.svg')   


#%%Calcula la fraccion de congelacion F(t)
fig,ax = plt.subplots(nrows=3,ncols=3,figsize =(30,15))

fig.tight_layout(pad=10.0)
columna = 0
fila = 0
Datos_rime=[]
for mvD in ['4','5','6']:
    fila = 0
    for lwC in ['A','B','C']:
        
        nombre_archivo ='C://Users//gonzalezvm//programa hielo//P35_0_grados_rime//ensayo_'+mvD+lwC+'.txt'
        file = open(nombre_archivo,'r')
        lineas = file.readlines()
        file.close()
        LWC = lineas[1]
        for caracter in range(len(LWC)):
            if LWC[caracter]=='=':
                lwc = float(LWC[caracter+1:-1])
        MVD = lineas[2]
        for caracter in range(len(MVD)):
            if MVD[caracter]=='=':
                mvd = int(MVD[caracter+1:-1])
        df = pd.read_csv(nombre_archivo,header=3,sep='\t',index_col=False)
        if mvD+lwC=='6B':df['Tiempo (s)'] = df['Tiempo (s)']+20
        df = df[df['Tiempo (s)']<140]
        df.index =df['Tiempo (s)']
        x_experimental = np.linspace(0,0.039)
        T_experimental = np.linspace(0,0.039)
        LWC=lwc
        V=70
        ax[fila,columna].set_title(str(lineas[1])+ str(MVD))
        T_remanso=np.mean(df['FBG1'].loc[70:80])
        MVD = mvd
        alpha=0
        zona_estudio = 'extrados'
        (x_ext,f_sur_ext)=calcular_fraccion_congelacion(alpha,LWC,T_remanso,MVD,zona_estudio,V)
        # for i in X:
        #     x = X[i]
        #     for j in range(len(x_ext)-1):
        #         if x>= x_ext[j] and x<= x_ext[j+1] and x>=0 and T_sur_ext[j]>=-200:
        #             # ax[fila,columna].plot(x_ext[j],T_sur_ext[j],'ob')
        #             # ax[fila,columna].text(x_ext[j],T_sur_ext[j],i)
        #             #Datos_rime.append([x_ext[j],T_sur_ext[j],i,mvd,lwc])
        #             break

        
        zona_estudio = 'intrados'
        (x_int,f_sur_int)=calcular_fraccion_congelacion(alpha,LWC,T_remanso,MVD,zona_estudio,V)
        # for iteracion,T in enumerate(T_sur_int):
            
        #     if T <-200:T_sur_int[iteracion]=np.nan
        # # for i in X:
        #     x = X[i]
        #     for j in range(len(x_int)-1):
        #         if x<= x_int[j] and x>= x_int[j+1] and x<=0 and T_sur_int[j]>=-200:
        #             # ax[fila,columna].plot(x_int[j],T_sur_int[j],'ob')
        #             # ax[fila,columna].text(x_int[j],T_sur_int[j],i)
        #             Datos_rime.append([x_int[j],T_sur_int[j],i,mvd,lwc])
        #             break
      
        
        
        ax[fila,columna].plot(x_ext,f_sur_ext,'+k',label='analitic')
        ax[fila,columna].plot(x_int,f_sur_int,'+k')
        
        
        ax[fila,columna].set_xlabel('x(m)')
        ax[fila,columna].set_ylabel('f')
        ax[fila,columna].grid()
    
        fila = fila + 1
        
    ax[fila-1,columna].set_xlabel('x(m)')
    columna = columna + 1
plt.savefig('ensayos_0_grados_rime.svg')   


#%%

fig,ax = plt.subplots(nrows=3,ncols=3,figsize =(30,15))

fig.tight_layout(pad=10.0)
columna = 0
fila = 0
Datos_rime=[]
for mvD in ['4','5','6']:
    fila = 0
    for lwC in ['A','B','C']:
        
        nombre_archivo ='C://Users//gonzalezvm//programa hielo//ensayos_P35_0_grados_glaze//ensayo_'+mvD+lwC+'.txt'
        file = open(nombre_archivo,'r')
        lineas = file.readlines()
        file.close()
        LWC = lineas[1]
        for caracter in range(len(LWC)):
            if LWC[caracter]=='=':
                lwc = float(LWC[caracter+1:-1])
        MVD = lineas[2]
        for caracter in range(len(MVD)):
            if MVD[caracter]=='=':
                mvd = int(MVD[caracter+1:-1])
        df = pd.read_csv(nombre_archivo,header=3,sep='\t',index_col=False)
        if mvD+lwC=='6B':df['Tiempo (s)'] = df['Tiempo (s)']+20
        df = df[df['Tiempo (s)']<140]
        df.index =df['Tiempo (s)']
        x_experimental = np.linspace(0,0.039)
        T_experimental = np.linspace(0,0.039)
        LWC=lwc
        V=70
        ax[fila,columna].set_title(str(lineas[1])+ str(MVD))
        T_remanso=np.mean(df['FBG1'].loc[70:80])
        MVD = mvd
        alpha=0
        zona_estudio = 'extrados'
        (x_ext,f_sur_ext)=calcular_fraccion_congelacion(alpha,LWC,T_remanso,MVD,zona_estudio,V)
        # for i in X:
        #     x = X[i]
        #     for j in range(len(x_ext)-1):
        #         if x>= x_ext[j] and x<= x_ext[j+1] and x>=0 and T_sur_ext[j]>=-200:
        #             # ax[fila,columna].plot(x_ext[j],T_sur_ext[j],'ob')
        #             # ax[fila,columna].text(x_ext[j],T_sur_ext[j],i)
        #             #Datos_rime.append([x_ext[j],T_sur_ext[j],i,mvd,lwc])
        #             break

        
        zona_estudio = 'intrados'
        (x_int,f_sur_int)=calcular_fraccion_congelacion(alpha,LWC,T_remanso,MVD,zona_estudio,V)
        # for iteracion,T in enumerate(T_sur_int):
            
        #     if T <-200:T_sur_int[iteracion]=np.nan
        # # for i in X:
        #     x = X[i]
        #     for j in range(len(x_int)-1):
        #         if x<= x_int[j] and x>= x_int[j+1] and x<=0 and T_sur_int[j]>=-200:
        #             # ax[fila,columna].plot(x_int[j],T_sur_int[j],'ob')
        #             # ax[fila,columna].text(x_int[j],T_sur_int[j],i)
        #             Datos_rime.append([x_int[j],T_sur_int[j],i,mvd,lwc])
        #             break
      
        
        
        ax[fila,columna].plot(x_ext,f_sur_ext,'+k',label='analitic')
        ax[fila,columna].plot(x_int,f_sur_int,'+k')
        
        
        ax[fila,columna].set_xlabel('x(m)')
        ax[fila,columna].set_ylabel('f')
        ax[fila,columna].grid()
    
        fila = fila + 1
        
    ax[fila-1,columna].set_xlabel('x(m)')
    columna = columna + 1
plt.savefig('ensayos_0_grados_rime.svg')   
plt.show()