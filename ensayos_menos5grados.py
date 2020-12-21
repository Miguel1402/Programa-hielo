#%%
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
from Messinger.Procesado_termico import procesado_termico,procesado_termico_rapido
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pywt
import tkinter as tk
import numpy as np
from shapely.geometry import shape
import pandas as pd
from shapely.geometry import Point,LineString
from scipy.signal import find_peaks, peak_prominences
from scipy.optimize import curve_fit
import scipy
from tkinter import filedialog
from tkinter import *

V_inf = int(70)
os.chdir('Messinger')
#file = open('E://Ensayos_1_octubre//Ensayo_'+mvd+lwc+'_5.txt','r')
alpha=-5
(Elementos,Velocidades,Superficie,Puntos)=Procesado_aerodinamico(70,10)
print(os.getcwd())
Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie,Puntos)
Modelo_termico = fem_velocidades.analisis_termico()
Modelo_termico.punto_remanso(Modelo)
print(Modelo_termico.x_remanso)
df = pd.DataFrame(data =[[1,3,12.4,20.4,30.3,40.3,49,60],['A8','A7','A6','A5','A4','A3','A2','A1']])
df=df.T
df.columns = ['s(mm)','sensor']
print(df)
MVD=40
T_remanso=-5
LWC=0.9
V=int(70)
x_experimental=[0,0.01,0.02]
T_experimental=[0,0.0,0.0]
angulo_ataque = 10
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
df['x(m)']=X
print(df)
X={}
for i in range(len(df['x(m)'])-1):
    n_sensor=int(df['sensor'].loc[i][-1])
    X.update({'FBG'+str(n_sensor):-df['x(m)'].loc[i+1]})
    n_sensor=int(df['sensor'].loc[i][-1])+8
    X.update({'FBG'+str(n_sensor):df['x(m)'].loc[i+1]})
print(X)
(Elementos,Velocidades,Superficie,Puntos)=Procesado_aerodinamico(V,angulo_ataque)
Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie,Puntos)
x_experimental = [0,0.01,0.02,0.03,0.04]
T_experimental = [-2.5, -3,-5.5,-5.5,0]
def graficas_errores(sensor,stagnation):
    fig,ax = plt.subplots(nrows=3,ncols=3,figsize =(30,15))
    fig.tight_layout(pad=5.0)
    
    fig.suptitle(sensor)
    columna = 0
    fila = 0
    Datos_rime=[]
   
    output={}
    for mvD in ['4','5','6']:
        fila = 0
        for lwC in ['A','B','C']:
            # 'E://Angulo de ataque//RimeP7_'+mvD+lwC+'_5grados.txt'
            nombre_archivo='E://Angulo de ataque//RimeP7_'+mvD+lwC+'_5grados.txt'
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
            LWC=lwc
            V=70
            ax[fila,columna].set_title(str(lineas[1])+ str(MVD))
            T_remanso=np.mean(df['FBG1'].loc[110:120])
            MVD = mvd
                  
            error = []
            Temperatura_ref =float(np.mean(df[sensor].loc[85:95]))
            if Temperatura_ref >0:Temperatura_ref=0
            lwc_testing = np.arange(0.05,1.5,0.05)
            conv_coeff = np.arange(100,200,20)
            #se calcula el HCT real 
            mu = Modelo.viscosity((T_remanso+Temperatura_ref)/2+273.5)
            mu = 1e-5/(0.12764+124.38/((T_remanso+Temperatura_ref)/2+273.5))
           
            Modelo.set_presion_remanso(101325)
            Temperatura_estatica = T_remanso - V**2/(2*1004.5) +273.15
            P_st= Modelo. presion_estatica(Temperatura_estatica,V)
            rho_a = P_st/(287*Temperatura_estatica)
            k_a =-12.69 + 2.029 *np.sqrt(Temperatura_estatica) # cal/(hr m K)
            k_a = k_a*4.18/3600 #W m-2 K
            D= 0.0316*0.25
          
            Pr = 1004.5*mu/k_a
        
            h_c_remanso=1.14*(rho_a*V*D/mu)**0.5*Pr**0.4*k_a/D
            h_c_remanso=h_c_remanso/4.18
            
            output.update({mvD+lwC:(LWC,{})})
            colores = {20:'b',40:'g',70:'m'}
            print(h_c_remanso)
            
            for drop_size in [20,40,70]:
                K = 1000*(drop_size*1e-6)**2*V/(18*D*mu)
             
                Re_delta = V*drop_size*1e-6*rho_a/mu
                lambda_lambda_stokes = (0.8388+0.001483*Re_delta+.1847*Re_delta**0.5)**(-1)
                K_0 =1/8+lambda_lambda_stokes*(K-1/8)
                beta_0 = 1.4*(K_0-1/8)**0.84/(1+1.4*(K_0-1/8)**0.84)
                
                betha=pickle.load(open('Eficiencias_coleccion//betha'+str(V)+'_'+str(angulo_ataque)+'_'+str(drop_size)+'.p', "rb"))
                error = []
                for LWC_testing in lwc_testing: 
                    (x_int,T_sur_int)=procesado_termico_rapido(alpha,LWC_testing,T_remanso,drop_size,'extrados',V,Modelo,betha,h_c_remanso,stagnation,beta_0)
                    for i,x_intrados in enumerate(x_int):
        
                        if x_intrados>X[sensor]:
                            error.append(T_sur_int[i])
                            break
                ax[fila,columna].plot(lwc_testing,error,colores[drop_size]+'+-',label=f'MVD={drop_size} $\mu$m')
                linea = LineString([(lwc_testing[0],Temperatura_ref),(lwc_testing[-1],Temperatura_ref)])
                lineas_h_cte = LineString([(lwc_testing[i],error[i]) for i in range(len(error))])
                predicted_LWC= linea.intersection(lineas_h_cte)
                ax[fila,columna].plot([LWC,LWC],[min(error),max(error)],'k')
                ax[fila,columna].plot([min(lwc_testing),max(lwc_testing)],[Temperatura_ref,Temperatura_ref],'k')
                try:
                    output[mvD+lwC][1].update({drop_size:predicted_LWC.x})
                    ax[fila,columna].plot(predicted_LWC.x,predicted_LWC.y,colores[drop_size]+'o')
                    ax[fila,columna].plot([predicted_LWC.x,predicted_LWC.x],[min(error),max(error)],colores[drop_size])
                    
                except: 
                    output.update({mvD+lwC:(LWC,np.nan)})
            
            ax[fila,columna].legend()
            fila+=1
            print(LWC,MVD)
            

        ax[fila-1,columna].set_xlabel(f'LWC($g/m^3$)')
        columna = columna + 1
    
    return output  
Predicciones={}
for sensor in ['FBG16','FBG15','FBG14','FBG13']:
    Predicciones.update({sensor:graficas_errores(sensor,False)})
    
    break
print(Predicciones)
plt.show()  

#%%
from scipy.optimize import curve_fit
#Ahora se mira cual es coeficiente de transmision de calor convectivo
def T_equilibrium(LWC_testing,h_c_remanso):
    (x_int,T_sur_int)=procesado_termico_rapido(alpha,LWC_testing,T_remanso,drop_size,'extrados',V,Modelo,betha,h_c_remanso,stagnation,beta_0)
    for i,x_intrados in enumerate(x_int):     
        if x_intrados<X[sensor]:
            T=T_sur_int[i]
            
            break
    return T

stagnation = False
fig,ax = plt.subplots(nrows=1,ncols=1,figsize =(30,15))
fig.tight_layout(pad=10.0)
fig.suptitle(sensor)
columna = 0
fila = 0
Datos_rime=[]

output={}
for mvD in ['4','5','6']:
    fila = 0
    contenidos_agua_experimentales=[]
    Temperaturas_experimentales = []
    for lwC in ['A','B','C']:
        nombre_archivo='E://Angulo de ataque//RimeP7_'+mvD+lwC+'_5grados.txt'
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
        LWC=lwc
        contenidos_agua_experimentales.append(LWC)
        V=70
        T_remanso=np.mean(df['FBG1'].loc[110:120])
        MVD = mvd
        
        error = []
        Temperatura_ref =float(np.mean(df[sensor].loc[85:95]))
        if Temperatura_ref >0:Temperatura_ref=0
        Temperaturas_experimentales.append(Temperatura_ref)
    #se calcula el HCT real 
    mu = Modelo.viscosity((T_remanso+Temperatura_ref)/2+273.5)
    mu = 1e-5/(0.12764+124.38/((T_remanso+Temperatura_ref)/2+273.5))
    
    Modelo.set_presion_remanso(101325)
    Temperatura_estatica = T_remanso - V**2/(2*1004.5) +273.15
    P_st= Modelo. presion_estatica(Temperatura_estatica,V)
    rho_a = P_st/(287*Temperatura_estatica)
    k_a =-12.69 + 2.029 *np.sqrt(Temperatura_estatica) # cal/(hr m K)
    k_a = k_a*4.18/3600 #W m-2 K
    D= 0.0316*0.25
    Pr = 1004.5*mu/k_a
    colores = {20:'b',40:'g',70:'m'}
    drop_size=MVD
    betha=pickle.load(open('Eficiencias_coleccion//betha'+str(V)+'_'+str(angulo_ataque)+'_'+str(drop_size)+'.p', "rb"))
    K = 1000*(drop_size*1e-6)**2*V/(18*D*mu)
    Re_delta = V*drop_size*1e-6*rho_a/mu
    lambda_lambda_stokes = (0.8388+0.001483*Re_delta+.1847*Re_delta**0.5)**(-1)
    K_0 =1/8+lambda_lambda_stokes*(K-1/8)
    beta_0 = 1.4*(K_0-1/8)**0.84/(1+1.4*(K_0-1/8)**0.84)
    bounds =(0.03,1.5)
   
    coeficiente = np.arange(100,400,50)
    lwc_testing = lwc_testing = np.arange(0.05,1.5,0.05)
    
    ax.plot(lwc_testing,[T_equilibrium(lwc,100) for lwc in lwc_testing],colores[MVD]+'-',label=f'$MVD$={MVD} $\mu m$')
    ax.plot(contenidos_agua_experimentales,Temperaturas_experimentales,colores[MVD]+'o')
    
    print(MVD)
        

    ax.set_xlabel(f'LWC($g/m^3$)')
    columna = columna + 1
    
plt.show()
# %%
