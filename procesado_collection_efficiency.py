import pandas as pd
import numpy as np
import fem_velocidades  
import matplotlib.pyplot as plt
from matplotlib import cm
import pickle
import time
import seaborn as sns
from Procesado_aerodinamico import Procesado_aerodinamico
file = open('errores_betas.txt','w')
file.close()
for V_inf in [70]:

    for angulo_ataque in [-5]:
        alpha = np.pi/180*angulo_ataque
        #for mu_D in [50,60]:
        for mu_D in np.arange(20,40,70):
            try:
                print(mu_D)
                (Elementos,Velocidades,Superficie,Puntos)=Procesado_aerodinamico(70,angulo_ataque)
                Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie,Puntos)
                Modelo.set_T_remanso(273.15-5)
                Modelo.set_presion_remanso(1e5)
                x_0 = -0.1
                x_f = 0.08
                plt.figure(figsize =(15,10))
                # plt.plot(Modelo.x_superficie,Modelo.y_superficie,'o')
                beta_x = []
                U_d0 = V_inf*np.cos(alpha)
                V_d0 = V_inf*np.sin(alpha)
                #D =20e-6
                N=0
                print('calculando tiempo...')
                tiempo_inicial=time.time()
                # D_gotas.index = D_gotas['D(u)']
                for D in [mu_D]:
                    numero_gotas=500
                    D=D*10**-6    
                    n_nodos = 10000
                    N=0
                    for y_0 in np.linspace(-0.02,0.02,100):
                        print(N)
                        (x,y,U_d,V_d,t) = Modelo.trayectoria_gota(x_0,x_f,y_0,U_d0,V_d0,D,n_nodos)
                        
                        x_proyecccion =Modelo.proyeccion_gota(x_0,y_0,x_f,n_nodos)
                        ds = np.sqrt((y[-1]-y[0])**2+(x[-1]-x_proyecccion)**2)
                        dy = np.abs(y[-1]-y[0])
                        betha = dy/ds
                        if y[-1] >=0:
                            beta_x.append([D,x[-1],betha,'extrados',numero_gotas])
                        else:
                            beta_x.append([D,-x[-1],betha,'intrados',numero_gotas])
                        N=N+1
                # plt.savefig('trajectories_drops.eps')
                bethas = pd.DataFrame(data=beta_x,columns=['Diameter','x','beta','zona','n_gotas'])
                bethas=bethas.dropna()
                #bethas = bethas[(bethas['zona']=='intrados')&(bethas['x']>=0.001)]
                bethas['Volumen']=bethas['Diameter']**3
                bethas['Volumen']=bethas['Diameter']**3
                V_total = sum(bethas['Volumen'])
                bethas[abs(bethas.x)<0.06].plot.scatter(x='x',y='beta')
                bethas=bethas[abs(bethas.x)<0.06]
                betha_x = np.array(bethas[['x','beta']].sort_values(by=['x']))


                plt.plot(betha_x[:,0],betha_x[:,1])
                
                nombre='Eficiencias_coleccion//betha'+str(V_inf)+'_'+str(angulo_ataque)+'_'+str(mu_D)+'.p'
                pickle.dump( betha_x,open( nombre,"wb"))
                plt.savefig('Eficiencias_coleccion//Imagenes//betha'+str(V_inf)+'_'+str(angulo_ataque)+'_'+str(mu_D)+'.jpg')
                plt.close('all')
            except:
                file = open('errores_betas.txt','a')
                file.write(V_inf,angulo_ataque,mu_D)
                file.close()
        
