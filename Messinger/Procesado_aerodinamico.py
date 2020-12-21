

# %%
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
def Procesado_aerodinamico(V_inf,alpha):
    
    carpeta = 'Base_datos_NACA0012//'+str(int(V_inf))+'_'+str(int(round(alpha,0)))

    # Voy a la carpeta y leo los log para ver cual es la iteración idónea

    # %%
   
    file=open(carpeta+'//logs//UxFinalRes_0','r')
    UxRes=file.readlines()
    file.close()
    UxRes = np.loadtxt(carpeta+'//logs//UxFinalRes_0')[:,1]
    UyRes = np.loadtxt(carpeta+'//logs//UyFinalRes_0')[:,1]
    PRes = np.loadtxt(carpeta+'//logs//pFinalRes_0')[:,1]


    # %%
    Re = 1
    for i in range(UyRes.shape[0]):
        R = np.sqrt(UxRes[i]**2+UyRes[i]**2)
        if R<Re:
            iteracion = i+1
            
            Re=R

    iteracion=UyRes.shape[0]


    # %% [markdown]
    # ## Lectura de elementos

    # %%
    # os.listdir(carpeta+'//constant//polyMesh')
    # file = open(carpeta+'//constant//polyMesh')


    # %%
    file = open(carpeta+'//constant//polyMesh//boundary')
    lineas_contorno =file.readlines()
   
    file.close()
    numero_elementos =int(lineas_contorno[51][-11:-2])
    primer_elemento_superficie=int(lineas_contorno[52][-11:-2])

    # %% [markdown]
    # Lee todas las caras:

    # %%
    file = open(carpeta+'//constant//polyMesh//faces')
    lineas_contorno =file.readlines()
    file.close()
    # print(lineas_contorno[20+primer_elemento_superficie:20+primer_elemento_superficie+numero_elementos])
    puntos = []
    for cara in range(20+primer_elemento_superficie,20+primer_elemento_superficie+numero_elementos):
        elemento = lineas_contorno[cara][2:-2]+' '
        
        n=0
        for i in range(len(elemento)):
            if elemento[i]==' ':
                puntos.append(int(elemento[n:i]))
                n=i
                
    puntos_superficie =  set(puntos)  
    # print(set(puntos))


    # %%
    file = open(carpeta+'//constant//polyMesh//points')
    lineas_contorno =file.readlines()
    file.close()
    
    points =[]
    for j in puntos_superficie:
        punto = lineas_contorno[j+20][1:-2]+' '
        n=0
        coordenadas=[]
        for i in range(len(punto)):
            if punto[i]==' ':
                coordenadas.append(float(punto[n:i]))
                n=i
        points.append(coordenadas)
        
    points=pd.DataFrame(data =np.array(points),columns=['x','y','z'])
    points = points[(abs(points.y)<0.05)&(points.z>0)&(points.x<0.07)]
    points = points.sort_values(by=['y'])
    



    # %%
    import pickle 
    pickle.dump( points, open( "Superficie.p", "wb" ) )


    # %%
    file = open(carpeta+'//constant//polyMesh//faces')
    lineas_contorno =file.readlines()
    file.close()
    # print(lineas_contorno[20+primer_elemento_superficie:20+primer_elemento_superficie+numero_elementos])
    elementos=[]
 
    for cara in range(20,len(lineas_contorno)):
        if lineas_contorno[cara]==lineas_contorno[-4]:
            break
        elemento = lineas_contorno[cara][2:-2]+' '
        puntos = []
        n=0
        for i in range(len(elemento)):
            if elemento[i]==' ':
                puntos.append(int(elemento[n:i]))
                n=i
        elementos.append(puntos)
 


    # %%
    file = open(carpeta+'//constant//PolyMesh//owner')

    lineas_contorno =file.readlines()
    file.close()
    
    celdas=[]
    for linea in lineas_contorno[21:]:
        if linea==')\n':
            break
        celdas.append(int(linea))
   


    # %%
    Elementos = np.zeros((len(elementos),5))
    i=0
    for elemento in elementos:
        Elementos[i,4]=celdas[i]
        if len(elemento)==3:
            Elementos[i,3]=np.nan
        for j in range(len(elemento)):
            Elementos[i,j]= int(elemento[j])
        i+=1
   
    Elementos=pd.DataFrame(data=Elementos, columns=['punto_1','punto_2','punto_3','punto_4','celdas'])
  
    pickle.dump( points, open( "Elementos.p", "wb" ) )


    # %%
    pickle.dump( Elementos, open( "Elementos.p", "wb" ) )
    pickle.dump( points, open( "Puntos.p", "wb" ) )


    # %%
    file = open(carpeta+'//constant//polyMesh//points')
    lineas_contorno =file.readlines()
    file.close()
 
    points =[]
    for j in range(20,len(lineas_contorno)):
        if lineas_contorno[j]==lineas_contorno[-4]:
            break
        punto = lineas_contorno[j][1:-2]+' '
        n=0
        coordenadas=[]
        for i in range(len(punto)):
            if punto[i]==' ':
                coordenadas.append(float(punto[n:i]))
                n=i
        points.append(coordenadas)
        
    points=pd.DataFrame(data =np.array(points),columns=['x','y','z'])

    points = points[(points.z>0)]
    puntos_estudio = list(points.index)
   
    #(points.x,points.y,'o')
    Elementos = Elementos[(Elementos['punto_1'].isin(puntos_estudio))&(Elementos['punto_2'].isin(puntos_estudio))&(Elementos['punto_3'].isin(puntos_estudio))]
    # print(Elementos.head())
    # for i in Elementos.index:
        
    #     if points.loc[int(Elementos['punto_1'].loc[i])]['z']!=0:print('cuidado')
    #     if points.loc[int(Elementos['punto_2'].loc[i])]['z']!=0:print('cuidado')
    #     if points.loc[int(Elementos['punto_3'].loc[i])]['z']!=0:print('cuidado')
    #     if str(Elementos['punto_4'].loc[i])!='nan':
    #         if points.loc[int(Elementos['punto_4'].loc[i])]['z']!=0:print('cuidado')

    # %% [markdown]
    # # Procesado de Celdas

    # %%
    file = open(carpeta+'//constant//PolyMesh//owner')

    lineas_contorno =file.readlines()
    file.close()

    celdas=[]
    for linea in lineas_contorno[21:]:
        if linea==')\n':
            break
        celdas.append(int(linea))

    # %% [markdown]
    # # Procesado_velocidades

    # %%


    # %% [markdown]
    # Ahora lo que se hace es asignar las velocidades u y v a cada nodo. 

    # %%
    file = open(carpeta+'//'+str(iteracion)+'//U')
    lineas_contorno =file.readlines()
    file.close()
    # print(lineas_contorno[20+primer_elemento_superficie:20+primer_elemento_superficie+numero_elementos])
    elementos=[]
    # print(lineas_contorno[66873])
    for i in range(23,len(lineas_contorno)):
        
        elemento = lineas_contorno[i][1:-2]+' '
        puntos = []
        n=0
        if lineas_contorno[i][0]=='(' and lineas_contorno[i][1]!='\n':
            for i in range(len(elemento)):
                if elemento[i]==' ':
                    puntos.append(float(elemento[n:i]))
                    n=i
            elementos.append(puntos)

    elementos[-1]
    


    # %%
    Velocidades = pd.DataFrame(data=elementos,columns=['u','v','w'])
    Velocidades['id_elemento']=Velocidades.index
    #Velocidades=Velocidades[Velocidades['id_punto'].isin(puntos_estudio)]
    Velocidades
    pickle.dump( Velocidades, open( "Velocidades.p", "wb" ) )


    # %%
    Superficie = pickle.load( open( "Superficie.p", "rb" ) )
    Superficie['id_punto']=Superficie.index
    Superficie

    # %% [markdown]
    # # Dibujo de soluciones
    return (Elementos,Velocidades,Superficie,points)
