import Messinger
from Messinger.Procesado_aerodinamico import Procesado_aerodinamico
from Messinger import fem_velocidades
import pickle
import matplotlib.pyplot as plt
import numpy as np
def procesado_termico(angulo_ataque,LWC,T_remanso,MVD,zona,V):
   global pto_remanso 
   import __main__

   x_experimental=__main__.x_experimental
   (Elementos,Velocidades,Superficie,Puntos)=Procesado_aerodinamico(V,angulo_ataque)
   Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie,Puntos)
   Modelo.set_T_remanso(273.15+T_remanso)
   Modelo.set_presion_remanso(1e5)
   betha=pickle.load(open('Eficiencias_coleccion//betha'+str(V)+'_'+str(angulo_ataque)+'_'+str(MVD)+'.p', "rb"))
   Modelo_termico = fem_velocidades.analisis_termico()
   Modelo_termico.calculo_S(Modelo)
   Modelo_termico.set_zona_perfil(zona)
   Modelo_termico.set_T_remanso(T_remanso+273.15)
   Modelo_termico.x_experimental = __main__.x_experimental
   Modelo_termico.T_experimental = __main__.T_experimental
   Modelo_termico.set_recovery_factor(1)
   Modelo_termico.set_presion_remanso(101325)
   Modelo_termico.set_LWC(LWC)

   Modelo_termico.set_diametro_caracteristico(0.02)
   Modelo_termico.set_velocidad_flujo(V)
   Modelo_termico.set_freezing_fraction(0.8)
   Modelo_termico.set_flujo_masico_entrada(0)
   Modelo_termico.set_T_superficie_anterior(273.15)
   Modelo_termico.set_cp_ws_anterior(1.004)
   Modelo_termico.set_T_superficie(273.15)
   Modelo_termico.set_local_collection_efficiency(0.5)
   Modelo_termico.set_freezing_fraction(1)
   Modelo_termico.set_tamano_gota(20e-6)
   Modelo_termico.V_e = Modelo_termico.V
   Modelo_termico.T_estatica = Modelo_termico.T_remanso -Modelo_termico.V**2/2/1004.5
   Modelo_termico.set_coeficiente_convectivo(400)
   Modelo_termico.calculo_todos_calores()
   Modelo_termico.Modelo_CFD = Modelo
   pto_remanso=Modelo_termico.x_remanso
   Modelo_termico.betha_nodos = betha
   h_l=[Modelo_termico.coeficiente_convectivo_laminar(Modelo_termico.longitud_equivalente(i)) for i in x_experimental]
   # h_l[0]=h_l[1]
   
   
   Modelo_termico.coeficiente_convectivo_valores = h_l
   (x,T_sur)=Modelo_termico.calculate_T_sur()
  
   return (np.array(x)+pto_remanso,T_sur)

def calcular_fraccion_congelacion(angulo_ataque,LWC,T_remanso,MVD,zona,V):
   global pto_remanso 
   import __main__
 
   x_experimental=__main__.x_experimental
   (Elementos,Velocidades,Superficie,Puntos)=Procesado_aerodinamico(V,angulo_ataque)
   Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie,Puntos)
   Modelo.set_T_remanso(273.15+T_remanso)
   Modelo.set_presion_remanso(1e5)
   betha=pickle.load(open('Eficiencias_coleccion//betha'+str(V)+'_'+str(angulo_ataque)+'_'+str(MVD)+'.p', "rb"))
   Modelo_termico = fem_velocidades.analisis_termico()
   Modelo_termico.calculo_S(Modelo)
   Modelo_termico.set_zona_perfil(zona)
   Modelo_termico.set_T_remanso(T_remanso+273.15)
   Modelo_termico.x_experimental = __main__.x_experimental
   Modelo_termico.T_experimental = __main__.T_experimental
   Modelo_termico.set_recovery_factor(1)
   Modelo_termico.set_presion_remanso(101325)
   Modelo_termico.set_LWC(LWC)

   Modelo_termico.set_diametro_caracteristico(0.02)
   Modelo_termico.set_velocidad_flujo(V)
   Modelo_termico.set_freezing_fraction(0.8)
   Modelo_termico.set_flujo_masico_entrada(0)
   Modelo_termico.set_T_superficie_anterior(273.15)
   Modelo_termico.set_cp_ws_anterior(1.004)
   Modelo_termico.set_T_superficie(273.15)
   Modelo_termico.set_local_collection_efficiency(0.5)
   Modelo_termico.set_freezing_fraction(1)
   Modelo_termico.set_tamano_gota(20e-6)
   Modelo_termico.V_e = Modelo_termico.V
   Modelo_termico.T_estatica = Modelo_termico.T_remanso -Modelo_termico.V**2/2/1004.5
   Modelo_termico.set_coeficiente_convectivo(400)
   Modelo_termico.calculo_todos_calores()
   Modelo_termico.Modelo_CFD = Modelo
   pto_remanso=Modelo_termico.x_remanso
   Modelo_termico.betha_nodos = betha
   h_l=[Modelo_termico.coeficiente_convectivo_laminar(Modelo_termico.longitud_equivalente(i)) for i in x_experimental]
   h_l[0]=h_l[1]
   
   Modelo_termico.coeficiente_convectivo_valores = h_l
   (x,T_sur)=Modelo_termico.calculate_T_sur()
   f = Modelo_termico.Fraccion_congelacion
   for i in range(len(f)):
      if f[i]>1 or f[i]<0:
         f[i]=np.nan
   return (np.array(x)+pto_remanso,Modelo_termico.Fraccion_congelacion)

def procesado_termico_rime(angulo_ataque,LWC,T_remanso,MVD,zona,V,x_sensor):
   global pto_remanso 
   import __main__
   if MVD<10:
      MVD=10
   else:
      for mvd in np.arange(10,120,10):
         if np.abs(mvd-MVD)<=5:
            MVD=mvd
            break
   x_experimental=__main__.x_experimental
   (Elementos,Velocidades,Superficie,Puntos)=Procesado_aerodinamico(V,angulo_ataque)
   Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie,Puntos)
   Modelo.set_T_remanso(273.15+T_remanso)
   Modelo.set_presion_remanso(1e5)
   betha=pickle.load(open('Eficiencias_coleccion//betha'+str(V)+'_'+str(angulo_ataque)+'_'+str(MVD)+'.p', "rb"))
   Modelo_termico = fem_velocidades.analisis_termico()
   Modelo_termico.calculo_S(Modelo)
   Modelo_termico.set_zona_perfil(zona)
   Modelo_termico.set_T_remanso(T_remanso+273.15)
   Modelo_termico.x_experimental = __main__.x_experimental
   Modelo_termico.T_experimental = __main__.T_experimental
   Modelo_termico.set_recovery_factor(1)
   Modelo_termico.set_presion_remanso(101325)
   Modelo_termico.set_LWC(LWC)

   Modelo_termico.set_diametro_caracteristico(0.02)
   Modelo_termico.set_velocidad_flujo(V)
   Modelo_termico.set_freezing_fraction(0.8)
   Modelo_termico.set_flujo_masico_entrada(0)
   Modelo_termico.set_T_superficie_anterior(273.15)
   Modelo_termico.set_cp_ws_anterior(1.004)
   Modelo_termico.set_T_superficie(273.15)
   Modelo_termico.set_local_collection_efficiency(0.5)
   Modelo_termico.set_freezing_fraction(1)
   Modelo_termico.set_tamano_gota(20e-6)
   Modelo_termico.V_e = Modelo_termico.V
   Modelo_termico.T_estatica = Modelo_termico.T_remanso -Modelo_termico.V**2/2/1004.5
   Modelo_termico.set_coeficiente_convectivo(400)
   Modelo_termico.calculo_todos_calores()
   Modelo_termico.Modelo_CFD = Modelo
   pto_remanso=Modelo_termico.x_remanso
   Modelo_termico.betha_nodos = betha
   Modelo_termico.x_surface = x_sensor
   h_l=[Modelo_termico.coeficiente_convectivo_laminar(Modelo_termico.longitud_equivalente(i)) for i in x_experimental]
   h_l[0]=h_l[1]
   
   Modelo_termico.coeficiente_convectivo_valores = h_l
   (x,T_sur)=Modelo_termico.calculate_T_sur_rime()
  
   return (np.array(x),T_sur)


def procesado_termico_borde_ataque(angulo_ataque,LWC,T_remanso,MVD,zona,V,Modelo,betha):
   global pto_remanso 
   import __main__
   x_experimental=__main__.x_experimental
   Modelo.set_T_remanso(273.15+T_remanso)
   Modelo.set_presion_remanso(1e5)
   Modelo_termico = fem_velocidades.analisis_termico()
   Modelo_termico.calculo_S(Modelo)
   Modelo_termico.set_zona_perfil(zona)
   Modelo_termico.set_T_remanso(T_remanso+273.15)
   Modelo_termico.x_experimental = __main__.x_experimental
   Modelo_termico.T_experimental = __main__.T_experimental
   Modelo_termico.set_recovery_factor(1)
   Modelo_termico.set_presion_remanso(101325)
   Modelo_termico.set_LWC(LWC)
   Modelo_termico.set_diametro_caracteristico(0.02)
   Modelo_termico.set_velocidad_flujo(V)
   Modelo_termico.set_freezing_fraction(0.8)
   Modelo_termico.set_flujo_masico_entrada(0)
   Modelo_termico.set_T_superficie_anterior(273.15)
   Modelo_termico.set_cp_ws_anterior(1.004)
   Modelo_termico.set_T_superficie(273.15)
   Modelo_termico.set_local_collection_efficiency(0.5)
   Modelo_termico.set_freezing_fraction(1)
   Modelo_termico.set_tamano_gota(20e-6)
   Modelo_termico.V_e = Modelo_termico.V
   Modelo_termico.T_estatica = Modelo_termico.T_remanso -Modelo_termico.V**2/2/1004.5
   Modelo_termico.set_coeficiente_convectivo(400)
   Modelo_termico.calculo_todos_calores()
   Modelo_termico.Modelo_CFD = Modelo
   pto_remanso=Modelo_termico.x_remanso
   Modelo_termico.betha_nodos = betha
   h_l=[Modelo_termico.coeficiente_convectivo_laminar(Modelo_termico.longitud_equivalente(i)) for i in x_experimental]
   #h_l[0]=h_l[1]
   Modelo_termico.coeficiente_convectivo_valores = h_l
   (x,T_sur)=Modelo_termico.calculate_T_sur_borde_ataque()
   print(pto_remanso)
   return (np.array(x)+pto_remanso,T_sur)


def procesado_termico_rapido(angulo_ataque,LWC,T_remanso,MVD,zona,V,Modelo,betha,HCT,stagnation,beta_0):
   global pto_remanso 
   import __main__
   # if stagnation:
   #    beta = [beta_0 for i in range(len(betha[:,0]))]
   #    betha[:,1]=beta
   beta = [beta_0 for i in range(len(betha[:,0]))]
   betha[:,1]=beta
   x_experimental=__main__.x_experimental
   Modelo.set_T_remanso(273.15+T_remanso)
   Modelo_termico = fem_velocidades.analisis_termico()
   Modelo_termico.calculo_S(Modelo)
   Modelo_termico.set_zona_perfil(zona)
   Modelo_termico.set_T_remanso(T_remanso+273.15)
   Modelo_termico.x_experimental = __main__.x_experimental
   Modelo_termico.T_experimental = __main__.T_experimental
   Modelo_termico.set_recovery_factor(1)
   Modelo_termico.set_presion_remanso(101325)
   Modelo_termico.set_LWC(LWC)

   Modelo_termico.set_diametro_caracteristico(0.02)
   Modelo_termico.set_velocidad_flujo(V)
   Modelo_termico.set_freezing_fraction(0.8)
   Modelo_termico.set_flujo_masico_entrada(0)
   Modelo_termico.set_T_superficie_anterior(273.15)
   Modelo_termico.set_cp_ws_anterior(1.004)
   Modelo_termico.set_T_superficie(273.15)
   Modelo_termico.set_local_collection_efficiency(0.5)
   Modelo_termico.set_freezing_fraction(1)
   Modelo_termico.set_tamano_gota(20e-6)
   Modelo_termico.V_e = Modelo_termico.V
   Modelo_termico.T_estatica = Modelo_termico.T_remanso -Modelo_termico.V**2/2/1004.5
   Modelo_termico.set_coeficiente_convectivo(400)
   Modelo_termico.calculo_todos_calores()
   Modelo_termico.Modelo_CFD = Modelo
   pto_remanso=Modelo_termico.x_remanso
   Modelo_termico.betha_nodos = betha
   if stagnation:h_l=[HCT for i in x_experimental]
   else:
      h_l=[Modelo_termico.coeficiente_convectivo_laminar(Modelo_termico.longitud_equivalente(i)) for i in x_experimental]
   # h_l=[Modelo_termico.coeficiente_convectivo_laminar(Modelo_termico.longitud_equivalente(i)) for i in x_experimental] 
   Modelo_termico.coeficiente_convectivo_valores = h_l
   # plt.figure()
   # plt.plot(x_experimental,h_l)
   # plt.show()
   (x,T_sur)=Modelo_termico.calculate_T_sur()
   return (np.array(x)+pto_remanso,T_sur)