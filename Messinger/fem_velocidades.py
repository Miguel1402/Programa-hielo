import numpy as np
from shapely.geometry import Point,LineString
from shapely.geometry.polygon import Polygon
import pickle
import pandas as pd
from scipy import linalg
class modelo_fem(object):
    """[Hace un modelo Fem para integrar las trayectorias de las gotas]

    Args:
        object ([object]): []
    """    

    def __init__(self,Elementos,Velocidades,Superficie,points):

        self.x_nodo = np.array(points['x'])
        self.y_nodo = np.array(points['y'])
        self.Points = points
        self.u_nodo = np.array(Velocidades['u'])
        self.v_nodo = np.array(Velocidades['v'])
        self.puntos_singulares=[]
        self.Velocidades=Velocidades
        self.procesado_elementos(Elementos)
        self.x_superficie = np.array(Superficie['x'])
        self.y_superficie = np.array(Superficie['y'])
        self.rho_w = 1000 
        self.Superficie = Superficie
        self.perfil()
    
    def procesado_elementos(self,Elementos):
        """[Anade en elementos la x e y maxima y minima. 
        Asi luego se puede acceder de manera mas directa al elemento]

        Args:
            Elementos ([DataFrame]): [Los puntos que pertenecen a cada elemento]
        """        
        x_max =[]
        x_min=[]
        y_max =[]
        y_min=[]
        u=[]
        for indice in Elementos.index:
            element=Elementos.loc[indice]
            if str(element['punto_4'])!='nan':
                x = [self.Points['x'].loc[int(element['punto_1'])],
                    self.Points['x'].loc[int(element['punto_2'])],
                    self.Points['x'].loc[int(element['punto_3'])],
                    self.Points['x'].loc[int(element['punto_4'])]]
                y = [self.Points['y'].loc[int(element['punto_1'])],
                    self.Points['y'].loc[int(element['punto_2'])],
                    self.Points['y'].loc[int(element['punto_3'])],
                    self.Points['y'].loc[int(element['punto_4'])]]
            else:
                x = [self.Points['x'].loc[int(element['punto_1'])],
                            self.Points['x'].loc[int(element['punto_2'])],
                            self.Points['x'].loc[int(element['punto_3'])]]
                y = [self.Points['y'].loc[int(element['punto_1'])],
                    self.Points['y'].loc[int(element['punto_2'])],
                    self.Points['y'].loc[int(element['punto_3'])]]
            x_max.append(max(x))
            x_min.append(min(x))
            y_max.append(max(y))
            y_min.append(min(y))
           
        Elementos['x_max']=x_max
        Elementos['x_min']=x_min
        Elementos['y_max']=y_max
        Elementos['y_min']=y_min
     

        self.Elementos = Elementos

    def perfil(self):
        """[Cuando se llama a esta función se anade en el objeto la normal asociada a la superficie. También la velocidad exterior que puede valer para calcular el coeficiente convenctivo]
        """  
        distancia =0.0005      
        self.puntos_v_inf=[]
        x = np.array(self.Superficie['x'])
        y = np.array(self.Superficie['y'])
        coordenadas = []
        for i in range(len(x)):
            coordenadas.append((x[i],y[i]))
        self.polygon = Polygon(coordenadas)
        coordenadas = np.array(self.polygon.exterior.coords)
        self.exterior = LineString(coordenadas)
        normal = np.zeros((len(x),2))
        for i in range(1,len(x)):
            mod = np.sqrt((x[i]-x[i-1])**2+(y[i]-y[i-1])**2)
            if mod > 0:
                if y[i]<0:
                    
                    normal[i,1] = -(x[i]-x[i-1])/mod
                    normal[i,0] = (y[i]-y[i-1])/mod
                else:
                    mod = np.sqrt((x[i]-x[i-1])**2+(y[i]-y[i-1])**2)
                    normal[i,1] = (x[i]-x[i-1])/mod
                    normal[i,0] = -(y[i]-y[i-1])/mod

        self.normal=normal
        V_infinito = []
        for i in range(len(x)):
            if (normal[i,0]!=0) or (normal[i,1]!=0):
                try:
                    V=self.Velocidad(x[i]-normal[i,0]*distancia,y[i]-normal[i,1]*distancia)
                    V = np.sqrt(V[0]**2+V[1]**2)
                    if y[i]<0:
                        V_infinito.append([-x[i],V])
                    else:
                        V_infinito.append([x[i],V])

                    self.puntos_v_inf.append([x[i]-normal[i,0]*distancia,y[i]-normal[i,1]*distancia])
                except:
                    V=self.Velocidad(x[i]+normal[i,0]*distancia,y[i]+normal[i,1]*distancia)
                    V = np.sqrt(V[0]**2+V[1]**2)
                    if y[i]<0:
                        V_infinito.append([-x[i],V])
                    else:
                        V_infinito.append([x[i],V])
                    self.puntos_v_inf.append([x[i]+normal[i,0]*distancia,y[i]+normal[i,1]*distancia])
               
                
        self.V_infinito = np.array(V_infinito)


   
                
    def Velocidad(self,x_punto,y_punto):

        """[Calcula la velocidad en las coordenadas x,y]



        Arguments:

            x {[float]} -- [coordenada x (m)]

            y {[float]} -- [coordenada y (m)]



        Returns:

            [tuple] -- [(u,v) del flujo]

        """        

        Elementos=self.Elementos
        Elementos_filtrados=Elementos[(Elementos.x_max>=x_punto)&(Elementos.y_max>=y_punto)&(Elementos.x_min<=x_punto)&(Elementos.y_min<=y_punto)]       
        for indice in Elementos_filtrados.index:
            element=self.Elementos.loc[indice]
            if str(element['punto_4'])=='nan':
                x = [self.Points['x'].loc[int(element['punto_1'])],
                    self.Points['x'].loc[int(element['punto_2'])],
                    self.Points['x'].loc[int(element['punto_3'])]]
                y = [self.Points['y'].loc[int(element['punto_1'])],
                    self.Points['y'].loc[int(element['punto_2'])],
                    self.Points['y'].loc[int(element['punto_3'])]]
                coordenadas =[(x[i],y[i]) for i in range(3)]
                poly=Polygon(coordenadas)
                p1 = Point(x_punto,y_punto)
                if p1.within(poly):
                    celda = int(element['celdas'])
                    U = self.Velocidades['u'].loc[celda]
                    V=self.Velocidades['v'].loc[celda]
        return (U,V)


    def stream_line(self,x_0,y_0,x_f,n):
        """[Calcula las líneas de corriente del flujo]

        Arguments:
            x_0 {[float]} -- [coordenada x inicial en metros]
            y_0 {[float]} -- [coordenada y inicial en metros]
            x_f {[float]} -- [coordenada x final en metros]
            n {[int]} -- [numero de nodos]

        Returns:  (x,y)
            [tuple] -- [(x,y) del streamline]
        """        
        x = np.linspace(x_0,x_f,n)
        y = np.zeros(n)
        y[0] = y_0
        for i in range(1,n):
            try:
                (u,v) = self.Velocidad(x[i-1],y[i-1])
            except:
                x=x[:i]
                y=y[:i]
                break
            
            if u!=0:
                y[i] = y[i-1] + v/u*(x[i]-x[i-1])
            else:   y[1] = y[0]
        return (x,y)

    def viscosity(self,T):
        """[viscosidad en función de la temperatura]

        Arguments:
            T {[float]} -- [Temperatura en K]

        Returns:
            [mu] -- [viscosidad en el SI]
        """        
        return 5e-8*T+3.46e-5

    def set_T_remanso(self,T_remanso):
        """[Fija la temperatura de remanso del problema]

        Arguments:
            T_remanso {[float]} -- [Temperatura (K)]
        """        
        self.T_remanso = T_remanso
    
    def Temperatura(self,u,v):
        """[Temperatura estatica]

        Arguments:
            u {[float]} -- [velocidad horizontal (m/s)]
            v {[float]} -- [velocidad horizontal (m/s)]

        Returns:
            [float] -- [Temperatura estatica (K)]
        """        
        V = np.sqrt(u**2+v**2)
        return self.T_remanso - V**2/(2*1004.5)

    def set_presion_remanso(self,P_0):
        self.P_0 = P_0

    def presion_estatica(self,T,V):
        M =V/np.sqrt(1.4*287*T)
        gamma = 1.4
        return self.P_0/(1+(gamma -1)/2*M**2)**(gamma/(gamma-1))
    
    def densidad(self,T,P):
        return P/(287*T)

    def Reynolds(self,V,rho_a,mu_a,D):
        """[Devuelve el número de Reynolds, a partir del diámetro característico del perfil]

        Args:
            V ([float]): [Velocidad aerodinámica m/s^2]
            rho_a ([float]): [Densidad estática del aire (Kg/m^3)]
            mu_a ([float]): [Viscosidad el aire (N/m^2 s^-1)]
            D ([float]): [Diámetro característico del perfil[m]]

        Returns:
            [float]: [Número de Reynolds del flujo]
        """        
        return rho_a*V*D/mu_a
    
    def CD(self,Re):
        """[Coeficiente de resistencia  en función del número de Reynolds]

        Args:
            Re ([Float]): [número de Reynolds]

        Returns:
            [float]: [Coeficiente de resistencia]
        """        
        if (Re<=1000)and (Re >1): 
            c=  24./Re+6./(1+np.sqrt(Re))+.27
        elif ((Re<=3600.) and (Re>1000)):
            c=0.6649+0.2712e-3*Re+1.22e-7*Re**2.-10.919e-12*Re**3.  
        elif (Re <=1): c=0
        return c 

    def trayectoria_gota(self,x_0,x_f,y_0,U_d0,V_d0,D,n_nodos):
        """[Calcula la trayectoria de la gota desde un punto inicial hasta otra final]

        Args:
            x_0 ([float]): [X inicial en la que se empieza a integrar la trayectoria (m)]
            x_f ([float]): [X final en la que se empieza a integrar la trayectoria (m)]
            y_0 ([float]): [y inicial en la que se empieza a integrar la trayectoria (m)]
            U_d0 ([float]): [Velocidad inicial de la gota en la dirección x (m/s)]
            V_d0 ([float]): [Velocidad inicial de la gota en la dirección y (m/s)]
            D ([float]): [diámetro característico (m)]
            n_nodos ([int]): [número de nodos de integración ]

        Returns:
            [tuple]: [(x,y,U_d,V_d,t)]
                    x ([list]): [lista de x de la trayectoria]
                    y ([list]): [lista de y de la trayectoria]
                    U_d ([list]): [Velocidades de las trayectorias horizontales (m/s)]
                    V_d ([list]): [Velocidades de las trayectorias verticales (m/s)]
                    t ([list]): [lista de tiempos de la trayectoria]
                    
        """        
        #se inician las variables de la trayectoria
        y = [y_0]
        U_d = [U_d0]
        V_d = [V_d0]
        t=[0]
        delta_x=(x_f-x_0)/(n_nodos-1.)
        x=[x_0]
        #Se importa el polígono que define el perfil
        poly = self.polygon
        #Se itera de uno al número de nodos que se desea
        for i in range(1,n_nodos):
            #Se le añade a la nueva x una delta
            x.append(x[-1]+delta_x)
            #Se llama a la función Velocidad, obteniendo las componentes horizontal y vertical 
            #de la velocidad el flujo
            (u_a,v_a) = self.Velocidad(x[i-1],y[i-1])
            #Se calculan las variables aerodinámicas pertinentes
            V_a = np.sqrt(u_a**2 + v_a**2)
            T_a=self.Temperatura(u_a,v_a)
            P_a = self.presion_estatica(T_a,V_a)
            mu_a=self.viscosity(T_a)
            rho_a = self.densidad(T_a,P_a)
            #Se define la componente vertical de las fuerzas inerciales  (ver Zarling)
            fg=9.81*(rho_a/self.rho_w-1)
            #Se define la velocidad relativa entre el fluido y la gota
            uf_x = u_a - U_d[i-1]
            uf_y = v_a - V_d[i-1]
            V_f = np.sqrt(uf_x**2 + uf_y**2)
            Re = self.Reynolds(V_a,rho_a,mu_a,D)
            #Se calcula la fuerza de resistencia ejercida por el aire sobre la gota de agua 
            f_d=3./4.*self.CD(Re)/D*rho_a/self.rho_w*V_f    
            #u_x(i)=delta_x/u_x(i-1)*f_d*uf_x+u_x(i-1)
            #Se calcula el siguiente paso de la velocidad absoluta de la gota.
            #Se realiza para la integración un Euler en velocidad.
            #Para realizar la integracion se apliga la regla de la cadena

            U_d.append(delta_x/U_d[i-1]*f_d*uf_x+U_d[i-1])
            #u_y(i)=delta_x/u_x(i-1)*(f_d*uf_y+fg)+u_y(i-1)
            V_d.append(delta_x/U_d[i-1]*(f_d*uf_y+fg)+V_d[i-1])
            t.append((x[i]-x[i-1])/U_d[i-1]+t[i-1])
            y.append(V_d[i-1]*(t[i]-t[i-1])+y[i-1])
            # se define un objeto Point  de shapely
            p1 = Point(x[-1],y[-1])
            if p1.within(poly):
                #Se define una linea con el punto anterior
                linea = LineString([(x[-2],y[-2]),(x[-1],y[-1])])
                #Se mira si intesecciona con el exterior del perfil
                Interseccion=self.exterior.intersection(linea)
                try:
                    #Se define el ultimo punto como la intersección
                    x[-1] = Interseccion.x
                    y[-1] = Interseccion.y
                except:
                    
                    x[-1]=Interseccion[0].x
                    y[-1]=Interseccion[0].y
                break
        return (x,y,U_d,V_d,t)

    def proyeccion_gota(self,x_0,y_0,x_f,n_nodos):
        x= np.linspace(x_0,x_f,n_nodos)
        poly = self.polygon
        for i in range(n_nodos):
            p1 = Point(x[i],y_0)
            if p1.within(poly):
                break
        return x[i]


    def punto_remanso(self):
        """[Calcula el punto de remanso y lo guarda en la variable del objeto x_remanso]

        Args:
            Modelo_aerodinamico ([object]): [Modelo CFD aerodinámico resuelto con OpenFoam]
        """        
        V_infinito = pd.DataFrame(self.V_infinito,columns=['x','V'])
        self.x_remanso = float(V_infinito[V_infinito.V==min(self.V_infinito[:,1])]['x'])
       

    def calculo_S(self):
        """[Guarda en el objeto un array S donde la primera columna es la s la segunda el valor x y la tercera la y]

        Args:
            Modelo ([object]): [Objeto del CFD resulto con OpenFoam]
        """        
        self.punto_remanso()
        x_0 = self.x_remanso
       
        x = list(self.Superficie.x)
        y = list(self.Superficie.y)
        S=[]
        for i in range(len(self.Superficie.index)):
            
            if x_0<0:
                if x[i+1]<=abs(x_0) and x[i]>=abs(x_0):
                    S.append([0,x[i],y[i]])
                    n=i
                    break
            else:
                if x[i+1]>=x_0 and x[i]<=x_0:
                    S.append([0,x[i],y[i]])
                    n=i
                    break
            
        for i in range(n+1,len(self.Superficie.index)):
            S.append([0,x[i],y[i]])
            s=0
            for j in range(len(S)-1):
                s=s+np.sqrt((x[j+1]-x[j])**2+(y[j+1]-y[j])**2)
            S[-1][0]=s
        self.S=np.array(S)

class analisis_termico(object):
    def __init__(self):       
        """[Se comienza la función definiendo las constantes]
        """        
        self.cp_a = 0.24 #cal/g K
        self.tabla_presiones = pickle.load( open( "presiones_vapor.pickle", "rb" ) )
        self.t_f = 273.15 #punto de congelacion del agua
        self.tabla_hielo = pickle.load( open( "propiedades_hielo.pickle", "rb" ) )
        self.T_remanso=None
        self.T_superficie= None
        self.recovery_factor=None 
        self.p_remanso= None
        self.LWC=None
        self.cuerda= None
        self.d=None
        self.V= None
        self.n_0=None
        self.betha=None
        self.m_in=None
        self.T_s_anterior=None
        self.cp_ws_anterior=None

    def set_T_remanso(self,T_0):
        """[Fija la temperatura de remanso del objeto]

        Args:
            T_0 ([float]): [Temperatura de remanso (K)]
        """        
        self.T_remanso=T_0

    def set_recovery_factor(self,r):
        """[Fija el factor de recuperacion ]

        Args:
            r ([float]): [Factor de recuperacion]
        """        
        self.recovery_factor=r 

    def set_T_superficie(self,T_sur):
        """[Fija las temperaturas de la superficie]

        Args:
            T_sur ([array]): [Temperaturas de la superficie a lo largo de la cuerda]
        """        
        self.T_superficie=T_sur

    def set_presion_remanso(self,P_0):
        """[Fija la presión de remanso exterior]

        Args:
            P_0 ([float]): [presión de remanso (Pa)]
        """        
        self.p_remanso = P_0

    def set_LWC(self,LWC):
        """[Fija el contenido de agua liquida en g/m3]

        Args:
            LWC ([float]): [Contenido de agua liquida de la nube]
        """        
        self.LWC=LWC

    def set_cuerda(self,c):
        """[Fija un valor de la cuerda del perfil]

        Args:
            c ([float]): [cuerda en m]
        """        
        self.cuerda = c

    def set_diametro_caracteristico(self,d):
        self.d =d

    def set_velocidad_flujo(self,V):
        """[Fija la velocidad de flujo]

        Args:
            V ([float]): [Velocidad de flujo en el infinito en m/s]
        """        
        self.V=V

    def set_coeficiente_convectivo(self,h_c):
        """[Fija una lista de valores del coeficiente convectivo en los sitios donde se han tomado datos experimentalmente]

        Args:
            h_c ([float]): [coeficiente de calor convectivo cal/(K m^2)]
        """        
        self.h_c = h_c

    def set_freezing_fraction(self,n):
        self.n_0 = n

    def set_local_collection_efficiency(self,betha):
        """[Fija la eficiencia de colección local en cada punto de la cuerda del perfil]

        Args:
            betha ([float]): [array]
        """        
        self.betha = betha

    def set_flujo_masico_entrada(self,m_in):
        """[Fija el flujo másico de entrada para un determinado nodo]

        Args:
            m_in ([float]): [Gasto másico de agua que entra en el volumen de control (g/m^3)]
        """        
        self.m_in =m_in

    def set_T_superficie_anterior(self,T_sur_i_menos1):
        """[Fija la variable de la temperatura anterior]

        Args:
            T_sur_i_menos1 ([float]): [Temperatura de la superficie en el nodo anterior]
        """        
        self.T_s_anterior=T_sur_i_menos1

    def set_cp_ws_anterior(self,cp_ws_i_menos1):
        """[Fija la variable del calor específico en la superficie anterior]

        Args:
            cp_ws_i_menos1 ([float]): [calor específico en la superficie anterior]
        """        
        self.cp_ws_anterior = cp_ws_i_menos1

    def set_V_exterior(self,Ve):
        """[Fija el valor de la velocidad de flujo exterior en un determinado nodo]

        Args:
            Ve ([float]): [Valor de la velocidad exterior en m/s]
        """        
        self.V_e
    
    def __str__(self):
        texto ='Copyright INTA'
        texto = texto+'\nArea de Materiales compuestos, Departamento de Materiales 2020'
        texto = texto+'\nDesarrollador Miguel Gonzalez del Val'
        texto = texto+'\nReferencias:'
        texto = texto+'\n[1] Zarling 1981 Heat and mass transfer from freely falling drops at low temperatures'
        texto = texto+'\n[2] NASA, Manual of Scaling Methods'
        texto = texto+'\n[3] BERNARD L. MESSINGER, Equilibrium Temperature of an Unheated Icing Surface as a Function of Air Speed'
        return texto

    def set_calor_conductivo(self,calor_conductivo):
	    self.calor_conductivo = calor_conductivo
    

    def calor_convectivo(self,h_c,T_superficie,T_estatica,V):
        '''
        Calcula el calor convectivo en la linea de remanso (manual
        scaling methods eq. 3.53)
        variables:
            h_c: film coefficient (cal/(hr m^2 K)
            T_superficie: K
            T_estatica: K
            V:Velocidad del flujo
            cp_a: coeficiente calorifico a presion cte
            q_h (cal/hr m^2)
        '''
        
        return h_c*(T_superficie-T_estatica-V**2/(2*1004.5))/3600
    def film_coefficient(self,k_a,Nu_a,d):
        '''
        Calcula el coeficiente convectivo
        Variables:
            ka: coeficiente de conductivo del aire (cal/(m K hr)
            Nu_a: Nº de Nusselt del aire
            d: diametro caracteristico del perfil (borde de ataque) (metros)
        '''
        return k_a*Nu_a/d

    def Nusselt(self,Pr,Re):
        '''
        Ecuación (3.33) Manual de Scaling Methods (NASA)
        Pr: Nº de Prandtl
        Re: Nº De Reynolds
        '''
        return 1.14*Pr**0.4*Re**0.5
    def Pr(self,cp_a,mu_a,k_a):
        '''
        Numero de Prandtl (ecuacion 3.34)
        c_p: coeficiente calorifico a presion cte 
        mu_a  viscosidad
        k_a: conductividad termica
        '''
        return cp_a*mu_a/k_a

    def Reynolds(self,V,c,rho,mu_a):
        '''
        Numero de Reynolds (ecuacion 3.39)
        Se hace respecto a la cuerda
        rho_a: densidad
        c: cuerda del perfil
        mu_a viscosidad del aire
        V_velocidad del flujo
        '''
        return V*c*rho/mu_a
    def thermal_conductivity(self,T_film):
        '''
        Conductividad térmica del aire. 
        es una funcion de la temperatura (ver ecuacion A.3 de Manual Scaling Methods NASA)
        '''
        k_a =-12.69 + 2.029 *np.sqrt(T_film) # cal/(hr m K)
        return k_a

    def Temperatura_film(self,T_st,T_s):
        '''
        Ecuacion 3.35 de Manual Scaling methods
        '''
        return 0.5*(T_s+T_st)
    def convective_mass_coeff(self,h_c,Sc,Pr):
        '''
        Calcula el coeficiente de transmisión de masa
        Ecuacion A35 de manual scaling methods
        '''
        
        return h_c/self.cp_a*(Pr/Sc)**.67

    def difusividad(self,T_film,p_st):
        '''
        Calcula la difusividad en funcion de la temperatura. 
        Ecuacion A.4 Manual de Scaling Methods
        
        '''
        
        return 0.211*(T_film/273.15)**1.94*(101325/p_st)

    def Schmidt(self,mu_a,rho_a,Difusividad):
        '''
        Numero de Schmidt (Sc)
        Difusividad del aire, mu_a viscosidad y rho_a es la densidad
        Ecuacion  3.48 Manual Scaling NASA
        '''
        return mu_a/(rho_a*Difusividad)


    def presion_vapor(self,T):
        '''
        Presion de vapor del agua en funcion del ratio de humedad
        Datos de Beltramino, G., Rosso, L., Cuccaro, R., Tabandeh, S., Smorgon, D., & Fernicola, V. (2019).
        Accurate vapour pressure measurements of super cooled water in the temperature range between 252 K and 273 K. 
        The Journal of Chemical Thermodynamics, 105944. doi:10.1016/j.jct.2019.105944 
        
        '''
        Temperatura =self.tabla_presiones[:,0]
        P_vapor = self.tabla_presiones[:,2]
        p = P_vapor[-1]
        if T > Temperatura[-1]: p=P_vapor[-1]
        for i in range(1,len(Temperatura)):
            if (T>=Temperatura[i-1]) and (T<=Temperatura[i]):
                p = P_vapor[i-1]+ (P_vapor[i]-P_vapor[i-1])/(Temperatura[i]-Temperatura[i-1])*(T-Temperatura[i-1])
                break
        return p

    def mass_water_evaporation(self,p_ww,p_w,h_g,p_estatica,T_surface,T_remanso,p_remanso,n):
        '''
        p_ww,p_w,h_g,p_estatica,T_estatica,T_remanso,p_remanso,0.5
        Calcula la masa de aire que se evapora (g/(hr m^2)).
        Ecuacion 3.46 de Manual Scaling
        h_g constante de transmision de masa por evaporacion (#g/(m^2 hr) )
        p_ww Presion parcial de vapor del agua en la superficie
        p_w Presion parcial de vapor del agua en el flujo externo
        p_st presion estatica
        n: freezing factor
        Messiger  dice que si n=1 no hay evaporacion (m_e=0). 
        '''
       
        m = h_g*((p_ww/T_surface-p_remanso/T_remanso*p_w/p_estatica)/(p_remanso/T_remanso/.622-p_ww/T_surface))
        if n>=1 or m<0:m=0

        return m
    def heat_evaporation(self,m_e,Lambda_v):
        '''
        Calcula el flujo de calor por evaporacion (cal/(m^2 hr)).
        m_e : la masa de aire que se evapora (g/(hr m^2))
        Lambda_v = calor latente de evaporacion  (cal/g)
        Se toma como incompresible
        '''
        return m_e*Lambda_v

    def viscosity(self,T):
        '''
        #g/(cm s)
        Viscosidad del aire en funcion de la temperatura
        '''    
        return 10**(-4)/(.12764+124.38/T)#g/(cm s)

    def Latent_heat_vaporisation(self,T):
            '''
            calor latente de vaporizacion (cal/g) en funcion de la Temperatura
            de la superficie (glaze T=273.15 K) 
            T en Kelvin
            '''
            E =  0.197 + 3.670e-4*T
            return 597.3*(273.15/T)**E
        

    def cp_ws(self,T_superficie):
        '''
        Calcula el calor especifico del agua en la superficie 
        unidades cal/(g K)
        cp = 8.29e-5*(T_superficie.273.15)**2
        '''
        cp = 1.0074+8.29e-5*(T_superficie-273.15)**2
        return cp

    def impinging_mass_water(self,LWC,V,betha_0):
        '''
        impinging_mass_water Masa de agua que choca
        Ecuacion 3.50 Manual scaling methods
        LWC: liquid water content
        V: Velocidad del flujo
        betha_0 = catch efficiency at stagnation (eficiencia de adhesion en remanso)
        '''
        return LWC*V*betha_0

    def modified_inertia_parameter(self,LAMBDA_LAMBDAstokes,K):
        '''
        Langmuir and Blodgett’s expression for modified inertia
        parameter (eq. 3.8 Scaling methods)
        '''
        if K>1/8:return 1/8+ LAMBDA_LAMBDAstokes*(K-1/8)
        else:print('K no valida')
    def inertia_parameter(self,rho_w,delta,V,d,mu_a):
        '''
        K=rho_w*delta**2*V/(18*d*mu_a)
        K is the non-dimensional inertia parameter
        defined by Langmuir and Blodgett  Eq. 3.5
        rho_w, densidad del agua
        delta: diametro gota
        d: diametro curvatura del borde de ataque del perfil
        mu_a: viscosidad
        '''
        return rho_w*delta**2*V/(18*d*mu_a)
        
    def dimensionless_range_parameter(self,V,delta,rho,mu_a):
        '''
        define el parametro lambda/lambda_stokes (equation 3.9 Manual scaling methods)
        delta (m): tamaño de la gota
        V (m/s): velocidad
        mu_a (g/m s) viscosidad
        rho (g/m^3)
        '''
        Re = self.Reynolds(V,delta,rho,mu_a)
        return 1/(.8388+0.001483 *Re +.1847*Re**0.5)

    def catch_efficiency_stagnation(self,K_0):
        '''
        Ecuacion 3.13 del Manual of Scaling Methods
        (betha_0)
        K_0 es el parametro de inercia (ver inertia_parameter)
        '''
        return 1.4*(K_0-1/8)**.84/(1.4*(K_0-1/8)**.84+1)

    def Sensible_Heat_Water(self,cp_ws,dot_m,t_st):
        ''' 
        Ecuacion 5 del apartado 3.5 Manual Scaling methods
        cp,ws: calor especifico del agua a la temperatura de la superficie,
        dot_m: flujo masico del agua que impregna la superficie
        t_st temperatura estatica
        '''
        return dot_m*cp_ws*(self.t_f-t_st)


    def heat_kinetic(self,V,dm_dt):
        '''
        calor debido a Energia cinetica del flujo (cal/(m^2 s))S
        ver termino (10) del apartado 3.5
        dm_dt (kg/(s m^2):impinging_mass_water
        '''
        return dm_dt*V**2/2/4.1868 #cal/m^2 s

    def Temperatura_pared(self,r,V,T_inf):
        '''
        FUente: http://www.thermopedia.com/content/291/
        Se supone flujo laminar
            Pr:Numero de Prandtl
            T_inf: temperatura estatica del flujo exterior
            V Velocidad del flujo
             For some simple cases, its value can be estimated as follows: 
                At the front stagnation point of bodies in the flow, r = 1;
                in a laminar boundary layer on a plane plate, r =Pr^0.5 for Prandtl numbers 0.5 < Pr <10; 
                in a turbulent boundary layer on a plate, r =Pr^0.33 for Prandtl numbers close to 1
        '''
        return T_inf+r*V**2/(2*self.cp_a)
    
    def latent_heat_freezing(self,T_superficie):
        '''
        Ecuacion A19 Manual Scaling Methods
        Unidades: cal/g
        T_superficie (K): en glaze es 0 K
        79.7 + .485*(T_superficie-273.15)-2.5e-3*(T_superficie-273.15)**2 #cal/g
        '''
        return 79.7 + .485*(T_superficie-273.15)-2.5e-3*(T_superficie-273.15)**2 #cal/g
    def relative_heat_factor(self,LWC,cp_ws,betha_0,V,h_c):
        '''
        Ecuacion 3.55 del Manual Scaling Methods
        LWC (g/m^3)
        V (m/s)
        cp_ws (calor especifico del agua en la superficie) (cal/g)
        h_c: cal/(s m^2) coeficiente de transmision de calor convectivo
        b=LWC*V*betha_0*cp_ws/h_c
        '''
        return LWC*V*betha_0*cp_ws/h_c

    def drop_energy_transfer(self,T_estatica,V,cp_ws):
        '''
        ecuacion 3.56
        phi = T_f - T_estatica - V**2/(2*cp_ws)
        V(m/s)
        Temperaturas en K
        cp_ws (calor especigico del agua en la superficie J/(kg K))
        '''
        return self.t_f - T_estatica - V**2/(2*cp_ws)
        
    def air_energy_transfer(self,T_superficie,T_remanso,V,h_g,h_c,p_ww,p_w,p_estatica,p_remanso,T_estatica,Lambda_v):
        '''
        Ecuacion A52 
        Air energy transfer
        pww: Pa presion de vapor del agua en la superficie
        P_st: Pa presion estatica
        p_w: Pa p vapor del agua en el flujo
        p_remanso: Pa
        T_remanso:
        T_st
        T_superficie
        P_remanso
        h_g:coeficiente de transmision de masa por evaporacion (cal/(hr m^2 hr))
        h_c:coeficiente de transmision de calor por conveccion (cal/(hr m^2 hr))
        cp_a : J/Kg K
        Lambda_v: calor por evaporacion (cal/g)
        '''
        termino_1 = T_superficie - T_estatica-V**2/(2*1004.5) #K
        termino_2 = h_g/h_c*Lambda_v*((p_ww/T_estatica-p_remanso/T_remanso*p_w/p_estatica)/(
            p_remanso/T_remanso/.622-p_ww/T_estatica))
        return termino_1+termino_2

    def freezing_fraction_stagnation(self,cp_ws,Lamdbda_f,phi,theta,b):
        '''
        ecuacion 3.59
        n_0 es un numero adimensional que describe la cantidad de hielo que congela al llegar a la superficie
        Lambda_f es el calor latente de congelacion (cal/g)
        cp_ws calor especifico del agua en la superficie (cal/(K g))
        phi: drop_energy_transfer en Kelvin
        theta: air_energy_transfer en Kelvin
        b:relative_heat_factor adimensional
        '''
        return cp_ws/Lamdbda_f*(phi+theta/b)


    def calor_flujo_agua(self,n,dot_m,dot_me,cp_ws,T_superficie):
        '''
        Flujo de calor debido al flujo de agua que se escapa del volumen de control 
        unidades: cal/(m^2 s)
        Ecuacion 6 del apartado 3.5
        dot_m flujo de agua que impregna la superficie (g/(s m^2))
        dot_me flujo de agua que se evapora (g/(s m^2))
        cp_we: calor especifico del agua en la superficie (cal/(g K))
        n : ratio de congelacion
        T_superficie: K
        '''
        return ((1-n)*dot_m-dot_me)*cp_ws*(self.t_f-T_superficie)

    def cp_is(self,T_superficie):
        '''
        calor especifico del hielo (cal/g K)
        '''
        Temperatura =self.tabla_hielo[:,0]
        Cp = self.tabla_hielo[:,3]
        T=T_superficie-273.15
        cp = Cp[0]
        if T >0: cp = Cp[0]
        for i in range(1,len(Temperatura)):
            if (T>=Temperatura[i]) and (T<=Temperatura[i-1]):
                cp = Cp[i-1]+ (Cp[i]-Cp[i-1])/(Temperatura[i]-Temperatura[i-1])*(T-Temperatura[i-1])
                #cp = Cp[i-1]
                break
        return cp/4.1848

    def conductividad_hielo(self,T_superficie):
        '''
        calor especifico del hielo (cal/m s K)
        '''
        
        Temperatura =self.tabla_hielo[:,0]
        K = self.tabla_hielo[:,2]
        T=(T_superficie)/2-273.15
        k= K[0]
        for i in range(1,len(Temperatura)):
            if (T>=Temperatura[i]) and (T<=Temperatura[i-1]):
                k = K[i-1]+ (K[i]-K[i-1])/(Temperatura[i]-Temperatura[i-1])*(T-Temperatura[i-1])
                #cp = Cp[i-1]
                break
        return k/4.1848
    
    def heat_from_ice(self,dot_m,n,T_superficie,cp_is):
        '''
        flujo de calor del hielo
        q_i = dot_m*n*cp_is*(self.t_f-T_superficie)
        q_i: unidades 
        cal/(m^2 s)
        dot_m (g/s m^2)
        n adimensional
        T_superficie K
        cp_is:cal/(g K)
        '''
        return dot_m*cp_is*n*(self.t_f-T_superficie)

    def heat_freezing (self,dot_m,Lambda_f,n):
        '''
        ecuacion 8 de la seccion 3.5 manual scaling methods
        cal/(m^2 s)
        '''
        return n*dot_m*Lambda_f

    def heat_conduction(self,Delta,chi,l,k_i,T_surface,r,V,T_estatica):
        """[Calor perdido por el modelo por conducción]

        Arguments:
            Delta {[Float]} -- [Espesor de la capa de hielo (m)]
            chi {[Float]} -- [arco característico de la superficie de control (m)]
            l {[Float]} -- [longitud característica de la superficie de control (m)]
            k_i {[Float]} -- [Conductividad térmica hielo (cal /(m K s)]
            T_surface {[Float]} -- [Temperatura de la superficie (K)]
            r {[Float]} -- [factor de relajación]
            V {[Float]} -- [Temperatura de la superficie (K)]
            T_estatica {[Float]} -- [Temperatura de la superficie (K)]

        Returns:
            [q_cond] -- [calor de conduccion (cal/(m2 s))]
        """         
        if self.calor_conductivo:
            calor = k_i*Delta/chi/l*(T_surface-T_estatica-r*V**2/(1004.5*2))
        else: calor=0
        return calor

    def coeficiente_convectivo_lineal(self,a,b,T_film):
        return a*T_film+b
    def balance_calores(self,q_c,q_e,q_w,q_rb,q_f,q_i,q_k,q_in):
        '''
        Calculo de calores del Manual de Scaling methods
        ver seccion 3.5
        los calores deben ir en las mismas unidades
        predeterminadas: cal/(m^2 s)
        esta ecuacion deberia ser cero
        
        '''
        calores_perdidos = q_c +q_e +q_w+q_rb
        calores_ganados = q_f+q_i+q_k+q_in
        return calores_ganados-calores_perdidos

    def heat_water_into(self,dot_mi,T_superficie_anterior,cp_ws_anterior):
        """[Flujo de calor que entra en el volumen de control por el V.C anterior]

        Args:
            dot_mi ([float]): [flujo másico desde la superficie anterior]
            T_superficie_anterior ([float]): [Temperatura en K de la superficie anterior]
            cp_ws_anterior ([float]): [Calor específico d ela superficie anterior]

        Returns:
            [float]: [dot_mi*cp_ws_anterior*(T_superficie_anterior-273.15) ]
        """        
        return dot_mi*cp_ws_anterior*(T_superficie_anterior-273.15)

    def set_tamano_gota(self,D_gota):
        self.D_gota =D_gota


    def calculo_todos_calores(self):
        """[Calcula todos los calores del problema. Ver la siguiente bibliografia:
            NASA, Manual of Scaling Methods
            BERNARD L. MESSINGER, Equilibrium Temperature of an Unheated Icing Surface as a Function of Air Speed]

        Arguments:
            T_pared_inicial {[float]} -- [Temperatura de la pared antes de la nebulizacion (K)]
            T_modelo {[float]} -- [Temperatura de la pared durante la nebulizacion (K)]
            T_superficie {[float]} -- [Temperatura de la superficie del hielo durante la nebulizacion (K)]
            recovery_factor {[float]} -- [factor de recuperacion ]
            p_remanso {[float]} -- [Presión de remanso del problema (Pa) suele ser 10^5]
            LWC {[float]} -- [contenido de agua líquida (g/m^3)]
            cuerda {[float]} -- [cuerda del perfil (cm)]
            D_gota {[float]} -- [diametro de la gota (m)]
            d {[float]} -- [diametro caracteristico del perfil (cm)]
            Grosor_hielo {[float]} -- [Maximo espesor que se va a computar (m)]
            V {[float]} -- [Velocidad del flujo (m/s)]
            n_0 {[float]} -- [freezing fraction (adimensional)]

        Returns:
            [tuple] -- [(q_convectivo,q_evaporacion,q_w,q_rb,q_cond,q_f,q_i,q_k) calores (cal/(m^2 s))]
        """        
        T_remanso = self.T_remanso
        T_superficie = self.T_superficie
        
        recovery_factor=self.recovery_factor
        p_remanso=self.p_remanso
        LWC=self.LWC
        d=self.d
        V=self.V
        n_0=self.n_0
        betha=self.betha
        m_in=self.m_in
        T_s_anterior=self.T_s_anterior
        cp_ws_anterior = self.cp_ws_anterior
        V_e =self.V_e
        T_estatica =self.T_estatica
        p_estatica = p_remanso/(1+V**2/(2*287*T_estatica))#Pa
        rho_a = p_estatica/(.287*T_estatica) #g/m^3 (incompresible M<0.3)
        T_film = self.Temperatura_film(T_estatica,T_superficie)
        k_a = self.thermal_conductivity(T_film) # cal/(hr m K)
        mu_a = self.viscosity(T_estatica)
        mu_film=self.viscosity(T_film)
        Re_film=self.Reynolds(V, d, rho_a, mu_film)*10**-4
        Pr = self.Pr(self.cp_a, mu_film, k_a/360000)
        self.Prandtl = Pr
        Nu = self.Nusselt(Pr, Re_film)
        #h_c = self.film_coefficient(k_a, Nu, d/100) #cal/(m K hr)
        #h_c = h_c/3600 #cal/(m K s)
        h_c = self.h_c
        difusividad = self.difusividad(T_film,p_estatica) #cm^2/s
        difusividad=difusividad*10**(-4)
        Sc = self.Schmidt(mu_a*100,rho_a,difusividad)
        h_g = self.convective_mass_coeff(h_c, Sc, Pr) #g/(m^2 s)
        Lambda_v =  self.Latent_heat_vaporisation(T_superficie) #cal/g
        p_w = self.presion_vapor(T_estatica)
        p_ww = self.presion_vapor(T_superficie)
        dot_me=self.mass_water_evaporation(p_ww,p_w,h_g,p_estatica,T_superficie,T_remanso,p_remanso,n_0) #g/(m^2 s)
        self.m_e =dot_me
        dot_m = self.impinging_mass_water(LWC, V, betha) # g m^-2 s^-1
        alpha_a = k_a/3600/rho_a/self.cp_a
        self.m_c = dot_m
        cp_ws =self.cp_ws(T_superficie) ##cal/g K
        Lambda_f = self.latent_heat_freezing(T_superficie) #cal/g
        cp_is=self.cp_is(T_superficie) #(cal/g K)
        q_1 = dot_m*(cp_ws*(T_estatica-273.15)+V**2/2/4.1868/1000)
        q_2 = m_in*cp_ws_anterior*(T_s_anterior-273.15)
        q_3=0
        if T_superficie>=273.15:
            if dot_me >0:
               
                # q_3=dot_me*(cp_ws*(-T_superficie+273.15)+Lambda_v)
                q_3=dot_me*(Lambda_v)
        else:
            q_3 = 0
            
        q_3=0
        m_out = (1-n_0)*(dot_m+m_in)-dot_me
        if m_out < 0: 
            m_out=0
        self.m_out = m_out
        
        q_4 = m_out*cp_ws*(T_superficie-273.15)
        q_5 = n_0*(dot_m+m_in)*(cp_is*(T_superficie-273.15)-Lambda_f)
        q_6 = h_c*(T_superficie-T_estatica-recovery_factor*V_e**2/(2*1004.5))
        
        return q_1+q_2-q_3-q_4-q_5-q_6


   

    def calculo_factor_coleccion(self,D_gota):
        T_remanso = self.T_remanso
        T_superficie = self.T_superficie
        p_remanso=self.p_remanso
        LWC=self.LWC
        d=self.d
        V=self.V
        T_estatica =self.T_estatica
        p_estatica = p_remanso/(1+V**2/(2*287*T_estatica))#Pa
        rho_a = p_estatica/(.287*T_estatica) #g/m^3 (incompresible M<0.3)
        T_film = self.Temperatura_film(T_estatica,T_superficie)
        k_a = self.thermal_conductivity(T_film)
        mu_a = self.viscosity(T_estatica)
        mu_film=self.viscosity(T_film)
        Re_film=self.Reynolds(V, d, rho_a, mu_film)*10**-4
        Pr = self.Pr(self.cp_a, mu_film, k_a/360000)
        self.Prandtl = Pr
        K=self.inertia_parameter(1e6,D_gota, V, d, mu_a)
        LAMBDA_LAMBDAstokes= self.dimensionless_range_parameter(V,self.D_gota,rho_a,mu_a*100)
        K_0 = self.modified_inertia_parameter(LAMBDA_LAMBDAstokes, K)
        betha_0 =self.catch_efficiency_stagnation(K_0)
        self.betha_stagnation=betha_0
        return betha_0
    
    
    
    def coeficiente_convectivo_nodos(self,x,x_nodos,HCT):
        """[Interpola el valor del coeficiente de calor convectivo utilizando los valores previamente calculados en los puntos de toma de datos]

        Args:
            x ([float]): [posción de la cuerda]
            x_nodos ([list]): [Posiciones de los nodos de estudio]
            HCT ([list]): [coeficientes de calor convectivo en cada nodo]

        Returns:
            [float]: [Coeficiente de calor convectrivo en el punto de estudio]
        """        
        hct =HCT[0]
        for i in range(1,len(x_nodos)): 
            if x_nodos[i-1]<=x<=x_nodos[i]:
                hct=HCT[i-1]+(HCT[i]-HCT[i-1])/(x_nodos[i]-x_nodos[i-1])*(x-x_nodos[i-1])
        return hct
    
    def coll_efficiency(self,x,betha):
        """[Interpola en valor de la eficiencia de colección a partir de los datos de los nodos]

        Args:
            x ([type]): [description]
            betha ([type]): [description]

        Returns:
            [type]: [description]
        """        
        x_nodos = betha[:,0]
        bethas = betha[:,1]
        resultado=0
        for i in range(1,len(x_nodos)):
           
            if (x<=x_nodos[i]) and (x>=x_nodos[i-1]):
                resultado = (x-x_nodos[i-1])/(x_nodos[i]-x_nodos[i-1])*(bethas[i]-bethas[i-1])+bethas[i-1]
                break
        return resultado

    def set_zona_perfil(self,zona):
        """[Define la zona del perfil a estudiar]

        Args:
            zona ([char]): [Puede ser 'intrados' o 'extrados']
        """        
        
        if zona != 'intrados' or zona != 'extrados':
            self.Zona = zona
        else:
            print('no valido')
            ffffff

    
    def calculate_T_sur(self):
        """[Calcula el perfil de temperaturas a raíz de unos datos previamente fijados.
        Datos a fijar:
        
        (Elementos,Velocidades,Superficie,Puntos)=Procesado_aerodinamico(70,angulo_ataque)
        Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie,Puntos)
        Modelo.set_T_remanso(273.15+T_remanso)
        Modelo.set_presion_remanso(1e5)
        betha=pickle.load(open('betha'+str(V_inf)+'_'+str(angulo_ataque)+'_'+str(MVD)+'.p', "rb"))
        Modelo_termico = fem_velocidades.analisis_termico()
        Modelo_termico.calculo_S(Modelo)
        Modelo_termico.set_zona_perfil(zona)
        Modelo_termico.set_T_remanso(T_remanso+273.15)
        Modelo_termico.x_experimental = x_experimental
        Modelo_termico.T_experimental = T_experimental
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
        Modelo_termico.betha_nodos = betha
        h_l=[Modelo_termico.coeficiente_convectivo_laminar(Modelo_termico.longitud_equivalente(i)) for i in x_experimental]
        h_l[0]=h_l[1]
        
        Modelo_termico.coeficiente_convectivo_valores = h_l

        Returns:
            [tuple]: [Tupla con los valores de x,T_sur (ºC)]
        """        
        zona = self.Zona
        if zona =='intrados':
            zona_coeff = -1
        if zona =='extrados':
            zona_coeff = 1
        T_experimental = self.T_experimental
        x_experimental = self.x_experimental
        HCT = self.coeficiente_convectivo_valores 
        Modelo = self.Modelo_CFD
        betha = self.betha_nodos
        self.punto_remanso(self.Modelo_CFD)
        
        x =np.linspace( self.x_remanso,zona_coeff*x_experimental[-1],100)
        x_experimental = x_experimental[:]*zona_coeff
        col_eff = betha[:,1]
        n = np.zeros(len(x))
        iter_max =10000
        T_sur  = np.zeros(len(x))
        Delta_T_sur  = np.zeros(len(x))
        freezing_fraction = np.empty(len(x))
        for nodo in range(len(x)):
            self.set_coeficiente_convectivo(self.coeficiente_convectivo_nodos(x[nodo],x_experimental,HCT))
            for i in range(1,len(Modelo.V_infinito[:,0])):
                if (x[nodo]<= Modelo.V_infinito[i,0]) and  (x[nodo]>= Modelo.V_infinito[i-1,0]):
                    Ve=Modelo.V_infinito[i,1]
            self.T_estatica = self.T_remanso -Ve**2/2/1004.5
            self.V_e = Ve
            self.set_T_superficie(273.15)
            self.set_local_collection_efficiency(self.coll_efficiency(x[nodo],betha))
            T_sur_0 =  self.T_remanso
            if nodo !=0:
                r = np.sqrt(self.Prandtl)
                T_sur_0 = self.T_estatica+r*Ve**2/2/1004.5
                self.set_recovery_factor(np.sqrt(self.Prandtl))
        
            n_0 =0.1
            n_1 = 0.2
            for i in range(iter_max):
                self.set_freezing_fraction(n_0)
                f_0 =self.calculo_todos_calores()
                self.set_freezing_fraction(n_1)
                f_1 =self.calculo_todos_calores()
                df = (f_1-f_0)/(n_1-n_0)
                if df < 1e-6:
                    n_1 = 1
                    break
                n_2 = n_1 -f_1/df
                if abs(f_1)<0.1:
                    break
                n_0 = n_1
                n_1 =n_2
            if i ==iter_max-1:
        
                n_1 =1
                break
            if np.abs(n_1)>1 or n_1<0:
                n_1 =1
                self.set_freezing_fraction(1)
                #print(self.n_0)
                T_0 =273.15
                T_1 =260
                for i in range(iter_max): 
                    self.set_T_superficie(T_0)
                    f_0 =self.calculo_todos_calores()
                    self.set_T_superficie(T_1)
                    f_1 =self.calculo_todos_calores()
                    df = (f_1-f_0)/(T_1-T_0)
                    T_2 = T_1 -f_1/df
                    if abs(f_1)<0.1:
                        break
                    T_0 = T_1
                    T_1 =T_2
            if n_1>=0 and n_1 <=1: 
                freezing_fraction[nodo]=n_1
            else: 
                freezing_fraction[nodo]=np.nan
            self.set_flujo_masico_entrada(self.m_out)
            
            if nodo ==0: self.set_flujo_masico_entrada(self.m_out/2)
            
            self.set_T_superficie_anterior(self.T_superficie)
            self.set_cp_ws_anterior(self.cp_ws(self.T_superficie))
            T_sur[nodo] = self.T_superficie
            Delta_T_sur[nodo]=T_sur[nodo]-T_sur_0
            
        for i,T in enumerate(T_sur):
            if T==0:
                T_sur[i]= np.nan
            else:
                T_sur[i]= T-273.15
        self.Fraccion_congelacion = freezing_fraction
        return (x,T_sur)

    def punto_remanso(self,Modelo_aerodinamico):
        """[Calcula el punto de remanso y lo guarda en la variable del objeto x_remanso]

        Args:
            Modelo_aerodinamico ([object]): [Modelo CFD aerodinámico resuelto con OpenFoam]
        """        

        
        for i in range(len(Modelo_aerodinamico.V_infinito[:,1])):
            if Modelo_aerodinamico.V_infinito[i,1]==min(Modelo_aerodinamico.V_infinito[:,1]):
                self.x_remanso = Modelo_aerodinamico.V_infinito[i,0]
                break


    def calculo_S(self,Modelo):
        """[Guarda en el objeto un array S donde la primera columna es la s la segunda el valor x y la tercera la y]

        Args:
            Modelo ([object]): [Objeto del CFD resulto con OpenFoam]
        """        
        self.punto_remanso(Modelo)
        x_0 = self.x_remanso
       
        x = list(Modelo.Superficie.x)
        y = list(Modelo.Superficie.y)
        S=[]
        for i in range(len(Modelo.Superficie.index)):
            
            if x_0<0:
                if x[i+1]<=abs(x_0) and x[i]>=abs(x_0):
                    S.append([0,x[i],y[i]])
                    n=i
                    break
            else:
                if x[i+1]>=x_0 and x[i]<=x_0:
                    S.append([0,x[i],y[i]])
                    n=i
                    break
            
        for i in range(n+1,len(Modelo.Superficie.index)):
            S.append([0,x[i],y[i]])
            s=0
            for j in range(len(S)-1):
                s=s+np.sqrt((x[j+1]-x[j])**2+(y[j+1]-y[j])**2)
            S[-1][0]=s
        self.S=np.array(S)
        

    def longitud_equivalente(self,x):
        """[Devuelve el valor de s cuando se le proporciona un valor de la x]

        Args:
            x ([float]): [valor horizontal de la cuerda (m)]

        Returns:
            [float]: [longitud de arco desde el punto de remanso]
        """        
        x_perfil = self.S[:,1]

        S = self.S[:,0]
        
        for i in range(1,len(x_perfil)):
           
            if x>=x_perfil[i-1] and x<=x_perfil[i]:
                s = S[i-1]+(S[i]-S[i-1])/(x_perfil[i]-x_perfil[i-1])*(x-x_perfil[i-1])
                break
        return s 

    def V_S(self,s):
        """[Devuelve la velocidad exterior cuando se le proporciona un valor de arco desde el punto de remanso]

        Args:
            s ([float]): [longitud de arco desde el punto de remanso]

        Returns:
            [float]: [velocidad en ese punto (m/s)]
        """        
        Modelo=self.Modelo_CFD
        V = Modelo.V_infinito[:,1]
        x = Modelo.V_infinito[:,0]
        s_0 = 0
        for i in range(1,len(x)):
            if x[i]>=0:
                s_1 = self.longitud_equivalente(x[i])
                if s>=s_0 and s<=s_1:
                    v=(V[i]+V[i-1])/2
                   
                    break
                s_0=s_1
        return v

    def capa_limite_termica(self,s):
        """[Devuelve el espesor de la capa límite térmica. ]

        Args:
            s ([float]): [longitud de arco]

        Returns:
            [float]: [espesor de la capa límite térmica]
        """        
        S = self.S[:,0]
        rho_a = self.p_remanso/(.287*self.T_remanso)
        nu = self.viscosity(self.T_remanso)/rho_a
        
        Integral=0
        for i in range(1,len(S)):
            if s<=S[i] and s>=S[i-1]:
                Integral = Integral+(self.V_S(S[i])-self.V_S(S[i-1]))**1.87*(S[i]-S[i-1])/2 + (self.V_S(S[i-1]))**1.87*(S[i]-S[i-1])
            else:
                Integral = Integral+((self.V_S(s)-self.V_S(S[i-1])))**1.87*(s-S[i-1])/2 + (self.V_S(S[i-1]))**1.87*(s-S[i-1])
                break
        delta_T=np.sqrt(nu*46.72/self.V_S(s)**2.87*Integral)
        return delta_T

    def coeficiente_convectivo_laminar(self,s):
        """[Coeficiente de calor convectivo]

        Args:
            s ([float]): [longitud del arcotion]

        Returns:
            [h_l]: [cal/(K m^2)]
        """        
        delta_T=self.capa_limite_termica(s)
        k=self.thermal_conductivity(self.T_remanso)/3600 # cal/(s m K)
        return 2*k/delta_T

    def coeficiente_convectivo_turbulento(self,x):
        """[Calcula el coeficiente de calor convectivo turbulento. ]

        Args:
            x ([type]): [description]
        """        
    
    def calculo_thickness_rate(self,x_sensor,T_sensor,T_experimental,x_experimental,HCT,Modelo_CFD,T_remanso,zona,V):
        """[Calcula la betaLWC y la velocidad de crecimiento del hielo en el perfil]

        Args:
            x_sensor ([float]): [Posicion del sensor en metros]
            T_sensor ([float]): [Temperatura del sensor en metros]
            T_experimental ([list]): [nodos de T para el HCT]
            x_experimental ([list]): [nodos de x para el HCT]
            HCT ([list]): [Coeficiente de calor convectivo]
            Modelo_CFD ([object]): [Objeto del modelo CFD]
            T_remanso ([float]): [Temperatura de remanso]
            zona ([str]): [intrados/extrados]
            V ([float]): [Velocidad aerodinámica]

        Returns:
            [tuple]: [Tupla en la que el primer término es la accretion rate y el segundo (LWC beta)]
        """        
        zona = self.Zona
        if zona =='intrados':
            zona_coeff = -1
        if zona =='extrados':
            zona_coeff = 1
        Modelo = Modelo_CFD
        self.V = V
        self.T_superficie = T_sensor
        self.punto_remanso(self.Modelo_CFD)
        x_experimental = x_experimental[:]*zona_coeff
        # print(x,col_eff)
        # print(self.betha_stagnation,col_eff[0])
        self.set_coeficiente_convectivo(self.coeficiente_convectivo_nodos(x_sensor,x_experimental,HCT))
        for i in range(1,len(Modelo.V_infinito[:,0])):
            if (x_sensor<= Modelo.V_infinito[i,0]) and  (x_sensor>= Modelo.V_infinito[i-1,0]):
                Ve=Modelo.V_infinito[i,1]
        self.T_estatica = T_remanso -Ve**2/2/1004.5
        self.V_e = Ve
        h_c = self.h_c  #[cal/(K m^2)]
        rho_ice =0.917e6 #g/m3
        cp_ws= self.cp_ws(T_sensor) #cal/(g K)
        T_st = self.T_estatica #K
        Lambda_f = self.latent_heat_freezing(T_sensor) #cal/g
        cp_is = self.cp_is(T_sensor) #(cal/g K)
        #self.calculo_todos_calores()
        print(T_sensor,T_st,self.t_f)
        print(f'h_c = {h_c} cal/(K m^2)')
        print(f'rho_ice = {rho_ice} g/m3')
        print(f'cp_ws = {cp_ws} cal/(K g)')
        print(f'T_st = {T_st} cal/(K m^2)')
        print(f'cp_is = {cp_is} cal/(K g)')
        print(f'Lambda_f = {Lambda_f} cal/g')
        dm_dt = h_c*(T_sensor-T_st-0.85*V**2/2/1004.5)/(Lambda_f+V**2/2/4.1868/1000+(self.t_f-T_sensor)*cp_is-cp_ws*(self.t_f-T_st))
        thickness_rate=dm_dt/rho_ice
        return (thickness_rate,dm_dt/V)

    def calculate_T_sur_rime(self):
        """[Calcula el perfil de temperaturas a raíz de unos datos previamente fijados.
        Datos a fijar:
        
        (Elementos,Velocidades,Superficie,Puntos)=Procesado_aerodinamico(70,angulo_ataque)
        Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie,Puntos)
        Modelo.set_T_remanso(273.15+T_remanso)
        Modelo.set_presion_remanso(1e5)
        betha=pickle.load(open('betha'+str(V_inf)+'_'+str(angulo_ataque)+'_'+str(MVD)+'.p', "rb"))
        Modelo_termico = fem_velocidades.analisis_termico()
        Modelo_termico.calculo_S(Modelo)
        Modelo_termico.set_zona_perfil(zona)
        Modelo_termico.set_T_remanso(T_remanso+273.15)
        Modelo_termico.x_experimental = x_experimental
        Modelo_termico.T_experimental = T_experimental
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
        Modelo_termico.betha_nodos = betha
        h_l=[Modelo_termico.coeficiente_convectivo_laminar(Modelo_termico.longitud_equivalente(i)) for i in x_experimental]
        h_l[0]=h_l[1]
        
        Modelo_termico.coeficiente_convectivo_valores = h_l

        Returns:
            [tuple]: [Tupla con los valores de x,T_sur (ºC)]
        """        
        zona = self.Zona
        if zona =='intrados':
            zona_coeff = -1
        if zona =='extrados':
            zona_coeff = 1
        T_experimental = self.T_experimental
        x_experimental = self.x_experimental
        HCT = self.coeficiente_convectivo_valores 
        Modelo = self.Modelo_CFD
        betha = self.betha_nodos
        self.punto_remanso(self.Modelo_CFD)
        
        x =self.x_surface
        if zona=='intrados':
            x = list(filter(lambda number: number < self.x_remanso, x))
        else:
            x = list(filter(lambda number: number >= self.x_remanso, x))
        x_experimental = x_experimental[:]*zona_coeff
        col_eff = betha[:,1]
        n = np.zeros(len(x))
        iter_max =10000
        T_sur  = np.zeros(len(x))
        Delta_T_sur  = np.zeros(len(x))
        freezing_fraction = np.empty(len(x))
        for nodo in range(len(x)):
            self.set_coeficiente_convectivo(self.coeficiente_convectivo_nodos(x[nodo]-self.x_remanso,x_experimental,HCT))
            for i in range(1,len(Modelo.V_infinito[:,0])):
                if (x[nodo]<= Modelo.V_infinito[i,0]) and  (x[nodo]>= Modelo.V_infinito[i-1,0]):
                    Ve=Modelo.V_infinito[i,1]
            self.T_estatica = self.T_remanso -Ve**2/2/1004.5
            self.V_e = Ve
            self.set_T_superficie(273.15)
            self.set_local_collection_efficiency(self.coll_efficiency(x[nodo],betha))
            T_sur_0 =  self.T_remanso
            if nodo !=0:
                r = np.sqrt(self.Prandtl)
                T_sur_0 = self.T_estatica+r*Ve**2/2/1004.5
                self.set_recovery_factor(np.sqrt(self.Prandtl))
        
           
            n_1 =1
            self.set_freezing_fraction(1)
            #print(self.n_0)
            T_0 =273.15
            T_1 =260
            for i in range(iter_max): 
                self.set_T_superficie(T_0)
                f_0 =self.calculo_todos_calores()
                self.set_T_superficie(T_1)
                f_1 =self.calculo_todos_calores()
                df = (f_1-f_0)/(T_1-T_0)
                T_2 = T_1 -f_1/df
                if abs(f_1)<0.1:
                    break
                T_0 = T_1
                T_1 =T_2
            if n_1>=0 and n_1 <=1: 
                freezing_fraction[nodo]=n_1
            else: 
                freezing_fraction[nodo]=np.nan
            self.set_flujo_masico_entrada(self.m_out)
            
            if nodo ==0: self.set_flujo_masico_entrada(self.m_out/2)
            
            self.set_T_superficie_anterior(self.T_superficie)
            self.set_cp_ws_anterior(self.cp_ws(self.T_superficie))
            T_sur[nodo] = self.T_superficie
            Delta_T_sur[nodo]=T_sur[nodo]-T_sur_0
            
        for i,T in enumerate(T_sur):
            if T==0:
                T_sur[i]= np.nan
            else:
                T_sur[i]= T-273.15
        self.Fraccion_congelacion = freezing_fraction
        return (x,T_sur)

    def calculate_T_sur_rime_punto(self,x): 
        
        T_experimental = self.T_experimental
        x_experimental = self.x_experimental
        HCT = self.coeficiente_convectivo_valores 
        Modelo = self.Modelo_CFD
        betha = self.betha_nodos
        self.punto_remanso(self.Modelo_CFD)
        if x < self.x_remanso:
            zona = 'intrados'
        else: zona='extrados'
        if zona =='intrados':
            zona_coeff = -1
        if zona =='extrados':
            zona_coeff = 1
        x_experimental = x_experimental[:]*zona_coeff
        col_eff = betha[:,1]
        iter_max =10000
      
        self.set_coeficiente_convectivo(self.coeficiente_convectivo_nodos(x-self.x_remanso,x_experimental,HCT))
        for i in range(1,len(Modelo.V_infinito[:,0])):
            if (x<= Modelo.V_infinito[i,0]) and  (x>= Modelo.V_infinito[i-1,0]):
                Ve=Modelo.V_infinito[i,1]
        self.T_estatica = self.T_remanso -Ve**2/2/1004.5
        self.V_e = Ve
        self.set_T_superficie(273.15)
        self.set_local_collection_efficiency(self.coll_efficiency(x,betha))
        T_sur_0 =  self.T_remanso
        r = np.sqrt(self.Prandtl)
        T_sur_0 = self.T_estatica+r*Ve**2/2/1004.5
        self.set_recovery_factor(np.sqrt(self.Prandtl))
        n_1 =1
        self.set_freezing_fraction(1)
        T_remanso = self.T_remanso
        T_superficie = self.T_superficie
        
        recovery_factor=self.recovery_factor
        p_remanso=self.p_remanso
        LWC=self.LWC
        d=self.d
        V=self.V
        n_0=self.n_0
        betha=self.betha
        m_in=self.m_in
        T_s_anterior=self.T_s_anterior
        cp_ws_anterior = self.cp_ws_anterior
        V_e =self.V_e
        T_estatica =self.T_estatica
        p_estatica = p_remanso/(1+V**2/(2*287*T_estatica))#Pa
        rho_a = p_estatica/(.287*T_estatica) #g/m^3 (incompresible M<0.3)
        T_film = self.Temperatura_film(T_estatica,T_superficie)
        k_a = self.thermal_conductivity(T_film) # cal/(hr m K)
        mu_a = self.viscosity(T_estatica)
        mu_film=self.viscosity(T_film)
        Re_film=self.Reynolds(V, d, rho_a, mu_film)*10**-4
        Pr = self.Pr(self.cp_a, mu_film, k_a/360000)
        self.Prandtl = Pr
        Nu = self.Nusselt(Pr, Re_film)
        #h_c = self.film_coefficient(k_a, Nu, d/100) #cal/(m K hr)
        #h_c = h_c/3600 #cal/(m K s)
        h_c = self.h_c
        difusividad = self.difusividad(T_film,p_estatica) #cm^2/s
        difusividad=difusividad*10**(-4)
        Sc = self.Schmidt(mu_a*100,rho_a,difusividad)
        h_g = self.convective_mass_coeff(h_c, Sc, Pr) #g/(m^2 s)
        Lambda_v =  self.Latent_heat_vaporisation(T_superficie) #cal/g
        p_w = self.presion_vapor(T_estatica)
        p_ww = self.presion_vapor(T_superficie)
        dot_me=self.mass_water_evaporation(p_ww,p_w,h_g,p_estatica,T_superficie,T_remanso,p_remanso,n_0) #g/(m^2 s)
        self.m_e =dot_me
        dot_m = self.impinging_mass_water(LWC, V, betha) # g m^-2 s^-1
        alpha_a = k_a/3600/rho_a/self.cp_a
        self.m_c = dot_m
        cp_ws =self.cp_ws(T_superficie) ##cal/g K
        Lambda_f = self.latent_heat_freezing(T_superficie) #cal/g
        cp_is=self.cp_is(T_superficie) #(cal/g K)
        T_sur=dot_m*cp_ws*(T_estatica-273.15)
        T_sur+=dot_m*V**2/2/4.1868/1000
        T_sur+=dot_m*cp_is*273.15
        T_sur+=dot_m*Lambda_f
        T_sur+=h_c*(T_estatica+r*V**2/2/1004.5)
        T_sur=T_sur/(h_c+dot_m*cp_is)
        
        return (T_sur)

    def calculate_T_sur_borde_ataque(self):
        """[Calcula el perfil de temperaturas a raíz de unos datos previamente fijados.
        Datos a fijar:
        
        (Elementos,Velocidades,Superficie,Puntos)=Procesado_aerodinamico(70,angulo_ataque)
        Modelo = fem_velocidades.modelo_fem(Elementos,Velocidades,Superficie,Puntos)
        Modelo.set_T_remanso(273.15+T_remanso)
        Modelo.set_presion_remanso(1e5)
        betha=pickle.load(open('betha'+str(V_inf)+'_'+str(angulo_ataque)+'_'+str(MVD)+'.p', "rb"))
        Modelo_termico = fem_velocidades.analisis_termico()
        Modelo_termico.calculo_S(Modelo)
        Modelo_termico.set_zona_perfil(zona)
        Modelo_termico.set_T_remanso(T_remanso+273.15)
        Modelo_termico.x_experimental = x_experimental
        Modelo_termico.T_experimental = T_experimental
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
        Modelo_termico.betha_nodos = betha
        h_l=[Modelo_termico.coeficiente_convectivo_laminar(Modelo_termico.longitud_equivalente(i)) for i in x_experimental]
        h_l[0]=h_l[1]
        
        Modelo_termico.coeficiente_convectivo_valores = h_l

        Returns:
            [tuple]: [Tupla con los valores de x,T_sur (ºC)]
        """        
        zona = self.Zona
        if zona =='intrados':
            zona_coeff = -1
        if zona =='extrados':
            zona_coeff = 1
        T_experimental = self.T_experimental
        x_experimental = self.x_experimental
        HCT = self.coeficiente_convectivo_valores 
        Modelo = self.Modelo_CFD
        betha = self.betha_nodos
        self.punto_remanso(self.Modelo_CFD)
        
        x =np.linspace( self.x_remanso,zona_coeff*x_experimental[-1],100)
        x_experimental = x_experimental[:]*zona_coeff
        col_eff = betha[:,1]
        n = np.zeros(len(x))
        iter_max =10000
        T_sur  = np.zeros(len(x))
        Delta_T_sur  = np.zeros(len(x))
        freezing_fraction = np.empty(len(x))
        for nodo in range(len(x)):
            self.set_coeficiente_convectivo(self.coeficiente_convectivo_nodos(x[nodo],x_experimental,HCT))
            for i in range(1,len(Modelo.V_infinito[:,0])):
                if (x[nodo]<= Modelo.V_infinito[i,0]) and  (x[nodo]>= Modelo.V_infinito[i-1,0]):
                    Ve=Modelo.V_infinito[i,1]
            self.T_estatica = self.T_remanso -Ve**2/2/1004.5
            self.V_e = Ve
            self.set_T_superficie(273.15)
            self.set_local_collection_efficiency(self.coll_efficiency(x[nodo],betha))
            T_sur_0 =  self.T_remanso
            if nodo !=0:
                r = np.sqrt(self.Prandtl)
                T_sur_0 = self.T_estatica+r*Ve**2/2/1004.5
                self.set_recovery_factor(np.sqrt(self.Prandtl))
        
            n_0 =0.1
            n_1 = 0.2
            for i in range(iter_max):
                self.set_freezing_fraction(n_0)
                f_0 =self.calculo_todos_calores()
                self.set_freezing_fraction(n_1)
                f_1 =self.calculo_todos_calores()
                df = (f_1-f_0)/(n_1-n_0)
                if df < 1e-6:
                    n_1 = 1
                    break
                n_2 = n_1 -f_1/df
                if abs(f_1)<0.1:
                    break
                n_0 = n_1
                n_1 =n_2
            if i ==iter_max-1:
        
                n_1 =1
                break
            if np.abs(n_1)>1 or n_1<0:
                n_1 =1
                self.set_freezing_fraction(1)
                #print(self.n_0)
                T_0 =273.15
                T_1 =260
                for i in range(iter_max): 
                    self.set_T_superficie(T_0)
                    f_0 =self.calculo_todos_calores()
                    self.set_T_superficie(T_1)
                    f_1 =self.calculo_todos_calores()
                    df = (f_1-f_0)/(T_1-T_0)
                    T_2 = T_1 -f_1/df
                    if abs(f_1)<0.1:
                        break
                    T_0 = T_1
                    T_1 =T_2
            if n_1>=0 and n_1 <=1: 
                freezing_fraction[nodo]=n_1
            else: 
                freezing_fraction[nodo]=np.nan
            self.set_flujo_masico_entrada(self.m_out)
            
            if nodo ==0: self.set_flujo_masico_entrada(self.m_out/2)
            
            self.set_T_superficie_anterior(self.T_superficie)
            self.set_cp_ws_anterior(self.cp_ws(self.T_superficie))
            T_sur[nodo] = self.T_superficie
            Delta_T_sur[nodo]=T_sur[nodo]-T_sur_0
            break
            
        for i,T in enumerate(T_sur):
            if T==0:
                T_sur[i]= np.nan
            else:
                T_sur[i]= T-273.15
        self.Fraccion_congelacion = freezing_fraction
        return (x[0],T_sur[0])
    



        


        

        

    


 