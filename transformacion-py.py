import pandas as pd
import os
condicion={}
condicion.update({'4A':[0.39,20]})
condicion.update({'4B':[0.65,20]})
condicion.update({'4C':[0.95,20]})
condicion.update({'5A':[0.35,40]})
condicion.update({'5B':[0.62,40]})
condicion.update({'5C':[0.92,40]})
condicion.update({'6A':[0.33,70]})
condicion.update({'6B':[0.64,70]})
condicion.update({'6C':[0.93,70]})

for condition in condicion:
    file= open(f'Ensayos_p7_RIME_0_grados//ensayo P7 IWT {condition}.txt','r')
    primera_linea=int(file.readlines()[0])
    file.close()
    print(primera_linea)
    df = pd.read_csv(f'Ensayos_p7_RIME_0_grados//ensayo P7 IWT {condition}.txt',header=69,decimal=',',sep='\t')
    print(condition)
    print(df.columns)
    t= [.1*i for i in range(len(df.index))]
    file = open(f'Ensayos_p7_RIME_0_grados//ensayo_{condition}.txt','w')
    file.write(f'Tiempo inicio del ensayo: 2020-07-30 11:36:20.465682\nLWC (g/m3)={condicion[condition][0]}\nMVD (um)={condicion[condition][1]}\n')
    for i in range(1,17):
        file.write(f'FBG{i}\t')
    file.write(f'Tiempo (s)\n')
    for i in df.index:
        for c in range(1,9):
            file.write(str(round(df['A'+str(c)],2).loc[i])+'\t')
        for c in range(1,9):
            file.write(str(round(df['B'+str(c)],2).loc[i])+'\t')
        file.write(str(round(t[i],1))+'\n')    
    file.close()
        

    
    
