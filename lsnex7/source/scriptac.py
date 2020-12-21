import numpy as np
import statsmodels.tsa.stattools as st
import sys

ndatas=int(sys.argv[1])
par=['U','P']
phase=['solid','liquid','gas']

for p in par:
    for f in phase:
        try: data=np.loadtxt('../data/measures/'+p+'_inst.'+f, unpack=True, max_rows=ndatas)
        except: 
            print('file ../data/measures/'+p+'_inst.'+f+' not found')
            break
        np.savetxt('../data/measures/'+p+'_ac.'+f, st.acf(data, fft=False, nlags=ndatas-1), newline='\n')
