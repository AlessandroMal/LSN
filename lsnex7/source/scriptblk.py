import numpy as np

par=['U','P']
phase=['solid','liquid','gas']

for p in par:
    for f in phase:

        try: data=np.loadtxt('../data/measures/'+p+'_inst.'+f, unpack=True)
        except: 
            print('file ../data/measures/'+p+'_inst.'+f+' not found')
            break

        out=open('../data/measures/'+p+'_blockerr.'+f, 'w')
        for L in range(10,5000):
            blkmeasure=[]
            for iblk in range(int(len(data)/L)):
                blkmeasure.append(np.mean(data[iblk*L:(iblk+1)*L]))
            out.write(str(L) + ' ' + str( np.sqrt( np.mean(np.array(blkmeasure)**2) - np.mean(np.array(blkmeasure))**2 )/np.sqrt(int(len(data)/L)-1) ) +'\n' )

        out.close()
