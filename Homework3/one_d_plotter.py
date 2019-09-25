import numpy as np
import matplotlib.pyplot as plt

Dt = {"0.100000", "0.200000","0.300000","0.400000","0.500000"}
Diff = {"centered", "downstream", "upstream"}

#use_last = False
#if (use_last):
    #plt.figure(figsize=(10,10))
    #plt.plot(data[-1,1:-2])
    #plt.savefig("figures/ForwardEuler-SineWave_LastPoint.png", dpi=200)
    #plt.close()
#
#else:
    #for i in range(len(data[:])):
        #print(i)
        #plt.figure(figsize=(10,10))
        #plt.plot(data[i,1:-2])
        #plt.plot(sine_data[i,1:-2])
        #plt.savefig("figures/ForwardEuler-centered-SineWave-dt-0.1_"+str(i)+".png", dpi=300)
        #plt.close()

#plt.show()

for dt in Dt:
    print(dt)
    sine_data = np.genfromtxt("data_files/Sinewave-dt-"+dt+".dat")
    for diff in Diff:
        print("    "+diff)
        data = np.genfromtxt("data_files/ForwardEuler-"+diff+"-SineWave-dt-"+dt+".dat")
        for i in range(len(data))[0::50]:
            print("        "+str(i))
            plt.figure(figsize=(3,3))
            plt.title("Forward Euler,"+diff+","+dt)
            plt.plot(data[i,1:-1])
            plt.plot(sine_data[i,:])
            plt.savefig("figures/ForwardEuler/"+diff+"/"+dt+"/SineWave_"+str(i)+".png", dpi=200)
            plt.close()

