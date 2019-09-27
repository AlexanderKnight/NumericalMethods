import numpy as np
import matplotlib.pyplot as plt

Cf = ["0.100000", "0.200000","0.300000","0.400000","0.500000"]
Diff = ["centered", "downstream", "upstream"]
Method = ["ForwardEuler", "RungeKutta3"]

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
        #plt.savefig("figures/ForwardEuler-centered-SineWave-cf-0.1_"+str(i)+".png", dpi=300)
        #plt.close()

#plt.show()

for cf in Cf:
    print(cf)
    sine_data = np.genfromtxt("data_files/Sinewave-cf-"+cf+".dat")
    for method in Method:
        print("    "+method)
        for diff in Diff:
            print("        "+diff)
            data = np.genfromtxt("data_files/"+method+"-"+diff+"-SineWave-cf-"+cf+".dat")
            for i in range(len(data))[0::50]:
                print("            "+str(i))
                plt.figure(figsize=(3,3))
                plt.title(method+","+diff+","+cf)
                plt.plot(data[i,1:-1])
                plt.plot(sine_data[i,:])
                plt.savefig("figures/"+method+"/"+diff+"/"+cf+"/SineWave_"+str(i)+".png", dpi=200)
                plt.close()

