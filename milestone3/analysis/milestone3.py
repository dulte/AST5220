import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("darkgrid")
sns.set_context("talk")

def read_file(name):
    y = []

    with open(name) as f:
        for line in f.readlines():
            words = line.split()
            for w in words:
                y.append(float(w))
            

    return np.array(y,dtype=np.float)



if __name__ == "__main__":
    output_path = "../output/"
    report_path = "../report/"

    v = read_file(output_path+"Phi.dat")
    x = read_file(output_path+"x_modes.dat")

    modes = int((len(v)/len(x)))

    v = v.reshape(int((len(v)/len(x))),int(len(x)))

    print(v[0])
    for i in range(modes):
        plt.plot(x,v[i])
    plt.show()


