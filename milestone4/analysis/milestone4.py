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

    


    cl         = read_file(output_path+"cls.dat")
    ls         = read_file(output_path+"ls.dat")

    sjl        = read_file(output_path+"sjl.dat")
    x        = read_file(output_path+"x.dat")


    plt.plot(x,sjl/1e-3)
    plt.show()

    cl_norm = cl/(2.*np.pi)
    cl_norm = cl_norm/np.max(cl_norm)*5575

    plt.plot(ls,cl_norm)
    plt.ylim(0,6000)
    plt.show()
