import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("darkgrid")
sns.set_context("talk")

def read_file(name):
    x = []
    y = []

    with open(name) as f:
        for line in f.readlines():
            words = line.split()
            x.append(words[0])
            y.append(words[1])
            

    return np.array(x,dtype=np.float), np.array(y,dtype=np.float)




if __name__ == "__main__":
    output_path = "../output/"
    report_path = "../report/"

    x, X_e = read_file(output_path + "X_e.dat")

    z = np.exp(-x)-1

    plt.semilogy(z,X_e)
    plt.xlabel(r"$z$")
    plt.ylabel(r"$X_e$")
    plt.title("Electron Fraction")
    plt.xlim(1800,100)
    plt.savefig(report_path+"xe.png")
    plt.show()



    x, tau = read_file(output_path + "tau.dat")
    x, dtau = read_file(output_path + "dtau.dat")
    x, ddtau = read_file(output_path + "ddtau.dat")
    
    plt.semilogy(x,tau,label=r"$\tau$")
    plt.semilogy(x,abs(dtau),"--",label=r"$|\tau'|$")
    plt.semilogy(x,abs(ddtau),"--",label=r"$|\tau''|$")
    plt.xlabel(r"$x$")
    plt.ylabel(r"$\tau, |\tau'|, |\tau''|$")
    plt.title("Optical Depth")
    plt.legend()
    plt.xlim(-18,-0.5)
    #plt.ylim(1e-12,1e20)
    plt.savefig(report_path+"tau.png")
    plt.show()


    x, g = read_file(output_path + "g.dat")
    x, dg = read_file(output_path + "dg.dat")
    x, ddg = read_file(output_path + "ddg.dat")

    plt.plot(x,g,label=r"$\tilde{g}$")
    plt.plot(x,dg/10.,"--",label=r"$\tilde{g'}/10$")
    plt.plot(x,ddg/300.,".",label=r"$\tilde{g''}/300$")
    plt.legend()
    plt.xlabel(r"$x$")
    plt.ylabel(r"$\tilde{g}, \tilde{g'}/10, \tilde{g''}/300$")
    plt.title("Visibility Function")
    plt.xlim(-7.4,-6)
    plt.savefig(report_path+"g.png")
    plt.show()
