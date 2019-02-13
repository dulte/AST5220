import matplotlib.pyplot as plt
import numpy as np
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
    rapport_path = "../rapport/"

    """
    Reads and plots the density fractions
    """
    x, omega_b = read_file(output_path + "Omega_b.dat")
    x, omega_r = read_file(output_path + "Omega_r.dat")
    x, omega_m = read_file(output_path + "Omega_m.dat")
    x, omega_L = read_file(output_path + "Omega_L.dat")
    

    plt.plot(x,omega_b,label=r"$\Omega_{b}$")
    plt.plot(x,omega_r,label=r"$\Omega_{r}$")
    plt.plot(x,omega_m,label=r"$\Omega_{m}$")
    plt.plot(x,omega_L,label=r"$\Omega_{\Lambda}$")
    plt.xlim(x[0],x[-1])
    plt.legend(loc="best")
    plt.xlabel("Logarithmic Scale Factor")
    plt.ylabel("Density Parameter")
    plt.title("The Evolution of the Density Parameterss")
    plt.savefig(rapport_path+"Omega.png")
    plt.show()


    """
    Reads and plots the Hubble parameter
    """

    x, H_x = read_file(output_path + "H.dat")
    z, H_z = read_file(output_path + "H_z.dat")


    plt.semilogy(x,H_x)
    plt.xlim(x[0],x[-1])
    plt.xlabel("Logarithmic Scale Factor")
    plt.ylabel("H")
    plt.title("The Evolution of the Hubble Parameters H(x)")
    plt.savefig(rapport_path+"H.png")
    plt.show()

    plt.plot(z,H_z)
    
    ax = plt.gca()
    ax.invert_xaxis()

    plt.xlabel("Red Shift")
    plt.ylabel("H")
    plt.title("The Evolution of the Hubble Parameters H(z)")
    plt.savefig(rapport_path+"H_z.png")
    plt.show()


    """
    Reads and plots the conformal time
    """

    x, eta = read_file(output_path + "eta.dat")
  
    plt.semilogy(x,eta)
    plt.xlim(x[0],x[-1])
    plt.xlabel("Logarithmic Scale Factor")
    plt.ylabel(r"$\eta$")
    plt.title(r"The Conformal Time $\eta$")
    plt.savefig(rapport_path+"eta.png")
    plt.show()
    
