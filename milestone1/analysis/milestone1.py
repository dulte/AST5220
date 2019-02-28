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
    plt.xlabel(r"Logarithmic Scale Factor $x$")
    plt.ylabel("Density Parameter")
    plt.title("The Evolution of the Density Parameterss")
    plt.savefig(rapport_path+"Omega.png")
    plt.show()


    """
    Reads and plots the Hubble parameter
    """

    x, H_x = read_file(output_path + "H.dat")
    z, H_z = read_file(output_path + "H_z.dat")


    plt.plot(x,np.log(H_x))
    plt.xlim(x[0],x[-1])
    plt.xlabel(r"Logarithmic Scale Factor $x$")
    plt.ylabel(r"$\log H$ ")
    plt.title("The Evolution of the Hubble Parameters H(x)")
    plt.savefig(rapport_path+"H.png")
    plt.show()

    plt.plot(z,H_z)
    
    ax = plt.gca()
    ax.invert_xaxis()

    plt.xlabel("Red Shift")
    plt.ylabel(r"H $[s^{-1}]$")
    plt.title("The Evolution of the Hubble Parameters H(z)")
    plt.savefig(rapport_path+"H_z.png")
    plt.show()


    """
    Reads and plots the conformal time
    """

    x, eta = read_file(output_path + "eta.dat")
  
    plt.plot(x,np.log(eta))
    plt.xlim(x[0],x[-1])
    plt.xlabel(r"Logarithmic Scale Factor $x$")
    plt.ylabel(r"$\log \eta$")
    plt.title(r"The Conformal Time $\eta$")
    plt.savefig(rapport_path+"eta.png")
    plt.show()


    """
    Finds the slope of the logged eta to see if it holds
    up to the theory
    """

    rad_dominated = eta[x<-18]#[np.logical_and(x<-3,x>-5)]
    rad_x = x[x<-18]#[np.logical_and(x<-3,x>-5)]
    log_rad = np.log(rad_dominated)

    print("Slope for radiation dominated logged Universe:")
    print((log_rad[-1]-log_rad[3])/(rad_x[-1]-rad_x[3]))





    
