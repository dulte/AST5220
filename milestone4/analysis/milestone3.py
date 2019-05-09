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

def plot_quantity(y,x,k,nb_modes,title,xlab,ylab,name="",save=False,log_axis=False):
    for i in range(nb_modes):
        if log_axis:
            plt.semilogy(x,y[i],label=r"$kc /H_0$ = %.1f"%k[i])    
        else:
            plt.plot(x,y[i],label=r"$kc /H_0$ = %.1f"%k[i])

    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.legend(loc="best")

    if save:
        if name == "":
            print("Please give a name for saving the image!")
            exit(1)
        plt.savefig(name+".png")
    else:
        plt.show()



if __name__ == "__main__":
    output_path = "../output/"
    report_path = "../report/"

    


    Psi         = read_file(output_path+"Psi.dat")
    Phi         = read_file(output_path+"Phi.dat")

    v           = read_file(output_path+"v.dat")
    v_b         = read_file(output_path+"v_b.dat")

    delta       = read_file(output_path+"delta.dat")
    delta_b     = read_file(output_path+"delta_b.dat")

    theta       = read_file(output_path+"Theta.dat")

    

    
    k           = read_file(output_path+"k.dat")
    x           = read_file(output_path+"x_modes.dat")
    nb_modes    = int((len(v)/len(x)))

    Psi         = Psi.reshape(nb_modes,int(len(x)))
    Phi         = Phi.reshape(nb_modes,int(len(x)))

    v           = v.reshape(nb_modes,int(len(x)))
    v_b         = v_b.reshape(nb_modes,int(len(x)))

    delta       = delta.reshape(nb_modes,int(len(x)))
    delta_b     = delta_b.reshape(nb_modes,int(len(x)))
    
    
    
    theta = theta.reshape(int((len(theta)/(len(x)*7.))),7,int(len(x)))


    plot_quantity(Phi,x,k,nb_modes,r"$\Phi$ for different modes of k","x",r"$\Phi$")
    plot_quantity(Psi,x,k,nb_modes,r"$\Psi$ for different modes of k","x",r"$\Psi$")

    plot_quantity(delta,x,k,nb_modes,r"$\delta$ for different modes of k","x",r"$\delta$",log_axis=True)
    plot_quantity(delta_b,x,k,nb_modes,r"$\delta_b$ for different modes of k","x",r"$\delta_b$",log_axis=True)

    plot_quantity(v,x,k,nb_modes,r"$v$ for different modes of k","x",r"$v$")
    plot_quantity(v_b,x,k,nb_modes,r"$v_b$ for different modes of k","x",r"$v_b$")

    plot_quantity(theta[:,0],x,k,nb_modes,r"$\Theta_0$ for different modes of k","x",r"$\Theta_0$")