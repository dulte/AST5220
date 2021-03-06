import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import savgol_filter


sns.set_style("darkgrid")
sns.set_context("talk")

plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)

font_size = 35
legend_font_size = 25



def read_file(name):
    y = []

    with open(name) as f:
        for line in f.readlines():
            words = line.split()
            for w in words:
                y.append(float(w))
            

    return np.array(y,dtype=np.float)


def get_Plank_data():
    """
    The file I downloaded lacked a best fit curve, so I smoothed out the datapoints to
    make a curve. The values of the errorbars may thus be a bit wrongs, since they were meant for
    the raw datapoints.
    """
    filename = "../output/COM_PowerSpect_CMB-TT-full_R3.01.txt"

    Dl = []
    Dl_low = []
    Dl_high = []

    with open(filename) as f:
        
        for line in f.readlines():
            words = line.split()
            
            Dl.append(float(words[1]))
            Dl_low.append(float(words[2]))
            Dl_high.append(float(words[3]))


    Cl = np.array(Dl)[:1199]
    Cl = savgol_filter(Cl,19,1)
    norm_factor = 5775/np.max(Cl)
    Dl = np.array(Dl)[:1199]*norm_factor
    Dl_low = np.array(Dl_low)[:1199]*norm_factor
    Dl_high = np.array(Dl_high)[:1199]*norm_factor

    Dl_limits = np.zeros((2,len(Dl)))
    Dl_limits[0,:] = Dl_low
    Dl_limits[1,:] = Dl_high

    

    return Dl,Dl_limits

def normalize_cl(cl,cl_default):
    cl = cl[1:]

    cl_norm = cl/(2.*np.pi)
    cl_norm = cl_norm/np.max(cl_norm[100:])*5775
    
    return cl_norm[:]



def plot_default_againts_other(ls,cl_default,cl_up,cl_down,observed,observed_lim,up_label,down_label,title):
    
    
    plt.plot(ls,cl_default,"--",label="Default")
    plt.plot(ls,cl_up,label=up_label)
    plt.plot(ls,cl_down,label=down_label)
    plt.title(title,fontsize=font_size)
    plt.xlabel(r"$l$",fontsize=font_size)
    plt.ylabel(r"$C_l\cdot l(l+1)/2\pi [\mu K^2]$",fontsize=font_size)
    plt.legend(loc="best",fontsize=legend_font_size)
    plt.show()



if __name__ == "__main__":
    output_path = "../output/"
    report_path = "../report/"

    Dl,Dl_limits = get_Plank_data()

   

    

    


    cl_default = read_file(output_path+"cls_default.dat")
    cl_default_normed = normalize_cl(cl_default,cl_default)
    ls = read_file(output_path+"ls.dat")
    ls = ls[1:]

    # For default
    plt.errorbar(ls,Dl,yerr=Dl_limits,label="Observed")
    
    plt.plot(ls,cl_default_normed,label="default")
    plt.title(r"$C_l$ for simulation and observations",fontsize=font_size)
    plt.xlabel(r"$l$",fontsize=font_size)
    plt.ylabel(r"$C_l\cdot l(l+1)/2\pi [\mu K^2]$",fontsize=font_size)
    plt.legend(loc="best",fontsize=font_size)
    
    plt.show()


    """
    For plotting with different values of the cosmological parameters.
    Since the files used for this is not delivered, I've commented out the code.
    """

    """
    #For changes in n
    cl_n_up = read_file(output_path+"cls_n_up.dat")
    cl_n_down = read_file(output_path+"cls_n_down.dat")

    cl_n_up = normalize_cl(cl_n_up,cl_default)
    cl_n_down = normalize_cl(cl_n_down,cl_default)

    plot_default_againts_other(ls,cl_default_normed,cl_n_up,cl_n_down,Dl,Dl_limits,up_label=r"$n_s = 1.034$",down_label=r"$n_s = 0.9$",title=r"$C_l$ for Different $n_s$")

    #For changes in Omega_m
    cl_m_up = read_file(output_path+"cls_m_up.dat")
    cl_m_down = read_file(output_path+"cls_m_down.dat")

    cl_m_up = normalize_cl(cl_m_up,cl_default)
    cl_m_down = normalize_cl(cl_m_down,cl_default)

    plot_default_againts_other(ls,cl_default_normed,cl_m_up,cl_m_down,Dl,Dl_limits,up_label=r"$\Omega_m = 0.3$",down_label=r"$\Omega_m = 0.2$",title=r"$C_l$ for Different $\Omega_m$")


    #For changes in Omega_b
    cl_b_up = read_file(output_path+"cls_b_up.dat")
    cl_b_down = read_file(output_path+"cls_b_down.dat")

    cl_b_up = normalize_cl(cl_b_up,cl_default)
    cl_b_down = normalize_cl(cl_b_down,cl_default)

    plot_default_againts_other(ls,cl_default_normed,cl_b_up,cl_b_down,Dl,Dl_limits,up_label=r"$\Omega_b = 0.06$",down_label=r"$\Omega_b = 0.032$",title=r"$C_l$ for Different $\Omega_b$")


    #For changes in Omega_r
    cl_r_up = read_file(output_path+"cls_r_up.dat")
    cl_r_down = read_file(output_path+"cls_r_down.dat")

    cl_r_up = normalize_cl(cl_r_up,cl_default)
    cl_r_down = normalize_cl(cl_r_down,cl_default)

    plot_default_againts_other(ls,cl_default_normed,cl_r_up,cl_r_down,Dl,Dl_limits,up_label=r"$\Omega_r = 1.0\cdot 10^{-4}$",down_label=r"$\Omega_r = 6.6\cdot 10^{-5}$",title=r"$C_l$ for Different $\Omega_r$")

    #For changes in h
    cl_h_up = read_file(output_path+"cls_h_up.dat")
    cl_h_down = read_file(output_path+"cls_h_down.dat")

    cl_h_up = normalize_cl(cl_h_up,cl_default)
    cl_h_down = normalize_cl(cl_h_down,cl_default)

    plot_default_againts_other(ls,cl_default_normed,cl_h_up,cl_h_down,Dl,Dl_limits,up_label=r"$h = 0.8$",down_label=r"$h=0.6$",title=r"$C_l$ for Different $h$")


    #For best fit
    cl_best_1 = read_file(output_path+"cls_bf_1.dat")
    cl_best_1 = normalize_cl(cl_best_1,cl_default)

    cl_best_2 = read_file(output_path+"cls_bf_2.dat")
    cl_best_2 = normalize_cl(cl_best_2,cl_default)

    cl_best_6 = read_file(output_path+"cls_bf_6.dat")
    cl_best_6 = normalize_cl(cl_best_6,cl_default)

    cl_best_7 = read_file(output_path+"cls_bf_7.dat")
    cl_best_7 = normalize_cl(cl_best_7,cl_default)

    plt.errorbar(ls,Dl,yerr=Dl_limits,label="Observed")
    
    plt.plot(ls,cl_best_1,label=r"$\Omega_b = 0.06, \Omega_m = 0.2, n_s = 0.9$")
    plt.plot(ls,cl_best_2,label=r"$\Omega_b = 0.055, \Omega_m = 0.18, n_s = 0.9$")
    plt.plot(ls,cl_best_6,label=r"$\Omega_b = 0.06, \Omega_m = 0.18, n_s = 0.85$")
    plt.plot(ls,cl_best_7,label=r"$\Omega_b = 0.065, \Omega_m = 0.20, n_s = 0.80$")
    plt.title(r"$C_l$ for Best Fit",fontsize=font_size)
    plt.xlabel(r"$l$",fontsize=font_size)
    plt.ylabel(r"$C_l\cdot l(l+1)/2\pi [\mu K^2]$",fontsize=font_size)
    plt.legend(loc="best",fontsize=legend_font_size)
    plt.show()


    """


    
    
    # For plotting Source function times Bessel
    sjl = read_file(output_path+"sjl.dat")
    x = read_file(output_path+"x.dat")
    
    

    plt.plot(x,sjl/1e-3)
    plt.xlim(-8,0)
    plt.ylim(-3,2)
    plt.xlabel("x",fontsize=font_size)
    plt.ylabel(r"$\tilde{S}j_l(k(\eta_0-\eta))/10^{-3}$",fontsize=font_size)
    plt.show()


    #For plotting transfer function
    ks = np.linspace(0.1,1001,5000)
    c = 2.99792458e8
    h0 = 0.7
    Mpc = 3.08568025e22
    H_0 = h0 * 100.0 * 1e3 / Mpc
    
    



    for l in [6,200,400,800,1000,1200]:
        filename = "transfer_l_%s.dat" %l
        Theta = read_file(output_path+filename)
        plt.plot(ks[:2500],Theta[:2500],label="l=%s"%l)
    plt.legend(loc="best",fontsize=legend_font_size)
    plt.xlabel(r"$kc/H_0$",fontsize=font_size)
    plt.ylabel(r"$\Theta_l$",fontsize=font_size)
    plt.xlim(0,500)
    #plt.ylim(-.025,.025)
    plt.ylim(-.01,.03)
    plt.show()


    for l in [6,200,400,800,1000,1200]:
        filename = "transfer_l_%s.dat" %l
        Theta = read_file(output_path+filename)
        plt.plot(ks[:2000],l*(l+1)*Theta[:2000]**2/(ks[:2000]),label="l=%s"%l)
    plt.legend(loc="best",fontsize=legend_font_size)
    plt.xlabel(r"$kc/H_0$",fontsize=font_size)
    plt.ylabel(r"$l(l+1)\Theta_l^2/k \cdot H_0/c$",fontsize=font_size)
    plt.xlim(0,400)
    plt.ylim(0,0.1)
    plt.show()


    

    
