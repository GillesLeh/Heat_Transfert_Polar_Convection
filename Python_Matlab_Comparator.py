import numpy as np
import matplotlib.pyplot as plt

# chose which file
it_start    = 0         # first time step
it_end      = 175        # last time step
it_step     = 1

its = np.arange(it_start,it_end+1,it_step)

# ...needed for plotting...
rin  = 1                            # rayon interne
rout = 0.8387*rin                   # external radius
nnp  = 181                          # nombre de noeuds
Lp   = 2*np.pi + 2*np.pi/(nnp-1)    # longueur angulaire plus un noeuds pour les conditions polaires
dphi = Lp/(nnp-1)                   # espacement angulaire
nnr  = 101                          # nombre de noeuds
Lr   = rout                         # longueur radiale
drad = Lr/(nnr-1)                   # espacemennt radiale
radn1, phin1 = np.meshgrid(np.arange(-Lr/2, (Lr/2)+drad, drad) + (Lr/2 + rin), np.arange(-Lp/2, (Lp/2)+dphi/2, dphi)) 
radc1, phic1 = np.meshgrid(np.arange(-(Lr-drad)/2, (Lr-drad)/2+drad, drad) + (Lr/2 + rin), np.arange(-(Lp-dphi)/2, (Lp-dphi)/2+dphi, dphi))
radr1, phir1 = np.meshgrid(np.arange(-(Lr+drad)/2, (Lr+drad)/2+drad, drad) + (Lr/2 + rin), np.arange(-Lp/2, Lp/2+dphi/2, dphi))
radp1, phip1 = np.meshgrid(np.arange(-Lr/2, (Lr/2)+drad, drad) + (Lr/2 + rin), np.arange(-(Lp+dphi)/2, (Lp+dphi)/2+dphi, dphi))
radn = np.transpose(radn1)
phin = np.transpose(phin1)
radc = np.transpose(radc1)
phic = np.transpose(phic1)
radr = np.transpose(radr1)
phir = np.transpose(phir1)
radp = np.transpose(radp1)
phip = np.transpose(phip1)
x    = radn*np.cos(phin)               # conversion to cartesian coordinates x
y    = radn*np.sin(phin)    
xc   = radc*np.cos(phic)               # conversion to cartesian coordinates x
yc   = radc*np.sin(phic)    
xvr  = radr*np.cos(phir)               # conversion to cartesian coordinates x
yvr  = radr*np.sin(phir)    
xvp  = radp*np.cos(phip)               # conversion to cartesian coordinates x
yvp  = radp*np.sin(phip)    

plt.ion()
fig = plt.figure()

for it in its:
    # import
    Vr_python = np.genfromtxt("output_it_"+ str(it) +"_Vr.csv"  , delimiter=',')
    Vp_python = np.genfromtxt("output_it_"+ str(it) +"_Vp.csv"  , delimiter=',')
    P_python = np.genfromtxt("output_it_"+ str(it) +"_P.csv"  , delimiter=',')
    T_python = np.genfromtxt("output_it_"+ str(it) +"_T.csv"  , delimiter=',')
    Tau2_python = np.genfromtxt("output_it_"+ str(it) +"_Tau.csv"  , delimiter=',')
    T2_python = np.genfromtxt("output_it_"+ str(it) +"_T2.csv"  , delimiter=',')

    Vr_matlab = np.genfromtxt("output_it_"+ str(it) +"_Vr_M.csv", delimiter=',')
    Vp_matlab = np.genfromtxt("output_it_"+ str(it) +"_Vp_M.csv", delimiter=',')
    P_matlab = np.genfromtxt("output_it_"+ str(it) +"_P_M.csv", delimiter=',')
    T_matlab = np.genfromtxt("output_it_"+ str(it) +"_T_M.csv", delimiter=',')

    plt.clf() 

    # plot Vr
    ax1 = fig.add_subplot(3,5,1)
    pc1 = ax1.pcolor(xvr,yvr,Vr_python,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc1,ax=ax1)
    ax1.set_aspect('equal','box')
    ax1.set_ylabel('Python')
    ax1.set_title('Vr')

    ax2 = fig.add_subplot(3,5,6)
    pc2 = ax2.pcolor(xvr,yvr,Vr_matlab,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc2,ax=ax2)
    ax2.set_aspect('equal','box')
    ax2.set_ylabel('Matlab')

    ax3 = fig.add_subplot(3,5,11)
    pc3 = ax3.pcolor(xvr,yvr,Vr_python-Vr_matlab,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc3,ax=ax3)
    ax3.set_aspect('equal','box')
    ax3.set_ylabel('Python-Matlab')

    # plot Vp
    ax1 = fig.add_subplot(3,5,2)
    pc1 = ax1.pcolor(xvp,yvp,Vp_python,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc1,ax=ax1)
    ax1.set_aspect('equal','box')
    ax1.set_title('Vp')

    ax2 = fig.add_subplot(3,5,7)
    pc2 = ax2.pcolor(xvp,yvp,Vp_matlab,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc2,ax=ax2)
    ax2.set_aspect('equal','box')

    ax3 = fig.add_subplot(3,5,12)
    pc3 = ax3.pcolor(xvp,yvp,Vp_python-Vp_matlab,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc3,ax=ax3)
    ax3.set_aspect('equal','box')

    # plot P
    ax1 = fig.add_subplot(3,5,3)
    pc1 = ax1.pcolor(x,y,P_python,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc1,ax=ax1)
    ax1.set_aspect('equal','box')
    ax1.set_title('P')

    ax2 = fig.add_subplot(3,5,8)
    pc2 = ax2.pcolor(x,y,P_matlab,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc2,ax=ax2)
    ax2.set_aspect('equal','box')

    ax3 = fig.add_subplot(3,5,13)
    pc3 = ax3.pcolor(x,y,P_python-P_matlab,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc3,ax=ax3)
    ax3.set_aspect('equal','box')

    # plot T
    ax1 = fig.add_subplot(3,5,4)
    pc1 = ax1.pcolor(x,y,T_python,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc1,ax=ax1)
    ax1.set_aspect('equal','box')
    ax1.set_title('T')

    ax2 = fig.add_subplot(3,5,9)
    pc2 = ax2.pcolor(x,y,T_matlab,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc2,ax=ax2)
    ax2.set_aspect('equal','box')

    ax3 = fig.add_subplot(3,5,14)
    pc3 = ax3.pcolor(x,y,T_python-T_matlab,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc3,ax=ax3)
    ax3.set_aspect('equal','box')


    #Plot Tau2
    ax1 = fig.add_subplot(3,5,5)
    pc1 = ax1.pcolor(x,y,Tau2_python,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc1,ax=ax1)
    ax1.set_aspect('equal','box')
    ax1.set_title('Tau2')

    #plot T2
    ax2 = fig.add_subplot(3,5,10)
    pc2 = ax2.pcolor(x,y,T2_python,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc2,ax=ax2)
    ax2.set_aspect('equal','box')
    ax2.set_title('T2')


    '''ax3 = fig.add_subplot(3,5,15)
    pc3 = ax3.pcolor(x,y,Tau2_python,cmap='jet')#,vmin=-0.1,vmax=0.1)
    plt.colorbar(pc3,ax=ax3)
    ax3.set_aspect('equal','box')'''

    #fig.tight_layout()
    fig.suptitle("time step, it = "+ str(it))
    plt.pause(0.1)
    plt.draw()
    plt.show()

plt.ioff()
plt.show()
