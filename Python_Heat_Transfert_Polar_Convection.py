import numpy as np
import matplotlib.pyplot as plt
import time

path = "/Users/gilles/Desktop/TB/Output_python_v9/"                # path to the folder where to store the .csv files
#variables physiques 
rin      = 1                    #rayon interne
rout     = 0.8387*rin           #external radius
disc     = 50                   #155 est le nombre de cellule de la discontinuité du manteau
rho      = 1
rhog0_up = 1.0                  #densité du manteau supérieur
rhog0_dw = 1.0                  #densité du manteau inférieur
eta      = 1
eta_up   = 1                    #viscosité sup
eta_dw   = 1                    #viscosité infé
alph     = 1                    #coefficient d'expansion thermique
beta     = 0
Ttop     = -1                   #température surface du manteau
Tbot     = 1                    #température base du manteau
rhoCp    = 1                    #specific heat capacity
lam      = 1e-6                 #diffusivité thermique
ttime    = 1                    #temps total
s_ref    = 1e-10                  #Tau de référence afin d'avoir au début 1 avec la fraction pour power-law
n_exp    = 10                    #utile pour power-law, à faire varier pour trouver le début de l'activation de la power-law
rel      = 0.5
#variables numériques
#Angulaire (phi) 
nnp  = 181                         #nombre de noeuds
Lp   = 2*np.pi + 2*np.pi/(nnp-1)   #longueur angulaire plus un noeuds pour les conditions polaires
dphi = Lp/(nnp-1)                  #espacement angulaire
#Radiale (rad)
nnr  = 101                         #nombre de noeuds
Lr   = rout                        #longueur radiale
drad = Lr/(nnr-1)                  #espacemennt radiale
#densité
rhog_v = np.ones((nnr,nnp))         #densité initial                      #problème avec nnp.ones (pourtant documentation suivie)    
rhog_v [0:disc,:] = rhog0_dw        #partie inférieur avec différente densité 
rhog_v [disc:-1,:] = rhog0_up       #partie supérieur avec une densité
#viscosité avec stress angulaire et radiale
Eta = np.ones((nnr,nnp))           #matrice de viscosité initial pour les noeuds
Eta[0:disc,:] = eta_dw             #partie inférieur avec différente densité
Eta[disc:-1,:] = eta_up            #partie supérieur avec une densité
#print(Eta.shape)
#Contrainte de cisaillement
Eta_rp = np.zeros((nnr - 1, nnp - 1))  # matrice initial de cisaillement
# Staggered grid in 2 dimension rad and phi--------------------------------
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

#initialisation
P = (-rhog_v * radn)                   #pression lithostatique
taurr = np.zeros((nnr, nnp))
taupp = np.zeros((nnr, nnp))
taurp = np.zeros((nnr - 1, nnp - 1))

# Interpolation du centre aux noeuds
def c2n(A0):
    A1 = np.zeros((A0.shape[0]+1, A0.shape[1]))
    A1[:,:] = np.vstack((1.5*A0[0,:] - 0.5*A0[1,:], (A0[1:,:] + A0[:-1,:]) / 2, 1.5*A0[-1,:] - 0.5*A0[-2,:]))
    
    A2 = np.zeros((A1.shape[0], A1.shape[1]+1))
    A2[:,:] = np.hstack((1.5*A1[:,0].reshape(-1,1) - 0.5*A1[:,1].reshape(-1,1),(A1[:,1:] + A1[:,:-1]) / 2, 1.5*A1[:,-1].reshape(-1,1) - 0.5*A1[:,-2].reshape(-1,1)))
    return A2

#tictoc
def TicTocGenerator():
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)

Vr = np.zeros((nnr + 1, nnp))
Vp = np.zeros((nnr, nnp + 1))
#T = 0.1 * np.exp(-(radn - 3 / 2 * rin) ** 2 * 50 - phin ** 2 * 0. * radn ** 2) * np.sin(phin) + 1.0*(np.random.rand(nnr, nnp) - 0.5)
T = 0.1 * np.exp(-(radn - 3 / 2 * rin) ** 2 * 50 - phin ** 2 * 0. * radn ** 2) * np.sin(phin)
T_init = T   #ou T.copy()

#preprocessing
dtit = np.min([drad, dphi])**2 / (np.max(Eta) / np.max(rhog_v)) / 4 / 1e0  #timesteps
dtdiff = np.min([drad,dphi])**2 /(lam/rhoCp)/4/1e0
#ntime  = 21
ntime  = int(7e3)
#Conditions Mecaniques
bet    = 1e-2
#tol    = 5e-6
tol    = 5e-2
#Itération postprocessing
niter  = 1000
nout   = 3
k      = 1
finish = 4000
#Définiton de EtaL et Etapl
Eta_l = Eta
Eta_pl = Eta

#fonction pour les moyennes
def avr(U):
    U = 0.5 * (U[:-1, :] + U[1:, :])
    return U
def avp(U):
    U = 0.5 * (U[:, :-1] + U[:, 1:])
    return U
def av_N2C(U):
    U = avp(avr(U))
    return U


 # sets interactive mode on (needed for plot update)
plt.ion()
fig = plt.figure()

tic()
#boucle
for it in range(ntime):
     # Buoyancy
    rhofg = rhog_v * (1 - alph * T + beta * P)

    # Mechanics, must converge to 0
    error = 2 * tol
    iter = 0
    while error > tol or iter < 300:
        iter += 1
        Srr = -P[:, 1:-1] + taurr[:, 1:-1]
        Spp = -P[1:-1, :] + taupp[1:-1, :]
        deltaS = avr(taurr[:, 1:-1] - taupp[:, 1:-1])  #
        dVrdt = 1/rho * (np.diff(Srr, axis=0)/np.diff(radn[:, 1:-1], axis=0) +
                         (np.diff(taurp, axis=1)/np.diff(phic, axis=1) + deltaS) / radr[1:-1, 1:-1] -
                         (rhofg[:-1, 1:-1] + rhofg[1:, 1:-1])/2)
        dVpdt = 1/rho * (np.diff(taurp, axis=0)/np.diff(radc, axis=0) +
                         (np.diff(Spp, axis=1)/np.diff(phin[1:-1, :], axis=1) + 2 * avr(taurp)) / radp[1:-1, 1:-1])
        Vr[1:-1, 1:-1] += dVrdt * dtit
        Vp[1:-1, 1:-1] += dVpdt * dtit
        Vr[:, -1] = Vr[:, 1]
        Vr[:, 0]  = Vr[:, -2]
        Vp[:, -1] = Vp[:, 2]
        Vp[:, 0]  = Vp[:, -3] 
        
        # deviatoric strain rate
        Err = np.diff(Vr, axis=0)/np.diff(radr, axis=0)
        Epp = (np.diff(Vp, axis=1)/np.diff(phip, axis=1) + avr(Vr)) / radn
        Erp = ((np.diff(Vr[1:-1, :], axis=1)/np.diff(phir[1:-1, :], axis=1) - avr(Vp[:, 1:-1] )) / radc + np.diff(Vp[:, 1:-1], axis=0)/np.diff(radp[:, 1:-1], axis=0))/2
        # divergence of velocities
        divV = Err + Epp
        dPdt = -1/bet * divV
        P += dPdt * dtit

         #power-law viscosity
        tau2 = (0.5*(taurr**2+taupp**2)+(c2n(taurp)**2))**0.5

        if n_exp > 1:
            Eta_pl_it = Eta_pl
            Eta_pl    = Eta_l*((tau2/s_ref)**(1-n_exp))
            Eta_pl    = np.exp(np.log(Eta_pl)*rel+np.log(Eta_pl_it)*(1-rel))
            Eta       = 1/(1/Eta_l + 1/Eta_pl) 

        # not incompressible, compressible
        taurr = 2 * Eta * (Err - 1/3 * divV)
        taupp = 2 * Eta * (Epp - 1/3 * divV)
        Eta_rp = av_N2C(Eta)
        taurp = 2 * Eta_rp * Erp
        
        
        #print(tau2.shape)
       
        
        error = np.max([np.max(np.abs(dVrdt)), np.max(np.abs(dVpdt)), np.max(np.abs(dPdt))])
    print("time step = %d\terror = %10.3E\tTmin = %10.3E\tTmax = %10.3E\titerations = %d" % (it,error,T.min(),T.max(),iter)) # for print syntax, see: https://www.geeksforgeeks.org/python-output-formatting/
    if error > tol: 
        print('no convergence !! :\'(')
    
    # Thermics------------------------------------------------------------->we use the results from mechanicss to calculate thermics
    dt = np.min([drad / np.max(np.abs(Vr)) / 2 / 1e0, dtdiff])  # change of dt value if dt calculated is smaller than standard

    # Thermo
    dTdt_1 = - np.maximum(0, Vr[1:-2, 1:-1]) * np.diff(T[:-1, 1:-1], axis=0) / drad 
    dTdt_2 = - np.minimum(0, Vr[2:-1, 1:-1]) * np.diff(T[1:, 1:-1], axis=0) / drad 
    dTdt_3 = - np.maximum(0, Vp[1:-1, 1:-2]) * np.diff(T[1:-1, :-1], axis=1) / dphi / radn[1:-1, 1:-1]
    dTdt_4 = - np.minimum(0, Vp[1:-1, 2:-1]) * np.diff(T[1:-1, 1:], axis=1) / dphi / radn[1:-1, 1:-1]
    dTdt_5 = (np.diff(lam * np.diff(T[:, 1:-1], axis=0) / drad, axis=0) / drad) / rhoCp
    dTdt_6 = (lam * np.diff(avr(T[:, 1:-1]), axis=0) / drad / radn[1:-1, 1:-1]) / rhoCp
    dTdt_7 = (np.diff(lam * np.diff(T[1:-1, :], axis=1) / dphi, axis=1) / dphi / radn[1:-1, 1:-1] ** 2) / rhoCp
    dTdt = dTdt_1 + dTdt_2 + dTdt_3 + dTdt_4 + dTdt_5 + dTdt_6 + dTdt_7


    T[1:-1, 1:-1] = T[1:-1, 1:-1] + dTdt * dt

    # Polar boundaries conditions-------------------------------------------
    T[-1,:]            = Ttop                                                 # reset of the mantle surface T values
    T[0,  :]           = Tbot                                                 # reset of the mantle bottom T values
    T[:, -1]           = T[:, 1]                                             # reset of the last angular cell with the second cell (they are superimposed)
    T[:, 0]            = T[:, -2]                                           # reset of the second-to-last angular cell with the first cell (they are superimposed)
    
    # save variables--------------------------------------------------------
    #np.savetxt(path +"output_it_"+ str(it) +"_Vr.csv",Vr,delimiter=",")
    #np.savetxt(path +"output_it_"+ str(it) +"_Vp.csv",Vp,delimiter=",")
    #np.savetxt(path +"output_it_"+ str(it) +"_P.csv",P,delimiter=",")
    #np.savetxt(path +"output_it_"+ str(it) +"_T.csv",T,delimiter=",")
    #np.savetxt(path +"output_it_"+ str(it) +"_Tau.csv",tau2,delimiter=",")

    '''np.savetxt("output_it_"+ str(it) +"_Vr.csv",Vr,delimiter=",")
    np.savetxt("output_it_"+ str(it) +"_Vp.csv",Vp,delimiter=",")
    np.savetxt("output_it_"+ str(it) +"_P.csv",P,delimiter=",")
    np.savetxt("output_it_"+ str(it) +"_T.csv",T,delimiter=",")
    np.savetxt("output_it_"+ str(it) +"_Tau.csv",tau2,delimiter=",")'''
    # postprocessing=======================================================
    # polar to cartesian conversion-----------------------------------------
    x = radn*np.cos(phin)                                                      # conversion to cartesian coordinates x
    y = radn*np.sin(phin)                                                     # conversion to cartesian coordinates y
    Vx = Vr[1:,:]*np.cos(phir[1:,:]) - Vp[:,1:]*np.sin(phip[:,1:])         # conversion to cartesian coordinates Vx
    Vy = Vr[1:,:]*np.sin(phir[1:,:]) + Vp[:,1:]*np.cos(phip[:,1:])         # conversion to cartesian coordinates Vy

 
    if np.mod(it, nout) == 0: # postprocessing
        plt.clf() 
        ax1 = fig.add_subplot(2,2,1)
        pc1 = ax1.pcolor(x,y,T,cmap='jet',vmin=-0.1,vmax=0.1) 
        plt.colorbar(pc1,ax=ax1)
        ax1.set_aspect('equal','box')
        plt.title(it)
        ax1.set_title('T')

        ax1 = fig.add_subplot(2,2,2)
        pc1 = ax1.pcolor(xc,yc,taurp,cmap='jet')#,vmin=-0.1,vmax=0.1) 
        plt.colorbar(pc1,ax=ax1)
        ax1.set_aspect('equal','box')
        ax1.set_title('taurp')

        ax2 = fig.add_subplot(2,2,3)
        pc2 = ax2.pcolor(x,y,taupp,cmap='jet')#,vmin=-0.1,vmax=0.1) 
        plt.colorbar(pc2,ax=ax2)
        ax2.set_aspect('equal','box')
        ax2.set_title('taupp')

        ax2 = fig.add_subplot(2,2,4)
        pc2 = ax2.pcolor(x,y,taurr,cmap='jet')#,vmin=-0.1,vmax=0.1) 
        plt.colorbar(pc2,ax=ax2)
        ax2.set_aspect('equal','box')
        ax2.set_title('taurr')
        plt.pause(0.1)
        plt.draw()
toc()
 
plt.ioff()
plt.show()
