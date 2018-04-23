#Final Project simulation code for ME 520
#Josh Stout

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

#constant definitions
k = 10 #W/mK thermal conductivity
rho = 1 #kg/m^3 density
Cp = 1 #J/kgK specific heat
R = 400 # Ohm electrical resistance
Rp = 0.1 #Ohm/K electrical resistance thermal coefficient
L = 1 #m length
D = 0.1 #m diameter
eps = 1 #radiative emissivity
I = 0.1 #A electrical current
sig = 5.67e-8 #W/m^2K^4 steffan-boltzmann

S = np.pi*D**2/4    #cross section area
gamma = Cp*rho*L**2/(np.pi**2*k) #theoretical thermal time constant
alpha = k/rho/Cp #theoretical thermal diffusivity

#spatial domain
nx = 25 #number cells
dx = L/nx #spacing
#time spacing, according to FTCS CFL
dt = dx*dx/2/alpha

#range of frequencies to simulate
w = np.linspace(1,81,20)

#cell center locations
x = np.linspace(0,L,nx,endpoint=False)+dx*0.5

#indexing
ind = np.arange(nx)
indp1 = np.roll(ind,-1)
indm1 = np.roll(ind,1)

#environment temperature
T0 = 275
#convenient coefficient groupings:
g = 4*I**2/(L*np.pi**2*D**2*rho*Cp)
C1 = 4*eps*sig/D
beta = alpha*dt/dx/dx

#rms voltage storage for 3 omega signal at each frequency
V3rms = np.zeros(len(w))

#grouping of load terms for numerical linearization
def G(tbg,Tbg,freq):
    #Tbg is best guess temperature, t is time value at best guess, freq is load current frequency
    return g*np.sin(2*np.pi*freq*tbg)**2*(R+Rp*Tbg) \
        - C1*(Tbg*Tbg*Tbg*Tbg+4*Tbg*Tbg*Tbg*T0+6*Tbg*Tbg*T0*T0+4*Tbg*T0*T0*T0)

#temperature rate of change of load terms
def dGdT(tbg,Tbg,freq):
    #rate of change of G with respect to temperature
    return g*np.sin(2*np.pi*freq*tbg)**2*Rp \
        - C1*(4*Tbg*Tbg*Tbg+12*Tbg*Tbg*T0+12*Tbg*T0*T0+4*T0*T0*T0)

#plt.ion()
#loop through each frequency
for k in range(len(w)):
#for k in range(19,20):
    wi = w[k] #current frequency, in Hz
    #reset intial time
    t = 0 #start time
    #final time: long enough for 5 load oscillations, 10 temperature oscillations, 15 3 \omega voltage oscillations
    oc = 5 #number of load oscillations to simulate
    tf = oc/wi  #stop time
    nt = 0 #time loop count
    nf = int(tf/dt) #number of timesteps required
    ntl = nf - int(1/wi/dt) #timestep to start recording on for last oscillation
    n_out = 100 #plot output step interval
    tol = 1e-6 #iterative residual tolerance level
    T_0 = np.zeros(nx)  #old time Temperature
    T_1 = np.copy(T_0)  #new time Temperature
    T_1s = np.copy(T_0) #new Temperature best guess
    T_max = np.zeros(nf+1) #storage for maximum temperature for each timestep
    T_ave = np.zeros(nf+1) #storage for average temperature for each timestep
    t_plot = np.linspace(0,tf,nf+1) #time array for plotting

    #storage for current best guess G and derivative
    G_s = np.zeros(nx)
    dG_s = np.zeros(nx)
    #spatial integration of temperature variation, for 1 omega worth of timesteps
    dR = np.zeros(nf-ntl)
    #resulting time dependence of voltage for last omega worth of steps
    V = np.zeros(nf-ntl)
    #3 omega fit to full voltage response
    V3 = np.zeros(nf-ntl)
    t_last = np.zeros(nf-ntl)

#    fig1 = plt.figure(figsize=(16,8))
#    ax1 = plt.subplot(111)
#    dat = ax1.plot(x,T_1,label='T(x,t)')
#    ax1.set_ylim([0,1])
#    ax1.set_xlabel('x (m)')
#    ax1.set_ylabel('T (C)')
#    ax1.set_title('For $\omega$ = %1.0f, \t T(x,t = %1.2f)'%(wi,t))

    T_max[0] = max(T_0)
    T_ave[0] = np.average(T_0)
    #simulation loop, while current count less than max count
    while nt < nf:
        t += dt
        nt +=1
        errmax = 1 #reset iterative error to arbitrary>tol
        #update T_0 with last loop's T_1
        T_0 = np.copy(T_1)
        while abs(errmax) > tol:
            #update current best guess with last iteration's T_1
            T_1s = np.copy(T_1)
            #update best guess G
            G_s = G(t,T_1s,wi)
            dG_s = dGdT(t,T_1s,wi)

            #solve for next best guess of each cell
            for i in range(nx):
                ip1 = indp1[i]
                im1 = indm1[i]
                if i == 0: #left boundary
                    T_1[i] = 1/(1-dt*dG_s[i])*(T_0[i]+beta*(T_0[ip1]-3*T_0[i]) +dt*(G_s[i]-dG_s[i]*T_1s[i]) )
                elif i == nx-1: #right boundary
                    T_1[i] = 1/(1-dt*dG_s[i])*(T_0[i]+beta*(-3*T_0[i]+T_0[im1]) +dt*(G_s[i]-dG_s[i]*T_1s[i]) )
                else: #non boundary cells
                    T_1[i] = 1/(1-dt*dG_s[i])*(T_0[i]+beta*(T_0[ip1]-2*T_0[i]+T_0[im1]) +dt*(G_s[i]-dG_s[i]*T_1s[i]) )

            err_array = T_1-T_1s
            errmax = max(err_array)
    #        print(errmax)
#            print("maximum error between T_1 and T_1s: %1.6f"%(errmax))
        #once errormax less than tol, leave iteration loop
        T_max[nt] = max(T_1)
        T_ave[nt] = np.average(T_1)

        #if currently in last load cycle, record voltage data for fitting
        if nt > ntl:
            #step in last cycle
            i2 = nt-ntl-1
            #calculate resistance variation for last load cycle
            #cell-centered riemann sum
            for i in range(nx):
                dR[i2] += T_1[i]*dx
            t_last[i2] = t
            dR[i2] = Rp/L*dR[i2]
            V[i2] = I*np.sin(2*np.pi*wi*t)*(R+dR[i2])

    #extract 3 omega component of voltage for last cycle
    guess_phase = np.pi/4
    guess_amp = 1
    #first, fit the primary \omega frequency
    optimize_func = lambda x: x[0]*np.sin(2*np.pi*wi*t_last+x[1])-V
    est_amp,est_phase = leastsq(optimize_func,[guess_amp, guess_phase])[0]
    #then subtract this from the total to retain just the 3\omega component
    V3t = V-est_amp*np.sin(2*np.pi*wi*t_last+est_phase)
    #optimize the amplitude and offset for this component
    opt2 = lambda x: x[0]*np.sin(2*np.pi*3*wi*t_last+x[1])-V3t
    amp2,phase2 = leastsq(opt2,[.1,np.pi/4])[0]
    V3 = amp2*np.sin(3*2*np.pi*wi*t_last+phase2)
    V3rms[k] = abs(amp2)/2**0.5*3

    print('3 w amplitude fit: %1.4f, \t 3 w phase fit: %1.4f'%(amp2,phase2))
        #every n_out steps, update plot
#        if nt%n_out ==0:
#            dat[0].set_ydata(T_1)
#            ax1.set_title('For $\omega$ = %1.0f, \t T(x,t = %1.2f)'%(wi,t))
#            fig1.canvas.draw()
#            fig1.canvas.flush_events()

#    dat[0].set_ydata(T_1)
#    ax1.set_title('For $\omega$ = %1.0f, \t T(x,t = %1.2f)'%(wi,t))
#    fig1.canvas.draw()
#    fig1.canvas.flush_events()

#    plt.ioff()

    fig2 = plt.figure(figsize=(8,4))
    ax2 = plt.subplot(111)
    ax2.plot(t_plot,T_max,t_plot,T_ave,t_plot,I*np.sin(2*np.pi*wi*t_plot))

    ax2.legend(['Maximum Temperature','Average Temperature','Load Current'],loc=3)
    ax2.set_title('$3\omega$ measurement Temperatures vs time, $\omega$ = %1.0f'%(wi))
    ax2.set_xlabel('t (s)')
    ax2.set_ylabel('T (C)')

#    fig3 = plt.figure(figsize=(16,8))
#    ax3 = plt.subplot(111)
#    ax3.plot(t_last,V)

#    ax3.set_title('Final Paper Temperatures Full Voltage oscillation in last cycle, $\omega$ = %1.0f'%(wi))
#    ax3.set_xlabel('t (s)')
#    ax3.set_ylabel('V (V)')

#    fig4 = plt.figure(figsize=(16,8))
#    ax4 = plt.subplot(111)
#    ax4.plot(t_last,V3t, t_last, V3)

#    ax4.set_title('Final Paper Temperatures 3 $\omega$ Voltage oscillation in last cycle, $\omega$ = %1.0f'%(wi))
#    ax4.set_xlabel('t (s)')
#    ax4.set_ylabel('$V_{3\omega}$ (V)')

#now fit the V3rms value to the theoretical curve:

guess_k = 9
guess_gamma = .01

opt3 = lambda x: 4*(I/np.sqrt(2))**3*L*R*Rp/(np.pi**4*x[0]*S*np.sqrt(1+(2*2*np.pi*w*x[1])**2)) - V3rms
k_ap, gamma_ap = leastsq(opt3,[guess_k,guess_gamma])[0]
#print(k_ap,gamma_ap)
Cp_ap = np.pi**2*gamma_ap*k_ap/rho/L**2
print('Apparent thermal conduction: %1.4f, \t Apparent Time Constant: %1.4f'%(k_ap, gamma_ap))
print('Apparent Specific Heat: %1.4f'%(Cp_ap))
V3rmsfit = 4*(I/np.sqrt(2))**3*L*R*Rp/(np.pi**4*k_ap*S*np.sqrt(1+(2*2*np.pi*w*gamma_ap)**2))


fig5 = plt.figure(figsize=(8,4))
ax5 = plt.subplot(111)
ax5.plot(w,V3rms,w,V3rmsfit)
ax5.legend(['Computed V3rms','Fitted V3rms'])
ax5.set_title('$V_{3\omega ,rms}$ vs Frequency, $k_{fit}$ = %1.4f, $Cp_{fit}$ = %1.4f'%(k_ap,Cp_ap))
ax5.set_xlabel('Frequency, Hz')
ax5.set_ylabel('$V_{3\omega ,rms}$')


#plt.ioff()
plt.show()
