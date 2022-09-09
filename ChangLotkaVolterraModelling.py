# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 00:06:20 2021

@author: Edward Chang
"""

import numpy as np
import matplotlib.pyplot as plt


#                        Part 1: The Basic Outline

#Start with initial equation developed independently by Lotka and Volterra
#From https://web.ma.utexas.edu/users/davis/375/popecol/lec10/lotka.html
#Lotka, A. J. 1925. Elements of physical biology. Baltimore: Williams & Wilkins Co.
#Volterra, V. 1926. Variazioni e fluttuazioni del numero d'individui in specie animali conviventi. Mem. R. Accad. Naz. dei Lincei. Ser. VI, vol. 2.

H0 = 2.  #initial prey population density
P0 = 2.  #initial predator population density
r = 1.  #rate of prey population increase
a = 1.  #predation rate coefficient
b = 1.  #reproduction rate of predators per prey eaten
m = 1.  #predator mortality rate
dt = .005  #time step
t0 = 0.  #start time (can't change)
tf = 100.  #finishing time

H = np.array([])  #empty arrays to be used for plotting
P = np.array([])
t = np.array([])

def h(H,P):  #dH/dt
    return H*(r - a*P)
def p(H,P):  #dP/dt
    return P*(b*H - m)

fig,ax = plt.subplots(2,1)

while t0 < tf:  #runge-kutta 4th order method
    k1 = dt*h(H0,P0)  #change in prey population
    k2 = dt*h(H0 + k1/2,P0)
    k3 = dt*h(H0 + k2/2,P0)
    k4 = dt*h(H0 + k3,P0)
    Hf = H0 + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    H0 = Hf
    H = np.append(H,Hf)
    
    k1 = dt*p(H0,P0)  #change in predator population
    k2 = dt*p(H0,P0 + k1/2)
    k3 = dt*p(H0,P0 + k2/2)
    k4 = dt*p(H0,P0 + k3)
    Pf = P0 + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    P0 = Pf
    P = np.append(P,Pf)
    
    if t0 > 0. and t0 < .20:
        ax[0].plot(Hf,Pf,'bo')
    if t0 > .20 and t0 < .206:
        ax[0].plot(Hf,Pf,'ro')
    
    t = np.append(t,t0)
    t0 += dt

ax[0].plot(H,P,'magenta')
ax[0].set(title='Predators vs Prey', xlabel='H (prey pop. density)', ylabel='P (pred. pop. density)')
ax[0].set_aspect('equal')
ax[0].grid()
ax[0].set_xlim([0,np.amax(H) + 1])
ax[0].set_ylim([0,np.amax(P) + 1])

ax[1].plot(t,H,'b')
ax[1].plot(t,P,'r')
ax[1].set(title='Predators (Red) and Prey (Blue) vs time', xlabel='Time', ylabel='Population Density')
ax[1].grid()
ax[1].set_xlim(0,100)
if np.amax(P) >= np.amax(H):
    ax[1].set_ylim(0,np.amax(P) + 1)
else:
    ax[1].set_ylim(0,np.amax(H) + 1)
plt.tight_layout()
plt.show()

'''
#                       Part 2: Adding Complexity

#Let's add another prey and another predator

H10 = 1.2  #initial prey population density
H20 = 1.1
P10 = .9  #initial predator population density
P20 = .8
r1 = 1.  #rate of prey population increase
r2 = 1.
a1 = 1.  #predation rate coefficient
a2 = 1.
b1 = 1.  #reproduction rate of predators per prey eaten
b2 = 1.
m1 = 1.  #predator mortality rate
m2 = 1.
dt = .02  #time step
t0 = 0.  #start time (can't change)
tf = 1000.  #finishing time

H1 = np.array([])  #empty arrays to be used for plotting
H2 = np.array([])
P1 = np.array([])
P2 = np.array([])
t = np.array([])

hic = .1  #prey interaction coefficient  (for equal ic: 0 < ic < 19)
pic = .1  #predator interaction coefficient
cc12 = .5  #cross coefficient of prey 1 and predator 2  (for equal cc)
cc21 = .5  #cross coefficient of prey 2 and predator 1  (.32< cc <10.)
#                                                       (  equ =~ 1  )
def h1(H10,H20,P10,P20):  #dH1/dt
    return H10*(r1 - a1*P10 - hic*H20 - cc12*P20)
def h2(H10,H20,P10,P20):  #dH2/dt
    return H20*(r2 - a2*P20 - hic*H10 - cc21*P10)
def p1(H10,H20,P10,P20):  #dP1/dt
    return P10*(b1*H10 - m1 - pic*P20 + cc21*H20)
def p2(H10,H20,P10,P20):  #dP2/dt
    return P20*(b2*H20 - m2 - pic*P10 + cc12*H10)

H1 = np.append(H1,H10)  #initial values appended
H2 = np.append(H2,H20)
P1 = np.append(P1,P10)
P2 = np.append(P2,P20)
t = np.append(t,t0)

fig,ax = plt.subplots()

while t0 < tf:  #runge-kutta 4th order method
    k1 = dt*h1(H10,H20,P10,P20)  #change in prey 1
    k2 = dt*h1(H10 + k1/2,H20,P10,P20)
    k3 = dt*h1(H10 + k2/2,H20,P10,P20)
    k4 = dt*h1(H10 + k3,H20,P10,P20)
    H1f = H10 + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    H10 = H1f
    H1 = np.append(H1,H1f)
    
    k1 = dt*h2(H10,H20,P10,P20)  #change in prey 2
    k2 = dt*h2(H10,H20 + k1/2,P10,P20)
    k3 = dt*h2(H10,H20 + k2/2,P10,P20)
    k4 = dt*h2(H10,H20 + k3,P10,P20)
    H2f = H20 + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    H20 = H2f
    H2 = np.append(H2,H2f)
    
    k1 = dt*p1(H10,H20,P10,P20)  #change in predator 2
    k2 = dt*p1(H10,H20,P10 + k1/2,P20)
    k3 = dt*p1(H10,H20,P10 + k2/2,P20)
    k4 = dt*p1(H10,H20,P10 + k3,P20)
    P1f = P10 + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    P10 = P1f
    P1 = np.append(P1,P1f)
    
    k1 = dt*p2(H10,H20,P10,P20)  #change in predator 2
    k2 = dt*p2(H10,H20,P10,P20 + k1/2)
    k3 = dt*p2(H10,H20,P10,P20 + k2/2)
    k4 = dt*p2(H10,H20,P10,P20 + k3)
    P2f = P20 + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    P20 = P2f
    P2 = np.append(P2,P2f)
    
    t = np.append(t,t0)
    t0 += dt

ax.plot(t,H1,'b')
ax.plot(t,H2,'cyan')
ax.plot(t,P1,'r')
ax.plot(t,P2,'orange')
ax.set(title='Predators (Red) and Prey (Blue) vs time', xlabel='Time', ylabel='Population Density')
ax.grid()
ax.set_xlim(0,tf)
maxpopara = np.array([np.amax(H1),np.amax(H2),np.amax(P1),np.amax(P2)])
maxpop = np.amax(maxpopara)
ax.set_ylim(0,maxpop + .5)
plt.show()

#Makes the assumption that a species can come back from extinction
'''
'''
#                     Part 3: Population Limitation

#Problems with the Lotka/Volterra model of predator-prey relationships?
#Several assumptions are made
#Prey have no limits to their reproduction
#Predators have no limits to their consumption
#Both reproduce continuously

#This time, I'll implement a logistical growth term in order
#    to limit the population sizes and simulate overcrowding

H10 = .5  #initial prey population density
H20 = .4
P10 = .3  #initial predator population density
P20 = .2
r1 = 1.  #rate of prey population increase
r2 = 1.
a1 = 1.  #predation rate coefficient
a2 = 1.
b1 = 1.  #reproduction rate of predators per prey eaten
b2 = 1.
m1 = 1.  #predator mortality rate
m2 = 1.
dt = .02  #time step
t0 = 0.  #start time (can't change)
tf = 1000.  #finishing time

H1 = np.array([])  #empty arrays to be used for plotting
H2 = np.array([])
P1 = np.array([])
P2 = np.array([])
t = np.array([])

hic = .1  #prey interaction coefficient  (for equal ic: 0 < ic < 19)
pic = .1  #predator interaction coefficient
cc12 = .5  #cross coefficient of prey 1 and predator 2  (for equal cc)
cc21 = .5  #cross coefficient of prey 2 and predator 1  (.32< cc <10.)
#                                                       (  equ =~ 1  )
Kh1 = 1.
Kh2 = 1.
Kp1 = 1.
Kp2 = 1.

def h1(H10,H20,P10,P20):  #dH1/dt
    return H10*(r1 - a1*P10 - hic*H20 - cc12*P20)*(1 - H10/Kh1)
def h2(H10,H20,P10,P20):  #dH2/dt
    return H20*(r2 - a2*P20 - hic*H10 - cc21*P10)*(1 - H20/Kh2)
def p1(H10,H20,P10,P20):  #dP1/dt
    return P10*(b1*H10 - m1 - pic*P20 + cc21*H20)*(1 - P10/Kp1)
def p2(H10,H20,P10,P20):  #dP2/dt
    return P20*(b2*H20 - m2 - pic*P10 + cc12*H10)*(1 - P20/Kp2)

H1 = np.append(H1,H10)  #initial values appended
H2 = np.append(H2,H20)
P1 = np.append(P1,P10)
P2 = np.append(P2,P20)
t = np.append(t,t0)

fig,ax = plt.subplots()

while t0 < tf:  #runge-kutta 4th order method
    k1 = dt*h1(H10,H20,P10,P20)  #change in prey 1
    k2 = dt*h1(H10 + k1/2,H20,P10,P20)
    k3 = dt*h1(H10 + k2/2,H20,P10,P20)
    k4 = dt*h1(H10 + k3,H20,P10,P20)
    H1f = H10 + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    H10 = H1f
    H1 = np.append(H1,H1f)
    
    k1 = dt*h2(H10,H20,P10,P20)  #change in prey 2
    k2 = dt*h2(H10,H20 + k1/2,P10,P20)
    k3 = dt*h2(H10,H20 + k2/2,P10,P20)
    k4 = dt*h2(H10,H20 + k3,P10,P20)
    H2f = H20 + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    H20 = H2f
    H2 = np.append(H2,H2f)
    
    k1 = dt*p1(H10,H20,P10,P20)  #change in predator 2
    k2 = dt*p1(H10,H20,P10 + k1/2,P20)
    k3 = dt*p1(H10,H20,P10 + k2/2,P20)
    k4 = dt*p1(H10,H20,P10 + k3,P20)
    P1f = P10 + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    P10 = P1f
    P1 = np.append(P1,P1f)
    
    k1 = dt*p2(H10,H20,P10,P20)  #change in predator 2
    k2 = dt*p2(H10,H20,P10,P20 + k1/2)
    k3 = dt*p2(H10,H20,P10,P20 + k2/2)
    k4 = dt*p2(H10,H20,P10,P20 + k3)
    P2f = P20 + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    P20 = P2f
    P2 = np.append(P2,P2f)
    
    t = np.append(t,t0)
    t0 += dt

ax.plot(t,H1,'b')
ax.plot(t,H2,'cyan')
ax.plot(t,P1,'r')
ax.plot(t,P2,'orange')
ax.set(title='Predators (Red) and Prey (Blue) vs time', xlabel='Time', ylabel='Population Density')
ax.grid()
ax.set_xlim(0,tf)
maxpopara = np.array([np.amax(H1),np.amax(H2),np.amax(P1),np.amax(P2)])
maxpop = np.amax(maxpopara)
ax.set_ylim(0,maxpop + .5)
plt.show()
'''