# test the integration
import matplotlib.pyplot as plt
from models_eq import SinglePendulum as SP, DoublePendulum as DP, CartSinglePendulum as CSP, CartDoublePendulum as CDP, PAR
from plotting import animate_pendulum, animate_double_pendulum, animate_cart_single, animate_cart_double
from utils import *
from tqdm import tqdm

FPS = 60.0


# simulation parameters
dt = 0.001 # time step
T  = 15 # simulation time
nt = int(T / dt) # number of time steps

uλ = lambda i: 0.2*sin(0.3*2*π*T*i/nt+π/2) # control input
# uλ = lambda i: 0 # control input

us = vec([uλ(i) for i in range(nt)]) # control input

# single pendulum
ss = zeros((nt, 2))
ss[0] = vec([0.1, 0.0]) # initial state
# integrate
for i in tqdm(range(1, nt)):
    ss[i] = step(SP.f, ss[i-1], uλ(i), dt)
# plot
plt.figure(figsize=(10, 2))
plt.plot(ss)
plt.plot(us)
plt.legend(SP.LAB + ['u'])
plt.title('Single Pendulum')
# animation
anim = animate_pendulum(ss, us, dt, PAR.l1, figsize=(8,8), title='Single Pendulum')
plt.show()

# double pendulum
ss = zeros((nt, 4))
ss[0] = vec([0.1, 0.1, 0.0, 0.0]) # initial state
# integrate
for i in tqdm(range(1, nt)):
    ss[i] = step(DP.f, ss[i-1], uλ(i), dt)
# plot
plt.figure(figsize=(10, 4))
plt.plot(ss)
plt.plot(us)
plt.legend(DP.LAB + ['u'])
plt.title('Double Pendulum')
# animation
anim = animate_double_pendulum(ss, us, dt, PAR.l1, PAR.l2, figsize=(8,8), title='Double Pendulum')
plt.show()

# cart single pendulum
ss = zeros((nt, 4))
ss[0] = vec([0.01, 0.0, 0.0, 0.0]) # initial state
# integrate
for i in tqdm(range(1, nt)):
    ss[i] = step(CSP.f, ss[i-1], uλ(i), dt)
# plot
plt.figure(figsize=(10, 4))
plt.plot(ss)
plt.plot(us)
plt.legend(CSP.LAB + ['u'])
plt.title('Cart Single Pendulum')
# animation
anim = animate_cart_single(ss, us, dt, PAR.l1, figsize=(8,8))
plt.show()

# cart double pendulum
ss = zeros((nt, 6))
ss[0] = vec([0.01, -0.01, 0.0, 0.0, 0.0, 0.0]) # initial state
# integrate
for i in tqdm(range(1, nt)):
    ss[i] = step(CDP.f, ss[i-1], uλ(i), dt)
# plot
plt.figure(figsize=(10, 6))
plt.plot(ss)
plt.plot(us)
plt.legend(CDP.LAB + ['u'])
plt.title('Cart Double Pendulum')
# animation
anim = animate_cart_double(ss, us, dt, PAR.l1, PAR.l2, fps=FPS, figsize=(8,8))
plt.show()




