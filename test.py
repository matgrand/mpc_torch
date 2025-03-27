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

# single pendulum
ss = zeros((nt, 2))
ss[0] = vec([0.1, 0.0]) # initial state
# integrate
for i in tqdm(range(1, nt)):
    ss[i] = step(SP.f, ss[i-1], 0, dt)
# plot
plt.figure(figsize=(10, 2))
plt.plot(ss[:,0], label='s')
plt.plot(ss[:,1], label='ds')
plt.legend()
plt.title('Single Pendulum')
# animation
anim = animate_pendulum(ss, zeros(nt), dt, PAR.l1, figsize=(8,8), title='Single Pendulum')
plt.show()

# double pendulum
ss = zeros((nt, 4))
ss[0] = vec([0.1, 0.1, 0.0, 0.0]) # initial state
# integrate
for i in tqdm(range(1, nt)):
    ss[i] = step(DP.f, ss[i-1], 0, dt)
# plot
plt.figure(figsize=(10, 4))
plt.plot(ss[:,0], label='s1')
plt.plot(ss[:,1], label='ds1')
plt.plot(ss[:,2], label='s2')
plt.plot(ss[:,3], label='ds2')
plt.legend()
plt.title('Double Pendulum')
# animation
anim = animate_double_pendulum(ss, zeros(nt), dt, PAR.l1, PAR.l2, figsize=(8,8), title='Double Pendulum')
plt.show()

# cart single pendulum
ss = zeros((nt, 4))
ss[0] = vec([0.01, 0.0, 0.0, 0.0]) # initial state
# integrate
for i in tqdm(range(1, nt)):
    ss[i] = step(CSP.f, ss[i-1], 0, dt)
# plot
plt.figure(figsize=(10, 4))
plt.plot(ss[:,0], label='s')
plt.plot(ss[:,1], label='ds')
plt.plot(ss[:,2], label='x')
plt.plot(ss[:,3], label='dx')
plt.legend()
plt.title('Cart Single Pendulum')
# animation
anim = animate_cart_single(ss, zeros(nt), dt, PAR.l1, figsize=(8,8))
plt.show()

# cart double pendulum
ss = zeros((nt, 6))
ss[0] = vec([0.01, -0.01, 0.0, 0.0, 0.0, 0.0]) # initial state
# integrate
for i in tqdm(range(1, nt)):
    ss[i] = step(CDP.f, ss[i-1], 0, dt)
# plot
plt.figure(figsize=(10, 6))
plt.plot(ss[:,0], label='s1')
plt.plot(ss[:,1], label='s2')
plt.plot(ss[:,2], label='s3')
plt.plot(ss[:,3], label='ds1')
plt.plot(ss[:,4], label='ds2')
plt.plot(ss[:,5], label='ds3')
plt.legend()
plt.title('Cart Double Pendulum')
# animation
anim = animate_cart_double(ss, zeros(nt), dt, PAR.l1, PAR.l2, fps=FPS, figsize=(8,8))
plt.show()




