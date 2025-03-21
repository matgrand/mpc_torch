# test the integration
import numpy as np
import matplotlib.pyplot as plt
from models_eq import SinglePendulum as SP, DoublePendulum as DP, CartSinglePendulum as CSP, CartDoublePendulum as CDP, PAR
from plotting import animate_pendulum, animate_double_pendulum, animate_cart_single, animate_cart_double
from tqdm import tqdm
FPS = 60.0

def euler_step(f, s, u, dt): # integrate one step with explicit euler
    n = s.shape[-1]
    s, ds = s[:n//2], s[n//2:] # split the state into position and velocity
    dds = np.array((f(*s, *ds, u)))
    ds = ds + dds * dt
    s = s + ds * dt
    return np.concatenate([s, ds])

def rk4_step(f, s, u, dt): # integrate one step with Runge-Kutta 4
    n = s.shape[-1]
    s, ds = s[:n//2], s[n//2:] # split the state into position and velocity
    dds1 = np.array(f(*s, *ds, u)) #k1
    ds1 = ds.copy()
    s1 = s.copy()
    ds2 = ds + dds1 * dt/2 # k2
    s2 = s + ds1 * dt/2
    dds2 = np.array(f(*s2, *ds2, u))
    ds3 = ds + dds2 * dt/2 # k3
    s3 = s + ds2 * dt/2
    dds3 = np.array(f(*s3, *ds3, u))
    ds4 = ds + dds3 * dt # k4
    s4 = s + ds3 * dt
    dds4 = np.array(f(*s4, *ds4, u))
    ds_new = ds + (dds1 + 2*dds2 + 2*dds3 + dds4) * dt/6 # weighted average
    s_new = s + (ds1 + 2*ds2 + 2*ds3 + ds4) * dt/6 
    return np.concatenate([s_new, ds_new])

# def step(f,s,u,dt): return euler_step(f,s,u,dt)
def step(f,s,u,dt): return rk4_step(f,s,u,dt)

# simulation parameters
dt = 0.0001 # time step
T  = 15 # simulation time
nt = int(T / dt) # number of time steps

# single pendulum
ss = np.zeros((nt, 2))
ss[0] = [0.1, 0.0] # initial state
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
anim = animate_pendulum(ss, np.zeros(nt), dt, PAR.l1, figsize=(8,8), title='Single Pendulum')
plt.show()

# double pendulum
ss = np.zeros((nt, 4))
ss[0] = [0.1, 0.1, 0.0, 0.0] # initial state
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
anim = animate_double_pendulum(ss, np.zeros(nt), dt, PAR.l1, PAR.l2, figsize=(8,8), title='Double Pendulum')
plt.show()

# cart single pendulum
ss = np.zeros((nt, 4))
ss[0] = [0.1, 0.0, 0.0, 0.0] # initial state
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
anim = animate_cart_single(ss, np.zeros(nt), dt, PAR.l1, figsize=(8,8))
plt.show()

# cart double pendulum
ss = np.zeros((nt, 6))
ss[0] = [0.1, 0.1, 0.0, 0.0, 0.0, 0.0] # initial state
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
plt.show()
# animation
anim = animate_cart_double(ss, np.zeros(nt), dt, PAR.l1, PAR.l2, fps=FPS, figsize=(8,8))
plt.show()




