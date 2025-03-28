## utils functions

## Select between numpy and torch
from numpy import array as vec, concatenate as cat, zeros, sin, cos # numpy version
# from torch import cat, tensor as vec, zeros, sin, cos # torch version

from time import time

# integration steps
def euler_step(f, s, u, dt): # integrate one step with explicit euler
    n = s.shape[-1]
    s, ds = s[:n//2], s[n//2:] # split the state into position and velocity
    dds = vec((f(*s, *ds, u)))
    ds = ds + dds * dt
    s = s + ds * dt
    return cat([s, ds])

def rk4_step(f, s, u, dt): # integrate one step with Runge-Kutta 4
    n = s.shape[-1]
    s, ds = s[:n//2], s[n//2:] # split the state into position and velocity
    dds1 = vec(f(*s, *ds, u)) #k1
    ds1 = ds
    ds2 = ds + dds1 * dt/2 # k2
    s2 = s + ds1 * dt/2
    dds2 = vec(f(*s2, *ds2, u))
    ds3 = ds + dds2 * dt/2 # k3
    s3 = s + ds2 * dt/2
    dds3 = vec(f(*s3, *ds3, u))
    ds4 = ds + dds3 * dt # k4
    s4 = s + ds3 * dt
    dds4 = vec(f(*s4, *ds4, u))
    ds_new = ds + (dds1 + 2*dds2 + 2*dds3 + dds4) * dt/6 # weighted average
    s_new = s + (ds1 + 2*ds2 + 2*ds3 + ds4) * dt/6 
    return cat([s_new, ds_new])

# def step(f,s,u,dt): return euler_step(f,s,u,dt)
def step(f,s,u,dt): return rk4_step(f,s,u,dt)