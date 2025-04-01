import matplotlib.pyplot as plt
from tqdm import tqdm
from utils import *
from plotting import animate_cart_single
from models_eq import CartSinglePendulum as CSP, PAR

A, B, C, LAB, f = CSP.A, CSP.B, CSP.C, CSP.LAB, CSP.f

FPS = 60.0

# simulation parameters
TS = 0.025 # controller time step / sampling time
STSM = 100 # simulation time step multiplier
DT = TS/STSM # sim time step
T = 50 # [s] simulation time
nt = int(T / DT) # number of time steps (simulation)
print(f"Simulation time steps: {nt} ({T} s)")

θ0 = deg2rad(4) # [rad] initial angle
x0 = -1 # [m] initial position

class PID():
    def __init__(self, kpx=-1, kix=0.0, kdx=2, kpθ=-30, kixθ=0.0, kdxθ=5):
    # def __init__(self, kpx=-0.01, kix=-0.0, kdx=+0.03, kpθ=-1, kixθ=0.0, kdxθ=0.1):
        self.kpx, self.kix, self.kdx = kpx, kix, kdx
        self.kpθ, self.kixθ, self.kdxθ = kpθ, kixθ, kdxθ
        self.exi, self.eθi = 0, 0
        self.x_prev, self.θ_prev = 0, 0
    
    def get_control(self, y, r):
        θ, x = y # system output
        rθ, rx = r # reference
        epθ = rθ - θ
        epx = rx - x
        self.exi += epx * TS
        self.eθi += epθ * TS
        edx = (x - self.x_prev) / TS
        edθ = (θ - self.θ_prev) / TS
        self.x_prev, self.θ_prev = x, θ

        u = self.kpx * epx + self.kix * self.exi + self.kdx * edx + \
            self.kpθ * epθ + self.kixθ * self.eθi + self.kdxθ * edθ

        return u
    

# MPC parameters
N = 100 # prediction horizon
Q = vec([[50,0,0,0],
        [0,100,0,0],
        [0,0,1,0],
        [0,0,0,1]])
R = vec([[0.1]])
S = Q
MAXU = 20 # [N] max control input

class MPC(): # linear mpc
    def __init__(self, A, B, Q, R, S, f, F, TS, N):
        self.Ac, self.Bc = A, B
        self.Q, self.R = Q, R
        self.N = N
        n, m = A.shape[0], B.shape[1] # state and input dimensions
        # discretize the system, assume C = identity, D = 0
        from scipy.signal import cont2discrete as c2d
        import numpy as np
        Ad, Bd, Cd, _, _ = c2d((A, B, np.eye(n), 0), TS, method='zoh') # zero order hold

        # create condensed/calligrafic matrices
        M = np.zeros((n, n))
        cms = create_condensed_matrices(Ad, Bd, Cd, M, Q, R, S, F, f, N)
        self.cAC, self.cBC, self.cMC, self.cQ, self.cR, self.cF, self.cf = cms

        self.Hqp = None

    def get_control(self, x, r):
        assert x.shape == (4,), f'x.shape: {x.shape}'
        assert r.shape == (4,), f'r.shape: {r.shape}'
        


# cart single pendulum
ss = zeros((nt, 4)) # [ θ, x, dθ , dx ]
ys = zeros((nt, 2)) # [ θ, x ]
us = zeros(nt) # [ u ]
ss[0] = vec([θ0, x0, 0.0, 0.0]) # initial state
ys[0] = C @ ss[0] # initial output

# integrate
u = 0 # control input
pid = PID()
for i in tqdm(range(1, nt)):
    if i % STSM == 0: # update control input
        u = pid.get_control(ys[i-1], vec([0.0, 0.0]))
        # u = 0
    ss[i] = step(f, ss[i-1], u, DT)
    ys[i] = C @ ss[i] # system output
    us[i] = u

# # plot
plt.figure(figsize=(10, 4))
plt.plot(ss)
plt.legend(LAB)
plt.title('Cart Single Pendulum')

# animation
anim = animate_cart_single(ss, us, DT, PAR.l1, figsize=(8,8))
plt.show()