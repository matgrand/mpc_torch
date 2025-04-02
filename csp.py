import matplotlib.pyplot as plt
from tqdm import tqdm
from utils import *
import numpy as np
from plotting import animate_cart_single
from models_eq import CartSinglePendulum as sys, PAR

# A, B, C, LAB, f = CSP.A, CSP.B, CSP.C, CSP.LAB, CSP.f

FPS = 60.0

# simulation parameters
TS = 0.025 # controller time step / sampling time
STSM = 100 # simulation time step multiplier
DT = TS/STSM # sim time step
T = 10 # [s] simulation time
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
N = 30 # prediction horizon
Q = vec([[100,0,0,0],
        [0,50,0,0],
        [0,0,1,0],
        [0,0,0,1]])
R = vec([[0.0]])
S = Q
# from scipy.linalg import solve_discrete_are as dare
# S = dare(sys.A, sys.B, Q, R) # terminal cost
MAXU = 10 # [N] max control input
F = vec([[1], [-1]]) # inequality constraints matrix
f = vec([MAXU, MAXU]) # inequality constraints vector


class MPC(): # linear mpc
    def __init__(self, A=sys.A, B=sys.B, Q=Q, R=R, S=S, f=f, F=F, TS=TS, N=N):
        self.Ac, self.Bc = A, B
        self.Q, self.R = Q, R
        self.N = N
        n, m = A.shape[0], B.shape[1] # state and input dimensions
        assert m == 1
        # discretize the system, assume C = identity, D = 0
        from scipy.signal import cont2discrete as c2d
        import numpy as np; np.set_printoptions(precision=3, suppress=True)
        Ad, Bd, Cd, _, _ = c2d((A, B, np.eye(n), 0), TS, method='zoh') # zero order hold
        print(f'Ad: {Ad.shape}, Bd: {Bd.shape}, Cd: {Cd.shape}')

        # create condensed/calligrafic matrices
        cms = create_condensed_matrices(Ad, Bd, Cd, np.zeros((n, n)), Q, R, S, F, f, N)
        self.cA, self.cB, self.cM, self.cQ, self.cR, self.cF, self.cf = cms
        self.n, self.m = n, m 
        assert self.cA.shape == (N*n, n), f'cA.shape: {self.cA.shape}'
        assert self.cB.shape == (N*n, N*m), f'cB.shape: {self.cB.shape}'
        assert self.cQ.shape == (N*n, N*n), f'cQ.shape: {self.cQ.shape}'
        assert self.cR.shape == (N*m, N*m), f'cR.shape: {self.cR.shape}'
        # print(f'cA: {self.cA.shape} \n{self.cA}')

    def get_control(self, x, r):
        # remove self for convenience
        cA, cB, cM, cQ, cR, cF, cf, n, m = self.cA, self.cB, self.cM, self.cQ, self.cR, self.cF, self.cf, self.n, self.m
        assert x.shape == (n,), f'x.shape: {x.shape}'
        assert r.shape == (n,), f'r.shape: {r.shape}'
        x = x.reshape(-1, 1) # make sure it's a column vector

        Yr = np.kron(np.ones(N), r).reshape(-1, 1) # reference
        assert Yr.shape == (N*n,1), f'Yr.shape: {Yr.shape}'

        # quadratic cost
        Hqp = cB.T @ cQ @ cB + cR # quadratic cost, [N*m, N*m]
        Hqp = (Hqp + Hqp.T) / 2 # make sure it's symmetric
        # assert Hqp.shape == (N*m, N*m), f'Hqp.shape: {Hqp.shape}'
        # fqp1 = cA @ x - Yr
        # assert fqp1.shape == (N*n, 1), f'fqp1.shape: {fqp1.shape}'
        # fqp2 = cB.T @ cQ
        # assert fqp2.shape == (N, N*n), f'fqp2.shape: {fqp2.shape}'
        # fqp3 = fqp2 @ fqp1
        # assert fqp3.shape == (N*m, 1), f'fqp3.shape: {fqp3.shape}'
        fqp = cB.T @ cQ @ (cA @ x - Yr)
        assert fqp.shape == (N*m, 1), f'fqp.shape: {fqp.shape}'
        fqp = fqp.flatten() # make sure it's 1D
        # constraints
        Aqp, bqp = self.cF, self.cf.flatten() # constraints
        # solve the quadratic program
        cU = quadprog(Hqp, fqp, Aqp, bqp)
        return cU[0] # return the first control input


# cart single pendulum
ss = zeros((nt, 4)) # [ θ, x, dθ , dx ]
ys = zeros((nt, 2)) # [ θ, x ]
us = zeros(nt) # [ u ]
ss[0] = vec([θ0, x0, 0.0, 0.0]) # initial state
ys[0] = sys.C @ ss[0] # initial output

# integrate
u = 0 # control input
pid = PID()
mpc = MPC()
for i in tqdm(range(1, nt)):
    if i % STSM == 0: # update control input
        # u = pid.get_control(ys[i-1], vec([0.0, 0.0]))
        u = mpc.get_control(ss[i-1], vec([0.0, 0.0, 0.0, 0.0]))
        assert -MAXU*1.01 <= u <= MAXU*1.01, f'u: {u}, MAXU: {MAXU}'
        # u = 0
    ss[i] = step(sys.f, ss[i-1], u, DT)
    ys[i] = sys.C @ ss[i] # system output
    us[i] = u

# # plot
plt.figure(figsize=(10, 4))
plt.plot(ss)
plt.legend(sys.LAB)
plt.title('Cart Single Pendulum')

# animation
anim = animate_cart_single(ss, us, DT, PAR.l1, figsize=(8,8))
plt.show()