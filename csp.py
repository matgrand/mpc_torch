import matplotlib.pyplot as plt
from tqdm import tqdm
from utils import *
from plotting import animate_cart_single
from models_eq import CartSinglePendulum as CSP, PAR

A, B, C, LAB, f = CSP.A, CSP.B, CSP.C, CSP.LAB, CSP.f

FPS = 60.0

# simulation parameters
DT = 0.01 # controller time step
STSM = 10 # simulation time step multiplier
SDT = DT/STSM # sim time step
T = 50 # [s] simulation time
nt = int(T / SDT) # number of time steps (simulation)
print(f"Simulation time steps: {nt} ({T} s)")

θ0 = 0.1 # [rad] initial angle
x0 = -1 # [m] initial position

class PID():
    def __init__(self, kpx=-0.01, kix=-0.0, kdx=+0.03, kpθ=-1, kixθ=0.0, kdxθ=0.1):
        self.kpx, self.kix, self.kdx = kpx, kix, kdx
        self.kpθ, self.kixθ, self.kdxθ = kpθ, kixθ, kdxθ
        self.exi, self.eθi = 0, 0
        self.x_prev, self.θ_prev = 0, 0
    
    def get_control(self, y, r):
        θ, x = y # system output
        rθ, rx = r # reference
        epθ = rθ - θ
        epx = rx - x
        self.exi += epx * DT
        self.eθi += epθ * DT
        edx = (x - self.x_prev) / DT
        edθ = (θ - self.θ_prev) / DT
        self.x_prev, self.θ_prev = x, θ

        u = self.kpx * epx + self.kix * self.exi + self.kdx * edx + \
            self.kpθ * epθ + self.kixθ * self.eθi + self.kdxθ * edθ

        return u
    
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
    ss[i] = step(f, ss[i-1], u, SDT)
    ys[i] = C @ ss[i] # system output
    us[i] = u

# # plot
plt.figure(figsize=(10, 4))
plt.plot(ss)
plt.legend(LAB)
plt.title('Cart Single Pendulum')

# animation
anim = animate_cart_single(ss, us, SDT, PAR.l1, figsize=(8,8))
plt.show()