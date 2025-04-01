## utils functions

## Select between numpy and torch
from numpy import array as vec, concatenate as cat, zeros, sin, cos, pi as π # numpy version
# from torch import cat, tensor as vec, zeros, sin, cos, pi as π # torch version

from time import time

# integration steps
def euler_step(f, s, u, dt): # integrate one step with explicit euler
    n = s.shape[-1]
    s, ds = s[:n//2], s[n//2:] # split the state into position and velocity
    dds = vec((f(*s, *ds, u)))
    s = s + ds * dt
    ds = ds + dds * dt
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
    s_new = s + (ds1 + 2*ds2 + 2*ds3 + ds4) * dt/6 
    ds_new = ds + (dds1 + 2*dds2 + 2*dds3 + dds4) * dt/6 # weighted average
    return cat([s_new, ds_new])

def step(f,s,u,dt): return euler_step(f,s,u,dt)
# def step(f,s,u,dt): return rk4_step(f,s,u,dt)

# linear system in the form x' = Ax + Bu
def lin_step_euler(A, B, s, u, dt): # integrate one step with explicit euler
    return s + (A @ s + B @ u) * dt

def lin_step_rk4(A, B, s, u, dt): # integrate one step with Runge-Kutta 4
    k1 = A @ s + B @ u
    k2 = A @ (s + dt/2 * k1) + B @ u
    k3 = A @ (s + dt/2 * k2) + B @ u
    k4 = A @ (s + dt * k3) + B @ u
    return s + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

def lstep(A,B,s,u,dt): return lin_step_euler(A,B,s,u,dt)
# def lstep(A,B,s,u,dt): return lin_step_rk4(A,B,s,u,dt)







## MATLAB like functions

def deg2rad(deg): return deg * π / 180 # degrees to radians
def rad2deg(rad): return rad * 180 / π # radians to degrees

# from scipy.signal import cont2discrete as c2d 
# # A, B, C, D, Ts = cont2discrete((A,B,C,D), Ts, method='zoh') # zero order hold


# creating the condensed matrices

def create_condensed_matrices(A, B, C, M, Q, R, S, F, f, N):
    """
    Create condensed matrices for Model Predictive Control.
    
    Parameters:
    -----------
    A, B, C : ndarray
        State space matrices
    M : ndarray
        Disturbance matrix
    Q, R, S : ndarray
        State, input, terminal cost matrices
    F, f : ndarray
        Inequality constraints
    N : int
        Prediction horizon
        
    Returns:
    --------
    tuple:
        (cAC, cBC, cMC, cQ, cR, cF, cf)
        Condensed matrices for the MPC problem
    """
    import numpy as np
    from numpy import kron, eye, zeros, ones
    # Get dimensions
    n, m = A.shape[0], B.shape[1]
    p, q = C.shape[0], M.shape[1]
    
    # Check dimensions
    assert A.shape == (n, n), f'A.shape: {A.shape}'
    assert B.shape == (n, m), f'B.shape: {B.shape}'
    assert C.shape == (p, n), f'C.shape: {C.shape}'
    assert M.shape == (n, q), f'M.shape: {M.shape}'
    assert Q.shape == (n, n), f'Q.shape: {Q.shape}'
    assert R.shape == (m, m), f'R.shape: {R.shape}'
    assert S.shape == (n, n), f'S.shape: {S.shape}'

    # Create standard condensed matrices
    cQ = kron(eye(N), Q)
    cR = kron(eye(N), R)
    cF = kron(eye(N), F)
    cf = kron(ones((N, 1)), f)

    # Pre-allocate matrices
    cAC = zeros((N*p, n))
    cBC = zeros((N*p, N*m))
    cMC = zeros((N*p, N*q))

    # Fill matrices
    for i in range(1, N+1):
        row_idx = (i-1)*p
        
        # cAC
        cAC[row_idx:row_idx+p, :] = C @ np.linalg.matrix_power(A, i)
        
        # cBC and cMC
        for j in range(1, N+1):
            col_idx_B = (j-1)*m
            col_idx_M = (j-1)*q
            
            if j <= i:
                A_power = np.linalg.matrix_power(A, i-j)
                cBC[row_idx:row_idx+p, col_idx_B:col_idx_B+m] = C @ A_power @ B
                cMC[row_idx:row_idx+p, col_idx_M:col_idx_M+q] = C @ A_power @ M
    
    # Add terminal cost
    # Update cQ with terminal cost
    cQ = np.block([
        [cQ[:(N-1)*p, :(N-1)*p], zeros(((N-1)*p, n))],
        [zeros((n, (N-1)*p)), S]
    ])
    
    # Create new matrices with terminal cost
    cAC_new = zeros((N*p, n))
    cAC_new[:(N-1)*p, :] = cAC[:(N-1)*p, :]
    cAC_new[(N-1)*p:, :] = np.linalg.matrix_power(A, N)
    cAC = cAC_new
    
    # Update cBC
    cBC_new = zeros((N*p, N*m))
    cBC_new[:(N-1)*p, :(N-1)*m] = cBC[:(N-1)*p, :(N-1)*m]
    
    # Fill in the bottom row for cBC
    for j in range(1, N+1):
        col_idx = (j-1)*m
        if j <= N:
            cBC_new[(N-1)*p:(N-1)*p+n, col_idx:col_idx+m] = np.linalg.matrix_power(A, N-j) @ B
    cBC = cBC_new
    
    # Update cMC
    cMC_new = zeros((N*p, N*q))
    cMC_new[:(N-1)*p, :(N-1)*q] = cMC[:(N-1)*p, :(N-1)*q]
    
    # Fill in the bottom row for cMC
    for j in range(1, N+1):
        col_idx = (j-1)*q
        if j <= N:
            cMC_new[(N-1)*p:(N-1)*p+n, col_idx:col_idx+q] = np.linalg.matrix_power(A, N-j) @ M
    cMC = cMC_new
    
    # Return tuple of matrices
    return cAC, cBC, cMC, cQ, cR, cF, cf


# quadprog as in matlab
def quadprog(H, f, A=None, b=None, Aeq=None, beq=None, lb=None, ub=None, x0=None, options=None):
    """
    Solve quadratic programming problem:
        min 0.5*x'*H*x + f'*x
        subject to:
            A*x <= b
            Aeq*x = beq
            lb <= x <= ub
    """
    import numpy as np
    from scipy.optimize import minimize
    
    # Make sure H is symmetric
    H = (H + H.T) / 2
    
    # Define the objective function - ensure it returns a scalar
    def objective(x):
        x = np.atleast_1d(x)
        return float(0.5 * x.T @ H @ x + f.T @ x)  # Ensure scalar output
    
    # Define the gradient of the objective function
    def gradient(x):
        x = np.atleast_1d(x)
        return H @ x + f.flatten()  # Ensure proper shape
    
    # Set up the constraints
    constraints = []
    
    # Inequality constraints: A*x <= b
    if A is not None and b is not None:
        for i in range(A.shape[0]):
            def constraint_func(x, i=i):
                x = np.atleast_1d(x)
                return float(b[i] - A[i, :] @ x)  # Ensure scalar output
                
            def constraint_grad(x, i=i):
                return -A[i, :]
                
            constraints.append({
                'type': 'ineq', 
                'fun': constraint_func, 
                'jac': constraint_grad
            })
    
    # Equality constraints: Aeq*x = beq
    if Aeq is not None and beq is not None:
        for i in range(Aeq.shape[0]):
            def eq_constraint_func(x, i=i):
                x = np.atleast_1d(x)
                return float(Aeq[i, :] @ x - beq[i])  # Ensure scalar output
                
            def eq_constraint_grad(x, i=i):
                return Aeq[i, :]
                
            constraints.append({
                'type': 'eq', 
                'fun': eq_constraint_func, 
                'jac': eq_constraint_grad
            })
    
    # Set up bounds
    bounds = None
    if lb is not None or ub is not None:
        n = H.shape[0]
        if lb is None:
            lb = np.full(n, -np.inf)
        if ub is None:
            ub = np.full(n, np.inf)
        bounds = list(zip(lb, ub))
    
    # Set default initial point if not provided
    if x0 is None:
        x0 = np.zeros(H.shape[0])
    
    # Parse options
    scipy_options = {}
    if options is not None:
        if 'Display' in options and options['Display'] == 'off':
            scipy_options['disp'] = False
    else:
        scipy_options['disp'] = False  # Default to no display
    
    # Call scipy.optimize.minimize
    result = minimize(
        objective, 
        x0, 
        method='SLSQP',
        jac=gradient,
        constraints=constraints,
        bounds=bounds,
        options=scipy_options
    )
    
    return result.x