clc; clear all; close all;
addpath('Sources')

% define the simple cruise control model
mass = 1e3; % [kg] mass of the car
damp = 10; % [Ns/m] damping coefficient
ms2kmh = 3.6; % [m/s] to [km/h]
dstrb = mass * 9.81 * sin(deg2rad(3)); % [N] disturbance force (slope of 5 degrees)
start_dstrb = 120; % [s] start of the disturbance

ref = 50; % [km/h] reference speed

u_max = 400; % [N] maximum control input

dt = .1; % [s] time step simulation
T = 500; % [s] simulation time

Ts = 1.0; % [s] sampling time 

td = linspace(0, T, T/Ts);

% continuous time state space model
Ac = -damp/mass;
Bc = 1/mass;
Cc = ms2kmh;

M = -1/mass; % disturbance matrix

% Exact Discretization using Matlab Toolbox 
sys_c = ss(Ac,Bc,Cc,0);
sys_d = c2d(sys_c, Ts);

% System
A = sys_d.A; B = sys_d.B; C = sys_d.C;

% define a function to perform a single simulation step, returns the new state and output
function [x1, y] = sim_step(x, u, d, A, B, C, M)
    x1 = A*x + B*u + M*d;
    y = C*x;
end

% define vector of states, outputs, inputs and disturbances
xs = zeros(length(td), 1); % states
ys = zeros(length(td), 1); % outputs
us = zeros(length(td), 1); % inputs
ds = dstrb * heaviside(td - start_dstrb); % disturbances

% simulate the system with a for loop, fixed control
for i = 2:length(td)
    % compute control input
    u = 900;
    
    % simulate the system
    [xs(i), ys(i)] = sim_step(xs(i-1), u, ds(i-1), A, B, C, M);
    us(i) = u;
end

% figure('Position', [0 0 2500 1500]); subplot(2,1,1); plot(td, xs); grid on; subplot(2,1,2); plot(td, us); grid on;

% return

%% MPC
N = 50; % prediction horizon
Q = 500; % state cost
R = 1; % input cost
S = 0; % terminal cost

F = [1.; -1.]; % input constraints
f = [u_max u_max]; % input constraints

% create the condensed matrices
[cA, cB, cM, cQ, cR, cF, cf] = BuildCondensedMPCmatrices(A, B, C, M, Q, R, S, F, f, N);

% define reference output and control input
cYr = ref * ones(length(td), 1); % reference output
cUr = zeros(N, 1); % reference input

% define vector of states, outputs, inputs and disturbances
xs = zeros(length(td), 1); % states
ys = zeros(length(td), 1); % outputs
us = zeros(length(td), 1); % inputs

% ds = dstrb * heaviside(td - start_dstrb); % disturbances
% ds = dstrb * ones(length(td)); % disturbances
ds = - dstrb * heaviside(td - start_dstrb) + .5*dstrb * heaviside(td - 2*start_dstrb); % 2 disturbances

% simulate the system with a for loop, MPC control
for i = 2:length(td)-N
    % mpc
    x = xs(i-1); d = ds(i-1);

    % cDk = ds(i-1)*ones(N, 1); % fixed disturbance for the prediction horizon
    cDk = ds(i-1:i+N-2)'; % varying disturbance for the prediction horizon

    Hqp = cB'*cQ*cB + cR;
    Hqp = (Hqp + Hqp')/2; % make sure it is symmetric

    % fqp = (cA*x + cM*cDk - cYr(i-1))' * cQ * cB; % - cUr' * cR
    fqp = cB'*cQ*(cA*x + cM*cDk - cYr(i-1)); % - cUr' * cR

    Aqp = cF;
    bqp = cf;

    % compute control input
    cU = quadprog(Hqp, fqp, Aqp, bqp, [],[],[],[],[], optimset('Display','off'));
    u = cU(1);
    
    % simulate the system
    [xs(i), ys(i)] = sim_step(x, u, d, A, B, C, M);
    us(i) = u;
end

% keep onlt the first length(td)-N elements
td = td(1:end-N);
xs = xs(1:end-N);
ys = ys(1:end-N);
us = us(1:end-N);
ds = ds(1:end-N);

figure('Position', [0 0 2500 1500]); subplot(3,1,1); 
plot(td, ys); hold on; plot(td, ref*ones(size(td)), 'r--'); hold off; grid on; title('Output');
subplot(3,1,2); plot(td, us); grid on; title('Input');
subplot(3,1,3); plot(td, ds); grid on; title('Disturbance');