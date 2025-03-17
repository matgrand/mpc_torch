clc
close all
% bdclose all
clear all

addpath('Sources')

%--------------------------------------------------------------------%
%% Step 0: Load Current Execise Data 
%         (needed to run the simulation)
%--------------------------------------------------------------------%
currentExercise=1;

% Possible simplifications (activate by setting the variable to 1):
removeContraints=0;

exerciseData=LoadMPCExercise(currentExercise,removeContraints);

% Relavant data for the controller design
exerciseData.u_max
exerciseData.u_min

%-----------------------------%
%%   MPC controller design   %%
%-----------------------------%
%% Step 1: upload a model of the system to be controlled
    % Possible simplifications (activate by setting the variable to 1):
    exactModel=1;
    SystemModel=LoadSystemModel(currentExercise,exactModel);
    
%% Step 2: use to model to create the condensed MPC matrices
% You are provided with the command
%    condensedMatrices=BuildCondendsedMPCmatrices(SystemModel,R,Q,F,f,PH); % Warning: this command will produce an error, since R,Q,F,f and PH are currently undefined.

%------------%
%% Option 2: Using a matlab for cycle:
%------------%
% Initalization
x_k=zeros(size(exerciseData.System.A,1),1);

u_memory=0;
x_memory=x_k;
y_memory=exerciseData.System.C*x_k;

% Basic Closed-Loop system simulation
for k=exerciseData.Ts:exerciseData.Ts:exerciseData.SimTime

    % Compute the the control action to be used using MPC
    u_k=1000; % in this example always u=1     

    % Simulate Plant Response
    x_kplus1=exerciseData.System.A*x_k+exerciseData.System.B*u_k; 
    y=exerciseData.System.C*x_k+exerciseData.System.D*u_k;
    
    % Save for plots
    u_memory(:,end+1)=u_k;
    x_memory(:,end+1)=x_k;
    y_memory(:,end+1)=y;

    % Update and repeat    
    x_k=x_kplus1;
end

% Basic Plot
figure
subplot(2,1,1)
    plot(0:exerciseData.Ts:exerciseData.SimTime,y_memory)
    xlabel('Time [s]')
    if currentExercise==1
        ylabel('Speed [km/h]')
    elseif currentExercise==3
        ylabel('Position [m]')
    end
    title('System Output')
subplot(2,1,2)
    plot(0:exerciseData.Ts:exerciseData.SimTime,u_memory)
    xlabel('Time [s]')
    ylabel('Force [m]')
    title('Control Action')
%------------%    