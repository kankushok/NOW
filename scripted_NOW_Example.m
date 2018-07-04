% Scripted example of Numerical Optimization of gradient Waveforms (NOW)
clear

% Change the parameters below to your liking. Those not
% specified are default-initialized as follows:
% Max gradient = 80 milliTesla/m
% Max slew rate = 100 milliTesla/m/milliSecond = 100 T/m/s
% Pulse-time = 50 milliSecond
% Eta (heat dissipation parameter) = 1
% Discretization points = 50
% Target tensor = eye(3)
% Initialguess = 'random'
% enforceSymmetry = false;
% redoIfFailed = true;
% name = 'NOW'
% 
% Written by Jens Sj�lund and Filip Szczepankiewicz


%%  PREP
% First, set up the optimization problem. Do this first to create a
% structure where fields are pre-specified. Note that some fields are
% read-only and that new fields cannot be created.
problem = optimizationProblem;

% Define the hardware specifications of the gradient system
problem.gMax =  70; % Maximal gradient amplitude, in [mT/m]
problem.sMax = 100; % Maximal gradient slew, in [T/(sm)]

% Request encoding and pause times based on sequence timing in [ms]
problem.durationFirstPartRequested    = 35;
problem.durationSecondPartRequested   = 28.5;
problem.durationZeroGradientRequested = 10;
problem.ecc_flag = 1;
problem.Ly = [20, 50, 150]'; % eddy current time constants on y axis in ms
problem.Lx = 150;
problem.by = [1, 1, 1/10]'; % tolerance on eddy current amplitude
problem.bx = [1/100]; 
% Define the b-tensor shape in arbitrary units. This example uses an
% isotropic b-tensor that results in spherical tensor encoding (STE).
problem.targetTensor = eye(3);

% Define the number of sample points in time. More points take longer to
% optimize but provide a smoother waveform that can have steeper slopes. 
% The basic code supports N = 50, 100, and 200. However, other values
% can be calculated using the createConstraintGradientFunction in the
% private folder.
problem.N = 77;

% Set the balance between energy consumption and efficacy
problem.eta = 0.5; %In interval (0,1]

% Make a new optimizationProblem object using the updated specifications.
% This explicit call is necessary to update all private variables.
problem = optimizationProblem(problem); 


%% PRINT REQUESTED AND TRUE TIMES
% Note that due to the coarse raster, the requested and actual times may
% differ slightly.
clc
fprintf(1, '------------ Requested timing parameters: ------------ \n');
fprintf(1, 'DurPre = %5.3f  DurPost = %5.3f  DurPi = %5.3f  [ms]\n\n', problem.durationFirstPartRequested, problem.durationSecondPartRequested, problem.durationZeroGradientRequested);
fprintf(1, '------------   Actual timing parameters:  ------------ \n');
fprintf(1, 'DurPre = %5.3f  DurPost = %5.3f  DurPi = %5.3f  [ms]\n\n', problem.durationFirstPartActual, problem.durationSecondPartActual, problem.durationZeroGradientActual);


%% RUN OPTIMIZATION
[result, problem] = NOW_RUN(problem);


%% PLOT RESULT
plot(0:problem.dt:problem.totalTimeActual,result.g)
xlabel('Time [ms]')
ylabel('Gradient amplitude [mT/m]')
measurementTensor = result.B


