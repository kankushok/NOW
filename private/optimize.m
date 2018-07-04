function [result, problem] = optimize(varargin)
% Script that optimizes an MR gradient waveform subject to a number of constraints.
% In particular it maximizes the b-value of a diffusion encoding
% pulse sequence subject to constraints on power, maximum gradient and maximum slew rate.
%
% If you use this in your research, please cite the following paper:
% Jens Sj�lund, Filip Szczepankiewicz, Markus Nilsson, Daniel Topgaard, Carl-Fredrik Westin, Hans Knutsson,
% "Constrained optimization of gradient waveforms for generalized diffusion encoding",
% Journal of Magnetic Resonance, Volume 261, December 2015, Pages 157-168, ISSN 1090-7807,
% http://dx.doi.org/10.1016/j.jmr.2015.10.012.
% (http://www.sciencedirect.com/science/article/pii/S1090780715002451)
%
% If you use asymmetric waveforms with Maxwell compensation, please cite the following abstract (or later paper):
% Filip Szczepankiewicz and Markus Nilsson
% "Maxwell-compensated waveform design for asymmetric diffusion encoding"
% ISMRM 2018, Paris, France
% Download PDF at: https://goo.gl/vVGQq2
%
%
% Written by Jens Sj�lund, jens.sjolund@liu.se / jens.sjolund@elekta.com
% Maxwell compensation by Filip Szczepankiewicz, filip.szczepankiewicz@med.lu.se

%% Initialize parameters
if nargin == 0
    problem = optimizationProblem();
elseif nargin == 1
    problem = optimizationProblem(varargin{1});
else
    error('Invalid number of input arguments')
end


%% Set optimization parameters
options = optimoptions('fmincon','Algorithm','sqp',...
    'DerivativeCheck','off','Display','off',...
    'GradObj','on','GradConstr','on','MaxFunEval',1e5,'MaxIter',5e3);
warning('off', 'optimlib:fmincon:ConvertingToFull'); %Disables warning when SQP converts sparse matrices to full

%% Set up constraints
% granty edit for eddy current nulled sequence...
[A, b, firstDerivativeMatrix, secondDerivativeMatrix] = defineLinearInequalityConstraints(problem.N, problem.gMaxConstraint, problem.sMaxConstraint, problem.useMaxNorm, problem.zeroGradientAtIndex, problem.dt, problem.ecc_flag, problem.Lx, problem.Ly, problem.Lz, problem.bx, problem.by, problem.bz);

[Aeq, beq] = defineLinearEqualityConstraints(problem.N, problem.zeroGradientAtIndex, problem.enforceSymmetry, firstDerivativeMatrix);

% Define nonlinear inequality constraints
nonlconFileName = getNonLinearConstraintsFileName(problem.N, problem.useMaxNorm);
if ~exist(nonlconFileName,'file')
    createConstraintGradientFunction(problem.N,problem.useMaxNorm); %Uses the symbolic toolbox to derive Jacobian ,SLOW!
end

%% Optimize
optimizationSuccess = false;
iter = 1;
while ~optimizationSuccess && iter <= 10
    
    x0 = getInitialGuess(problem, iter);
    dispInfo(problem, iter)
    
    tic
    
	[x,fval,exitflag,output,lambda,grad]  = fmincon(@(x) objFun(x), x0, A,b,Aeq,beq,[],[],@(x) feval(nonlconFileName,x,problem.tolIsotropy, ...
											problem.gMaxConstraint, problem.integralConstraint,problem.targetTensor, problem.tolMaxwell, ...
											problem.s_vec),options);
	
    optimizationTime = toc;
    
    disp(['Optimization took ' num2str(optimizationTime, 3) 's.']);
    
    optimizationSuccess = (exitflag > 0);
    
    iter = iter +1;
end

%% Evaluate and store results
gamma = 42.6e6*2*pi;
q = reshape(x(1:3*problem.N),[problem.N,3]);
g = [zeros(1,3);firstDerivativeMatrix*reshape(q,[problem.N 3]);zeros(1,3)]/problem.dt;
slew = [zeros(1,3);secondDerivativeMatrix*reshape(q,[problem.N 3]);zeros(1,3)]/(problem.dt)^2;
q = gamma*1e-6*q; %SI-units
q0 = reshape(x0(1:3*problem.N),[problem.N,3]);
q0 = gamma*1e-6*q0;
B = problem.dt*1e-3*(q'*q);
b = trace(B)*1e-6;%s/mm^2
C = b/(1e-6*gamma^2*(problem.gMaxConstraint*1e-3/problem.dt)^2*(problem.totalTimeActual*1e-3)^3);
kappa = 4*C;

etaOpt = problem.dt*max(diag(g'*g))/(problem.gMaxConstraint^2*problem.totalTimeActual);

result.q = q;
result.g = g;
result.slew = slew;
result.b = b;
result.B = B;
result.q0 = q0;
result.kappa = kappa;
result.etaOpt = etaOpt;
result.optimizerOutput = output;
result.optimizerProblem = problem;


% Format output in GWF format
rf = ones(size(g,1), 1);

if ~isempty(problem.zeroGradientAtIndex)
    zi = problem.zeroGradientAtIndex+1; 
    
    rf(zi) = 0;
    
    rf((max(zi)+1):end) = -rf((max(zi)+1):end);
end

result.rf  = rf;                                % Spin direction
result.gwf = g/1000 .* repmat(rf, 1, 3);        % T/m
result.dt  = result.optimizerProblem.dt/1000;   % s

end

function [A, b, firstDerivativeMatrix, secondDerivativeMatrix] = defineLinearInequalityConstraints(N, gMaxConstraint, sMaxConstraint, useMaxNorm,zeroGradientAtIndex, dt, ecc_flag, lx, ly, lz, bx, by, bz)
firstDerivativeMatrix = -diag(ones(N,1))+diag(ones(N-1,1),1); % Center difference, shifted forward by half a step. Ghost points implemented as zero rows.
firstDerivativeMatrix = firstDerivativeMatrix(1:end-1,:);
firstDerivativeMatrix = sparse(firstDerivativeMatrix); %SQP doesn't take advantage of this

if useMaxNorm == true %This is if we want to use max-norm on the gradients
    A1 = kron(eye(3),firstDerivativeMatrix);
    A1 = [A1 zeros(size(A1,1),1)]; %Add column of zeros for s
    b1 = gMaxConstraint*ones(size(A1,1),1);
else
    A1 = [];
    b1 = [];
end

%Constraint on change in gradients (slew rate)
secondDerivativeMatrix = diag(ones(N-1,1),-1)-2*diag(ones(N,1))+diag(ones(N-1,1),1);
secondDerivativeMatrix = sparse(secondDerivativeMatrix);%SQP doesnt take advantage of this
A2 = kron(eye(3),secondDerivativeMatrix);
A2 = [A2 zeros(size(A2,1),1)]; %Add column of zeros for s
b2 = sMaxConstraint*ones(size(A2,1),1);

%Granty edit to add eddy current nulled constraint
if(ecc_flag)
    A3 = [];
    b3 = [];
    % eddy current constraints at time constants 
%     bx = 1/10*ones(size(lx));
%     by = 1/10*ones(size(ly));
%     bz = 1/10*ones(size(lz));
    for i = 1:3
        ec = [];
        At = [];
        if(i==1&&~isempty(lx))       
            ec = exp(-bsxfun(@times,[0:N-1]*dt,1./lx));
            b3 = [b3 ; bx ;bx ;bx];
        end
        if(i==2&&~isempty(ly))
            ec = exp(-bsxfun(@times,[0:N-1]*dt,1./ly));
            b3 = [b3; by; by; by];
        end
        if(i==3&&~isempty(lz))
            ec = exp(-bsxfun(@times,[0:N-1]*dt,1./lz));
            b3 = [b3; bz; bz; bz];
        end

        if(~isempty(ec))
            H = fliplr(ec); % really?? I should just update the ec equation
            P = zeros(N,1); % we want the real slew rate as played on scanner
            P(1:zeroGradientAtIndex(1)) = 1;
            P(zeroGradientAtIndex(end):end) = -1;
            P = diag(P);
            At = H*P*secondDerivativeMatrix;
            tmp = zeros(1,3);
            tmp(i) = 1;
            At = kron(diag(tmp),At);
            At = [At zeros(size(At,1),1)]; %Add column of zeros for s
            
        end
        A3 = [A3 ; At];
    end
else
    A3 = [];
    b3 = [];
end

A = [A1;-A1;A2;-A2;A3;-A3]; %abs(Ax)<=b is equivalent to -Ax<=b && Ax<=b
b = [b1;b1;b2;b2;b3;b3];
end


function [Aeq, beq] = defineLinearEqualityConstraints(N, zeroGradientAtIndex, enforceSymmetry, firstDerivativeMatrix)
Aeq = zeros(2+length(zeroGradientAtIndex),N);
% Require start and end in q-space origin (echo condition)
Aeq(1,1)=1;
Aeq(2,N)=1;
% Require zero gradient at the specified indices
Aeq(2+(1:length(zeroGradientAtIndex)),:) = firstDerivativeMatrix(zeroGradientAtIndex,:);

% Enforce symmetry about zero gradient interval
if enforceSymmetry == true
    if isempty(zeroGradientAtIndex)
        indicesBefore = 1:floor(N/2);
        indicesAfter = floor(N/2) + (1:ceil(N/2));
    else
        indicesBefore = 1:zeroGradientAtIndex(1);
        indicesAfter = (zeroGradientAtIndex(end)+1):N;
    end
    
    assert(length(indicesBefore) == length(indicesAfter),'Cannot enforce symmetry since the number of time samples before and after zero gradient interval is not equal.')
    
    Nactivated = length(indicesBefore);
    Aeq = [Aeq;fliplr(eye(Nactivated)), zeros(Nactivated,N-2*Nactivated),-eye(Nactivated)]; %q(1) = q(end) and so on.
end

Aeq = kron(eye(3),Aeq);
Aeq = sparse([Aeq zeros(size(Aeq,1),1)]); %Add column of zeros for s
beq = zeros(size(Aeq,1),1);
end

function x0 = getInitialGuess(optimization, iter)
if strcmp(optimization.initialGuess,'user-provided')
    if iter > 1
        warning('Optimization with user-provided initial guess seems to have failed. Trying with random initialization instead.')
        x0 = randn(3*optimization.N+1,1);
    else
        x0 = optimization.x0; % Return user-specified initial guess.
    end
else
    if ~strcmp(optimization.initialGuess,'random') && iter == 1
        warning('Unrecognized initial guess string. Defaulting to random initialization.')
    end
    x0 = randn(3*optimization.N+1,1);
end
end

function dispInfo(problem, iter)
if iter == 1
    disp(['Optimizing ' problem.name])
else
    disp(['Optimizing ' problem.name ', attempt ' num2str(iter)])
end
end

