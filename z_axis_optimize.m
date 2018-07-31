function [ result] = z_axis_optimize( result,problem )
%z_axis_optimize Recomputes the third gradient axis of a two axis ecc sequence using a convex relaxation of the problem 
%   avoids issues with degenerate waveforms on the z-axis.
N = problem.N;
dt = problem.dt;
Gmax = problem.gMax;
SR = problem.sMax;
qx = result.q(:,1);
qy = result.q(:,2);
firstDerivativeMatrix = -diag(ones(N,1))+diag(ones(N-1,1),1); % Center difference, shifted forward by half a step. Ghost points implemented as zero rows.
firstDerivativeMatrix = firstDerivativeMatrix(1:end-1,:);
firstDerivativeMatrix = sparse(firstDerivativeMatrix); %SQP doesn't take advantage of this
A1 = firstDerivativeMatrix/dt/42.58/2/pi;
%Constraint on change in gradients (slew rate)
secondDerivativeMatrix = diag(ones(N-1,1),-1)-2*diag(ones(N,1))+diag(ones(N-1,1),1);
secondDerivativeMatrix = sparse(secondDerivativeMatrix);%SQP doesnt take advantage of this
A2 = secondDerivativeMatrix/dt^2/42.58/2/pi;
zGI = result.optimizerProblem.zeroGradientAtIndex;
ec = exp(-[0:N-1]*dt/10);
pol = ones(N,1);
pol(zGI(1):end) = -1;
cvx_begin 
    variables qzi(N,1)
    minimize sum(-qzi) 
    subject to
        % gradient heating limitations
        norm(A1*qzi)*dt<=sqrt(problem.eta*Gmax.^2*N*dt);
        % q-trajectory orthogonality
        norm(transpose(qx)*qzi) <= 200;
        norm(transpose(qy)*qzi) <= 200;
        % gradient is zero during RF pulse
        A1(zGI,:)*qzi == 0;
        % gradient limitation
        A1*qzi <= Gmax;
        A1*qzi >= -Gmax;
        % slew rate limitation
        A2*qzi <= SR;
        A2*qzi >= -SR;
        % q-space trajectory starts and ends at 0
        qzi(end) == 0;
        qzi(1) == 0;
        
cvx_end

qzf = qzi/sqrt((qzi'*qzi)/(qx'*qx));
result.q = [qx qy qzf];
result.g = [0 0 0; A1*result.q; 0 0 0];

end

