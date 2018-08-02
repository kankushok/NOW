function [f ,g] = objFun(x)
    f = -x(end); %fmincon minimizes the objective, so a minus sign is used to maximize instead 
    g = zeros(size(x));
    g(end) = -1;  
    
    currentTime = toc;
    if( currentTime > 1000)
        error('optimization timed out');
    end
