function [error] = cal_error(A,B,iopt)
% dataA: Simulation data
% dataB: Experimental data
% Calculate statistical error
if iopt==1 % sum_square
    error=sum(abs(A./B-1).^2);
elseif iopt==2 % sum_abs
    error=sum(abs(A./B-1).^1);
elseif iopt==3 % root_mean_square
    error=sqrt(mean((A./B-1).^2));
elseif iopt==4 % mean_abs
    error=mean(abs(A./B-1));
end
if isnan(error)
    error=100;
end
end

