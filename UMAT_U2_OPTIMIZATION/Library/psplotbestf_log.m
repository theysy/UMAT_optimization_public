function stop = psplotbestf(optimvalues,flag)
%PSPLOTBESTF PlotFcn to plot best function value.
%   STOP = PSPLOTBESTF(OPTIMVALUES,FLAG) where OPTIMVALUES is a structure
%   with the following fields:
% 
%   PATTERNSEARCH:
%              x: current point X
%      iteration: iteration count
%           fval: function value
%       meshsize: current mesh size
%      funccount: number of function evaluations
%         method: method used in last iteration
%         TolFun: tolerance on function value in last iteration
%           TolX: tolerance on X value in last iteration
%
%   FLAG: Current state in which PlotFcn is called. Possible values are:
%           init: initialization state
%           iter: iteration state
%           done: final state
%
%   STOP: A boolean to stop the algorithm.
%
%   See also PATTERNSEARCH, GA, OPTIMOPTIONS.


%   Copyright 2003-2018 The MathWorks, Inc.

stop = false;
switch flag
    case 'init'
        xlabel('Iteration','interp','none'); 
        ylabel('Function value','interp','none')
        if isscalar(optimvalues.fval)
            title(sprintf('Best Function Value: %g',optimvalues.fval),'interp','none');
        else
            % if we have a multi-objective call to this function, we
            % display a useful title message and disable functionality
            titleStr = getString(message('globaloptim:psplotcommon:PlotFcnUnavailableMultiObj', 'psplotbestf'));
            title(titleStr,'interp','none');
            return;
        end
        plotBest = plot(optimvalues.iteration,optimvalues.fval, '.b');
        set(plotBest,'Tag','psplotbestf');
    case 'iter'
        % no updates for multi-objective
        if ~isscalar(optimvalues.fval)
            return;
        end
        
        plotBest = findobj(get(gca,'Children'),'Tag','psplotbestf');
        newX = [get(plotBest,'Xdata') optimvalues.iteration];
        newY = [get(plotBest,'Ydata') optimvalues.fval];
        set(plotBest,'Xdata',newX, 'Ydata',newY);
        set(get(gca,'Title'),'String',sprintf('Best Function Value: %g',optimvalues.fval));
        set(gca, 'YScale', 'log');
end
