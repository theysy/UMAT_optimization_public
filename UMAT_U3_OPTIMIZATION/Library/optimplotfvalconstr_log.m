function stop = optimplotfvalconstr(~,optimValues,state)
% optimplotfvalconstr Plot value of the objective function vs iteration.
% Infeasible points are marked red. If no fval is available, this function
% will plot constraint violation (infeasibility).
%

%   Copyright 2019 The MathWorks, Inc.

persistent plotBest plotBestInfeas legendHndl legendStr legendHndlInfeas ... 
    legendStrInfeas tolCon feasLegend infeasLegend

stop = false;

if strcmpi(state,'init')
    plotBest = []; 
    plotBestInfeas = [];     
    legendHndl = []; legendStr = {};
    legendHndlInfeas = []; legendStrInfeas = {};
    feasLegend = []; infeasLegend = [];
end

if optimValues.funccount == 0 || (isempty(optimValues.fval) && isempty(optimValues.constrviolation))
    % no function evals or none of the trials are successfully evaluated; no plots.
    return;
end

haveConstr = isfield(optimValues, 'constrviolation');

if isempty(optimValues.fval)    
    % Not worth an error message but protect against misuse
    assert(haveConstr)
    % feasibility problem
    best = optimValues.constrviolation;
else
    best = optimValues.fval;
    % Not worth an error message but protect against misuse
    assert(isscalar(best))
end

if isempty(plotBest) && isempty(plotBestInfeas)
    xlabel('Iteration','interp','none');
    
    if isempty(optimValues.fval)
        ylabel('Infeasibility','interp','none');
        title('Current infeasibility: ','interp','none')
    else
        ylabel('Objective Function','interp','none');
        title('Best function value: ','interp','none')
    end
    hold on; grid on;
    
    % Use a tolerance to mark infeasible points.
    fig = findobj(0,'Type','figure','name','Optimization PlotFcns');
    if ~isempty(fig)
        options = get(fig,'UserData');
        if ~isempty(options)
            tolCon = options.ConstraintTolerance;
        else
            tolCon = 1e-6;
        end
    else
        tolCon = 1e-6;
    end
    
end

                        
if ~haveConstr || optimValues.constrviolation <= tolCon
    % Plot feasible points
    if isempty(plotBest)
        plotBest = plot(optimValues.iteration,best,'b.');
        set(plotBest,'Tag','plotbestf');
        legendHndl(end+1) = plotBest;
        legendStr{end+1} = 'Best fval';
    else
        newX = [get(plotBest,'Xdata') optimValues.iteration];
        newY = [get(plotBest,'Ydata') best];
        set(plotBest,'Xdata',newX, 'Ydata',newY);
    end
else % Plot infeasible points
    if isempty(plotBestInfeas)
        plotBestInfeas = plot(optimValues.iteration,best,'r.');
        set(plotBestInfeas,'Tag','plotbestfinfeas');
        legendHndlInfeas(end+1) = plotBestInfeas;
        legendStrInfeas{end+1} = 'Best fval (Infeas)';
    else
        newX = [get(plotBestInfeas,'Xdata') optimValues.iteration];
        newY = [get(plotBestInfeas,'Ydata') best];
        set(plotBestInfeas,'Xdata',newX, 'Ydata',newY);
    end
end
if isempty(optimValues.fval)
    set(get(gca,'Title'),'String', sprintf('Current infeasibility: %g', best),'interp','none')
else
    set(get(gca,'Title'),'String', sprintf('Best function value: %g', best),'interp','none')
end



if isempty(feasLegend) && (~haveConstr || optimValues.constrviolation <= tolCon)
    % Feasible phase legends
    feasLegend = legend([legendHndlInfeas legendHndl], [legendStrInfeas legendStr]);
elseif isempty(infeasLegend)
    % Infeasible phase legend
    infeasLegend = legend(legendHndlInfeas, legendStrInfeas,'FontSize',8);
end
set(gca, 'YScale', 'log');
if strcmp(state, 'done')
    hold off;
end
end

