%--------------------------------------------------------------------------
% Jin-Hwan Kim, MML, GIFT, POSTECH
% Nov 7 2017 
% E-mail: jinkim@postech.ac.kr
% Ver.: 0.5

% Revised by Shinyeong Lee, MML, GIFT, POSTECH
% Aug 10 2019 
% E-mail: sylee413@postech.ac.kr
% Ver.: 1.0

% Revised by Seong-Yong Yoon, MML, GIFT, POSTECH
% Nov 26 2020 
% E-mail: theysy@postech.ac.kr
% Ver.: 2.0
% Main feature: Simulation is run based on UMAT subroutine.
%--------------------------------------------------------------------------
% This matlab script is written for indentification of
%   anisotropic yield function
clear;
close all;
clc;
addpath 'Data'
addpath 'Library'
addpath 'UMAT_props'
%%   Compilation of UMAT subroutine using MEX.
%    mex MML_UMAT.f
%%  Identification opt_par: Read PROPS.csv
props=csvread('props_AA2090_yld2k.csv');
optparam.props=props;
%%   SET plasticity model parameters
% initial guess for Anisotropic yield function
% Hill 1948
if props(2)==2
    X=[1 1 1 1 1 1];
    lb=[0 0 0 0 0 0];
    ub=[2 2 2 2 2 2];
% Yld2000-2d
elseif props(2)==3
    X=ones([1,8]);
    lb=ones(size(X)).*(0);
    ub=ones(size(X)).*(2);
%     X=props(33:40);
%     lb=X;
%     ub=X;
% Yld2004-18p
elseif props(2)==4
%     X=ones([1,18]);
%     lb=ones(size(X)).*(-0.1);
%     ub=ones(size(X)).*(2);
    X=props(33:50);
    lb=X;
    ub=X;
else
    X=[1 1 1 1 1 1 1 1];
    lb=X;
    ub=X;
end
xdim=max(size(X));
%% SET optimization parameters
if props(2)==4 || props(2)==1
    optparam.ntens=6;    
else
    optparam.ntens=3;
end
optparam.nstatv=70;
% optparam.dt=2.5e-3;
optparam.dt=1e-3;
optparam.xdim=xdim;
optparam.lb=lb;
optparam.ub=ub;
optparam.weight(1)= 1;
optparam.weight(2)= 0.1;
optparam.weight(3)= 0.0;
optparam.filemat{1} = 'AA2090_DATA';
optparam.mode=20; % Anisotropic yield function identification
optparam.optalg=1; % 1: Nelder-Mead | 2: Genetic | 3: Pattensearch | 4: Globalsearch
if optparam.optalg==1 || optparam.optalg==4
    optparam.plotmode=1; % 0: No plotting | 1: Plotting
else
    optparam.plotmode=0; % 0: No plotting | 1: Plotting
end
optparam.error=1; % 1: sum_sqaure | 2: sum_abs | 3: root_mean_square | 4: mean_abs
%% Reference stress state
%   ntens=3: s=[s11 s22 s12]
%   ntens=6: s=[s11 s22 s33 s12 s13 s23]
    s=zeros([1,optparam.ntens]);
    s(1)=1; % Uniaxial tension
%     s(3)=1; % Simple shear
optparam.s=s;
%% optimization algorithm
% Normalize input and boundaries 
nX=norm_coeff(X,lb,ub,1);
nlb=norm_coeff(lb,lb,ub,1);
nub=norm_coeff(ub,lb,ub,1);
% You can choose optimization algorithm
if optparam.optalg==1
    % 1. Nelder-Mead Simplex algorithm
    options=optimset('Display','iter','PlotFcns',@optimplotfval);
    options.TolFun=1e-4;
    options.TolX=1e-3;
    options.MaxIter=400;
    options.MaxFunEvals=1e+04;
    [Poptim,fval,exitflag,output] = fminsearchbnd(@(nX) UMAT_OPT(nX,optparam),nX,nlb,nub,options);
elseif optparam.optalg==2
    % 2. Genetic algorithm
    options=optimoptions('ga','Display','iter','PlotFcn',{@gaplotbestf_log});
    options.UseParallel=true; 
    options.MaxGenerations=20*xdim;
    options.MaxStallGenerations=10;
%     options.PopulationSize=500;
%     options.EliteCount=ceil(0.05*options.PopulationSize);
%     options.MutationFcn=@mutationadaptfeasible;
    options.FunctionTolerance=1e-5;
%     options.CrossoverFcn={@crossoverintermediate, 0.8};
    [Poptim,fval,exitflag,output,population,scores] = ga(@(nX) UMAT_OPT(nX,optparam),xdim,[],[],[],[],nlb,nub,[],options);
elseif optparam.optalg==3
    % 3. Pattensearch
    options = optimoptions('patternsearch','Display','iter','PlotFcn',{@psplotbestf_log, @psplotmeshsize});
%     options.MaxIterations=1000;
    options.InitialMeshSize=1;
%     options.MeshTolerance=1e-5;
    options.UseParallel=true;
    options.StepTolerance=1e-4;
    options.PollMethod='MADSPositiveBasis2N';
    [Poptim,fval] = patternsearch(@(nX)UMAT_OPT(nX,optparam),nX,[],[],[],[],nlb,nub,[],options);
elseif optparam.optalg==4
    % 5.1 Multi-start
%     ms = MultiStart('Display', 'iter','StartPointsToRun','bounds', 'PlotFcn', {@gsplotbestf; @gsplotfunccount});
%     ms.FunctionTolerance=2e-4;
%     ms.XTolerance=1e-3;
%     ms.UseParallel=false;
%     gs= GlobalSearch(ms);
    % 5.2 Globalsearch
    gs = GlobalSearch('StartPointsToRun','bounds','Display','iter','PlotFcn',{@gsplotbestf; @gsplotfunccount});
    gs.FunctionTolerance=1e-2;
    gs.XTolerance=1e-2;
    options = optimoptions(@fmincon,'Algorithm','interior-point','CheckGradients',true);
    options.FunctionTolerance=1e-2;
    options.StepTolerance=1e-4;
    problem = createOptimProblem('fmincon','x0',nX,'objective',@(nX)UMAT_OPT(nX,optparam),'lb',nlb,'ub',nub,'options',options);
    [Poptim,fval] = run(gs,problem);
end
%%
niter=size(optparam.filemat,2);
optparam.plotmode=1;
gap=UMAT_OPT(Poptim,optparam);
Poptim=norm_coeff(Poptim,lb,ub,2);
% Hill 1948
if props(2)==2
    props(32:32+xdim-1)=Poptim;
% Yld2000-2d
elseif props(2)==3
    props(33:33+xdim-1)=Poptim;
% Yld2004-18p
elseif props(2)==4
    props(33:33+xdim-1)=Poptim;
else
    props(33:33+xdim-1)=Poptim;
end
csvwrite('PROPS_FINAL.csv',props);
disp('-----------------------------------------------');
disp('Optimized coefficients:');
disp(Poptim');
disp('-----------------------------------------------');
%%
%     filename=strcat('curve',num2str(k),'.csv');
    %csvwrite('curve.csv', [strainUMAT(:,1),stressUMAT(:,1)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



