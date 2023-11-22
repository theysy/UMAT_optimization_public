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
%   (1) anisotropic for the forward-reverse loading
%   (2) dislcoation-based hardening model (RGBV)
clear;
close all;
clc;
addpath 'Data'
addpath 'Library'
addpath 'UMAT_props'
tic
%%   Compilation of UMAT subroutine using MEX.
%    mex MML_UMAT.f
%%  Identification opt_par: Read PROPS.csv
props=csvread('props_TR1180_HAH20.csv');
optparam.props=props;
%%   SET plasticity model parameters
% initial guess for Hardening models
% [0] Isotropic hardening coefficients
if props(1)==0
%     X=[2.5 38 85000 2.46e-10 180 30e-6 2.8 0.8 0 1e12];
%     lb=[2.24 0 85000 2.46e-10 15 1e-6 0.0 0.8 0 1e10];
%     ub=[3.04 100 85000 2.46e-10 580 60e-6 20.0 0.8 0 1e14];
    X=props(9:18); % K, K1, K2, K3, K4, K5
    lb=X;
    ub=X;
% [1] Chaboche
elseif props(1)==1
    X=props(21:29); % K, K1, K2, K3, K4, K5
    lb=props(21:29)*0.5;
    ub=props(21:29)*2.0;
% [2] Yoshida-Uemori
elseif floor(props(1))==2
    X=[1 1 1 1 1 1 1];
%     lb=props(21:21+xdim-1)*0.1;
%     ub=props(21:21+xdim-1)*2.0;
    lb=X;
    ub=X;
% [3] HAH11
elseif props(1)==3
    X=[200, 120, 0.35, 0.9, 15]'; % K1, K2, K3, K4, K5
     lb=[50 50 0.1 0.8 1];
     ub=[200 200 1 1 30];
%     lb=X;
%     ub=X;
% [4] HAH14
elseif props(1)==4
    X=[100, 150, 0.55, 0.85, 5]'; % K1, K2, K3, K4, K5
%     lb=[50 50 0.3 0.8 1.0];
%     ub=[200 200 1.0 1.0 30];
    lb=X;
    ub=X;
% [5] HAH20
elseif floor(props(1))==5
    X=[300, 300, 0.5, 0.93, 15]'; % K1, K2, K3, K4, K5
    lb=[50, 50, 0.3, 0.9, 1];
    ub=[450 450 1.0 0.95 60];
%     X=props(23:27); % K1, K2, K3, K4, K5
%     lb=X;
%     ub=X;
end
xdim=max(size(X));
%% SET optimization parameters
if props(2)==4
    optparam.ntens=6;    
else
    optparam.ntens=3;
end
optparam.nstatv=70;
% optparam.dt=2.5e-3;
optparam.dt=2.5e-3;
optparam.xdim=xdim;
optparam.lb=lb;
optparam.ub=ub;
optparam.weight(1)= 1.0;
optparam.weight(2)= 5.0;
optparam.filemat{1} = 'TR1180_EXP_TCT5P';
optparam.mode=11; % Tension-Compression (cyclic)
optparam.optalg=1; % 1: Nelder-Mead | 2: Genetic | 3: Pattensearch | 4: Globalsearch
if optparam.optalg==1 || optparam.optalg==4
    optparam.plotmode=1; % 0: No plotting | 1: Plotting
else
    optparam.plotmode=0; % 0: No plotting | 1: Plotting
end
optparam.error=2; % 1: sum_sqaure | 2: sum_abs | 3: root_mean_square | 4: mean_abs
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
if props(1)==0 % Isotropic hardening
%     props(10:10+xdim-1)=Poptim;
    props(9:9+xdim-1)=Poptim;
elseif props(1)==1 % Chaboche
    props(21:21+xdim-1)=Poptim;
elseif floor(props(1))==2 % Yoshida-Uemori
    props(21:21+xdim-1)=Poptim;
elseif props(1)==3 % HAH11
    props(23:23+xdim-1)=Poptim;
elseif props(1)==4 % HAH14
    props(23:23+xdim-1)=Poptim;
elseif floor(props(1))==5 % HAH20
    props(23:23+xdim-1)=Poptim;
end
csvwrite('PROPS_FINAL.csv',props);
disp('-----------------------------------------------');
disp('Optimized coefficients:');
disp(Poptim);
disp('-----------------------------------------------');
%%
%     filename=strcat('curve',num2str(k),'.csv');
    %csvwrite('curve.csv', [strainUMAT(:,1),stressUMAT(:,1)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc



