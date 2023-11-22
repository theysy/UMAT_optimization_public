function gap = UMAT_OPT(nX,optparam) 
% display the current K parameters
X=norm_coeff(nX,optparam.lb,optparam.ub,2);
format shortEng
disp(X);

%% [0] Define parameters
filemat=optparam.filemat;
weight=optparam.weight;
ntens=optparam.ntens;
error_par=optparam.error;
props=optparam.props;
xdim=optparam.xdim;
if ntens==3
    ndi=2;
else
    ndi=3;
end
dt=optparam.dt;
mode=optparam.mode;
nstatv=optparam.nstatv;
s=optparam.s; % Reference stress state
plotmode=optparam.plotmode;
% Index for reference stress state
for i=1:ntens
    if s(i)~=0
        rindx=i;
    end
end
nsim=size(filemat,2);
enu=props(8);
gap=0;
%% [1] Reverse loading effects optimization: [k1-k5]
if mode==11
%--------------------------------------------------------------------------
% [1.0] DEFINE PROPS
if props(1)==0
    props(9:9+xdim-1)=X;
elseif props(1)==1
    props(21:21+xdim-1)=X;
elseif floor(props(1))==2
    props(21:21+xdim-1)=X;
elseif props(1)==3
    props(23:23+xdim-1)=X;
elseif props(1)==4
    props(22:22+xdim-1)=X;
elseif floor(props(1))==5
    props(23:23+xdim-1)=X;
end

%--------------------------------------------------------------------------
for k=1:nsim % Number of simulations
%--------------------------------------------------------------------------
% [1.1] SET UP BOUNDARY CONDITION: DSTRAN
    load(filemat{k});
    nstep=max(size(bc));
    time=nstep/dt;
    bc0=exp(bc)-1;
    stran=zeros([time,ntens]);
    for i=1:nstep
        if i==1
            stran(1:1/dt,rindx)=linspace(0,bc0(i),1/dt);
        else
            stran((i-1)/dt+1:i/dt,rindx)=linspace(bc0(i-1),bc0(i),1/dt);
        end
    end
    stran=log(1+stran);
    if rindx<=ndi
        for m=1:ndi
            if m==rindx
                stran(:,m)=stran(:,rindx);
            else
                stran(:,m)=-enu.*stran(:,rindx);
            end
        end
    else
        for m=ndi+1:ntens
            if m==rindx
                stran(:,m)=stran(:,rindx);
            else
                stran(:,m)=0;
            end
        end
    end
    for m=1:time % Local -> Global
        stran(m,:)=rotmat(ang, stran(m,:), 1);
    end
    dstran=zeros(size(stran));
    for i=2:time
        for j=1:ntens
            dstran(i,j)=stran(i,j)-stran(i-1,j);
        end
    end
%--------------------------------------------------------------------------
% [1.2] Run UMAT subroutine
    statev=zeros([1,nstatv]);
    [strainUMAT, stressUMAT]=UMAT_MEX(props, statev, dstran, ang, s);
    % Global -> Local coordinate system
    for i=1:time
        strainUMAT(i,:)=rotmat(-ang, strainUMAT(i,:),1);
        stressUMAT(i,:)=rotmat(-ang, stressUMAT(i,:),2);
    end
    [stressSIM0, indx]=unique(stressUMAT(:,1),'stable');
    stressSIM=zeros(size(stressSIM0));
    strainSIM=zeros(size(stressSIM0));
    for i=1:size(indx,1)
        stressSIM(i)=stressUMAT(indx(i),1);
        strainSIM(i)=strainUMAT(indx(i),1);
    end
%--------------------------------------------------------------------------
% [1.3] Descretize the simulation results
    indx2=1;
    indx=ones([nstep,1]);
    Tol=1e-5;
    for i=1:nstep
        % Find index for the given boundary conditions
        for j=indx2:max(size(strainSIM))
            if abs(bc(i)-strainSIM(j))<=Tol
                indx(i)=j;
                break;
            end
        end
        indx2=indx(i);
    end
    if k==1
        EXPdata=cell([nstep,nsim]);
        SIMdata=cell([nstep,nsim]);
    end
%--------------------------------------------------------------------------
% [1.4] Interpolation of simulated stress with regard to experimental strain data points
    gap2=zeros([nstep+1,1]);
    for n=1:nstep
        EXPdata{n,k}=MATdata{n};
        if n==1
            SIMdata{n,k}(:,1)= interp1(strainSIM(1:indx(n)),strainSIM(1:indx(n)),EXPdata{n,k}(:,1));
            SIMdata{n,k}(:,2)= interp1(strainSIM(1:indx(n)),stressSIM(1:indx(n)),EXPdata{n,k}(:,1));
        else
            SIMdata{n,k}(:,1)= interp1(strainSIM(indx(n-1):indx(n)),strainSIM(indx(n-1):indx(n)),EXPdata{n,k}(:,1));
            SIMdata{n,k}(:,2)= interp1(strainSIM(indx(n-1):indx(n)),stressSIM(indx(n-1):indx(n)),EXPdata{n,k}(:,1));
        end
% [1.5] Calculation of statistical error
        npoints= round(max(size(SIMdata{n,k}))*0.3,0);
        gap2(n)= weight(1)*cal_error(EXPdata{n,k}(3:end-2,2),SIMdata{n,k}(3:end-2,2),error_par);
        gap2(end)=weight(2)*cal_error(EXPdata{n,k}(end-npoints+1:end,2),SIMdata{n,k}(end-npoints+1:end,2),error_par); % End points
    end

    weigt2=1/nsim;
    gap=gap+weigt2*(sum(gap2));
end
%--------------------------------------------------------------------------
%   [1.6] Plot the simulation results
if plotmode==1
    figure(52)
    hold off
    for k=1:nsim
        EXPresult=cat(1,EXPdata{:,k});
        SIMresult=cat(1,SIMdata{:,k});        
        plot(EXPresult(:,1),EXPresult(:,2),'.-b') % EXP: 3p
        if k==1
            hold on
            grid on
        end
        plot(EXPresult(:,1),SIMresult(:,2),'.-r') % HAH: 3p  
    end
    legend('EXP','SIM','Location','Best')
    axis([min(EXPresult(:,1))*1.1 max(EXPresult(:,1))*1.1 min(EXPresult(:,2))*1.1 max(EXPresult(:,2))*1.1])
    xlabel('True strain') % x-axis label
    ylabel('True stress (MPa)') % y-axis label
end
end
%% [2] Cross-loading
if mode==12
%--------------------------------------------------------------------------
% [2.0] DEFINE PROPS
if props(1)==1 % Kinemtaic hardeing
    props(21:21+xdim-1)=X;
elseif floor(props(1))==2 % Yoshida-Uemori
    props(21:21+xdim-1)=X;
elseif props(1)==3 % HAH11
    props(23:23+xdim-1)=X;
elseif props(1)==4 % HAH14
    props(28:28+xdim-1)=X;
elseif floor(props(1))==5 % HAH20
    props(28:28+xdim-1)=X;
end
%--------------------------------------------------------------------------
for k=1:nsim % Number of simulations
    load(filemat{k});
%--------------------------------------------------------------------------
%   [2.1] SET UP BOUNDARY CONDITION: DSTRAN
    nstep=size(ang,2);
    for i=1:nstep
        bc0=bc(i)*1.5;
        dis=linspace(0,bc0,1/dt);
        time= size(dis,2);        
        stran=zeros([time,ntens]);
        stran(1:time,rindx)=dis';
        stran=log(1+stran);
        for m=1:ndi
            if m==rindx
                stran(:,m)=stran(:,rindx);
            else
                stran(:,m)=-enu.*stran(:,rindx);
            end
        end
        for m=1:time
            stran(m,:)=rotmat(ang(i), stran(m,:), 1);
        end
        dstran=zeros(size(stran));
        for m=2:size(stran,1)
            for n=1:ntens
                dstran(m,n)=stran(m,n)-stran(m-1,n);
            end
        end
%   [2.2] RUN SIMULATION: Find true strain level corresponding to pre-strain
        if i==1
            statev=zeros([1,nstatv]);
        else
            statev=statevUMAT(indx,:);
        end
        [strainUMAT, stressUMAT, statevUMAT]=UMAT_MEX(props, statev, dstran, ang(i), s);
        for m=1:time
            if abs(statevUMAT(m,1)-bc(i)) < 1e-4
                indx=m;
                break;
            else
                indx=time;
            end
        end
    end
%--------------------------------------------------------------------------
%   [2.3] Simulation data post-processing
    % Global -> Local
    for i=1:indx
        stressUMAT(i,:)=trans_principal(stressUMAT(i,:),2);
        stressSIM0(i,:)=trans_principal(stressUMAT(i,:),1);
    end
    % Filter the redundant data
    [strainSIM0, indx2]=unique(statevUMAT(1:indx,1),'stable');
    stressSIM=zeros(size(strainSIM0));
    strainSIM=zeros(size(strainSIM0));
    for i=1:size(indx2,1)
        stressSIM(i)=stressSIM0(indx2(i),1);
        strainSIM(i)=strainSIM0(i,1);
    end
    if k==1
        EXPdata=cell([1,nsim]);
        SIMdata=cell([1,nsim]);
    end
%--------------------------------------------------------------------------
%   [2.4] Interpolation of simulated stress with regard to experimental strain data points
    EXPdata{k}=MATdata;
    SIMdata{k}(:,1)= interp1(strainSIM,strainSIM,EXPdata{k}(:,1)); % Interpolated strain
    SIMdata{k}(:,2)= interp1(strainSIM,stressSIM,EXPdata{k}(:,1)); % Interpolated stress
%--------------------------------------------------------------------------
%   [2.5] Calculation of statistical error
    npoints= round(max(size(SIMdata{k}))*0.3,0);
    gap1=weight(1)*cal_error(SIMdata{k}(4:npoints,2),EXPdata{k}(4:npoints,2),error_par); % Early stage of flow curve
    gap2=weight(2)*cal_error(SIMdata{k}(npoints+1:end-5,2),EXPdata{k}(npoints+1:end-5,2),error_par);
 
    weight3=1/nsim;
    gap=gap+weight3*(gap1+gap2);
%--------------------------------------------------------------------------
%   [2.6] Plot figure for compariosn between experimental and simulation data
if plotmode==1
    figure(52);
    set(gcf, 'Position',  [500, 400, 500*nsim, 400])
    subplot(1,nsim,k);
    hold off
    plot(EXPdata{k}(:,1),EXPdata{k}(:,2),'.-b') % Experimental data
    hold on
    grid on
    plot(EXPdata{k}(:,1),SIMdata{k}(:,2),'.-r') % Simulation data
    legend('EXP','SIM','Location','Best')
    axis([0.9*bc(1) 1.1*bc(2) 0.9*min(EXPdata{k}(:,2)) 1.1*max(EXPdata{k}(:,2))])
    title(filemat{k})
    xlabel('Effective plastic strain') % x-axis label
    ylabel('True stress (MPa)') % y-axis label
end
end
end
%% [3] Anistropic yield function identification
if mode==20
%   [3.0] DEFINE PROPS
if props(2)==2 % Hill 1948
    props(32:32+xdim-1)=X;
elseif props(2)==3 % Yld2000-2d
    props(33:33+xdim-1)=X;
elseif props(2)==4 % Yld2004-18p
    props(33:33+xdim-1)=X;
else
    props(33:33+xdim-1)=X;
end
load(filemat{1});
%--------------------------------------------------------------------------
%   [3.1] Simulate anisotropy of material
SIMdata=zeros([92,3]); % index / angle / flow_stress / r-value
for k=0:91
    if k <= 90
        ang=k;
    else
        ang=0;
        s(2)=1;
    end
    % [3.1.1] Set up boundary condition
    bc0=0.05;
    time=1/dt;
    stran=zeros([time,ntens]);
    statev=zeros([1, nstatv]);
    stran(1:time,rindx)=linspace(0,bc0,time);
    stran=log(1+stran);
    if k <= 90 % Uniaxial stress state
        for m=1:ndi
            if m==rindx
                stran(:,m)=stran(:,rindx);
            else
                stran(:,m)=-enu.*stran(:,rindx);
            end
        end
    else % Biaxial stress state
        stran(:,1)=stran(:,rindx);
        stran(:,2)=stran(:,rindx);
        if ntens==6
            stran(:,3)=-2*enu/(1-enu).*stran(:,rindx);
        end
    end
    for m=1:time % Local -> Global
        stran(m,:)=rotmat(ang, stran(m,:), 1);
    end
    dstran=stran(end,:)-stran(end-1,:);
    % [3.1.2] Run simulation
    [strainUMAT, stressUMAT, statevUMAT]=UMAT_MEX(props, statev, dstran, ang, s);
    % [3.1.3] Restore angle
    SIMdata(k+1,1)=ang;
    % [3.1.4] Normalized flow stress
    sig_bar=statevUMAT(end,2);
    ps=trans_principal(stressUMAT(end,:), 1);
    SIMdata(k+1,2)=ps(1)/sig_bar;
    % [3.1.5] R-value
    dfds=statevUMAT(end, 59:58+ntens);
    ple=rotmat(-ang,dfds,1);
    if k <= 90
        SIMdata(k+1,3)=-ple(2)/(ple(1)+ple(2));
    else
        SIMdata(k+1,3)=dfds(2)/dfds(1);
    end
end
%--------------------------------------------------------------------------
%   [3.3] Calculation of statistical error
nrows=max(size(EXPdata));
indx={};
if nrows==4 % 0, 45, 90, b
    indx{1}=[0, 45, 90, 91];
    indx{2}=[0, 45, 90, 91];
    indx{3}=[1, 2, 3, 4];
    indx{4}=[1, 2, 3, 4];
else
    indx{1}=[0, 45, 90, 91];
    indx{2}=[15, 30, 60, 75];
    indx{3}=[1, 4, 7, 8];
    indx{4}=[2, 3, 5, 6];
end
indx{1}=indx{1}+1;
indx{2}=indx{2}+1;

gap1= weight(1)*cal_error(SIMdata(indx{1},2), EXPdata(indx{3},2), error_par);
gap2= weight(2)*cal_error(SIMdata(indx{1},3), EXPdata(indx{3},3), error_par);
gap3= weight(1)*cal_error(SIMdata(indx{2},2), EXPdata(indx{4},2), error_par);
gap4= weight(2)*cal_error(SIMdata(indx{2},3), EXPdata(indx{4},3), error_par);
gap=(gap1+gap2)+weight(3)*(gap3+gap4);
%--------------------------------------------------------------------------
%   [3.4] Plot the simulation results
if plotmode==1
    figure(52);
    set(gcf, 'Position',  [500, 400, 500*nsim, 400])
    hold off
    yyaxis left
    h1=plot(EXPdata(1:end-1,1),EXPdata(1:end-1,2),'--o');
    axis([0 90 min(EXPdata(1:end-1,2))*0.95 max(EXPdata(1:end-1,2))*1.05]);
    set(h1, 'markerfacecolor', get(h1, 'color'),'markersize', 7);
    xticks([0:15:90]);
    xlabel('Angle') % x-axis label
    ylabel('Normalized flow stress (MPa)') % y-axis label
    hold on
    grid on
    plot(SIMdata(1:end-1,1),SIMdata(1:end-1,2),'-')

    yyaxis right
    hold off
    h2=plot(EXPdata(1:end-1,1),EXPdata(1:end-1,3),'--s');
    axis([0 90 min(EXPdata(1:end-1,3))*0.95 max(EXPdata(1:end-1,3))*1.05]);
    set(h2, 'markerfacecolor', get(h2, 'color') ,'markersize', 8);
    hold on
    plot(SIMdata(1:end-1,1),SIMdata(1:end-1,3),'-')
    ylabel('R-value') % y-axis label
    legend('EXP','SIM','Location','Best')
end
end
