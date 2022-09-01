function [varargout]= UMAT_MEX(props, statev, delas, ang, s0)
%%   Read state variables
    varargout=cell(1,3);
    [time, ntens]=size(delas);
    nstatv=max(size(statev));
    statev1=zeros(time,nstatv); % Statev(n+1)
    stress=zeros([1,ntens]); % Stress(n)
    stress1=zeros([time,ntens]); % Stress(n+1)
    stress2=zeros([time,3]); % Principal stress
    delas0=zeros(size(delas)); %Elastic strain increment on a local coord sys
    dplas=zeros(size(delas)); % Plastic strain increment
    dstran=zeros(size(delas)); % Total strain increment
    stran=zeros(size(delas)); % Total strain
%   strain increment
    if ntens==3
        ndi=2;
    else
        ndi=3;
    end
    for i=1:time % Global -> Local
        delas0(i,:)=rotmat(-ang,delas(i,:),1);
    end
%   non-zero stress component
    s=rotmat(ang, s0, 2).*1;
    [statev0, stress0]= MML_UMAT(props, statev, s, delas(end,:));
    dfds=statev0(59:58+ntens);
    dfds0=rotmat(-ang,dfds,1);
    sn=zeros([1, ntens]);
    for i=1:ntens
        if s0(i)==1
            rindx=i;
        end
    end
    for i=1:ntens
        if abs(s(i))<1e-8
            sn(i)=0;
        else
            sn(i)=1;
        end
    end
	scale=1/dfds0(rindx);
%%
    for i=1:time
        if statev(4)==1
            dplas(i,:)=delas0(i,rindx)*scale*dfds;
            if rindx>ndi
                dstran(i,:)=delas(i,:);
            else
                dstran(i,:)=dplas(i,:);
            end
        else
            dstran(i,:)=delas(i,:);
        end
        [statev1(i,:), stress1(i,:)]= MML_UMAT(props, statev, stress, dstran(i,:));
        if ang==0 || ang==90
            stress1(i,:)=stress1(i,:).*sn;
        else
            stress1(i,:)=trans_principal(stress1(i,:), 2); % Correction
            stress2(i,:)=trans_principal(stress1(i,:), 1); % Principal stress only for debugging
        end
        
        double=0;
        Tol=1e-5;
        for n=1:ntens
            double=double+abs(dstran(i,n));
        end
        if double<=Tol && i>1
            stress1(i,:)=stress1(i-1,:);
        end
        stress=stress1(i,:);
        statev=statev1(i,:);
        
        if i>1
            stran(i,:)=stran(i-1,:)+dstran(i,:);
        end
    end   

    varargout{1} = stran;
    varargout{2} = stress1;
    varargout{3} = statev1;

%     figure(51);
%     plot(statev1(:,1), statev1(:,40:42),'LineWidth',2);
%     hold on;
%     grid on;
%     plot(statev1(:,1),statev1(:,44:46),'LineWidth',2);
%     legend('\rho_{F}', '\rho_{R+}', '\rho_{R-}', '\rho_{0+}', '\rho_{0-}', '\rho_{TOT}','location','nw');
%     xlabel('Effective strain') % x-axis label
%     ylabel('Dislocation density (m^{-2})') % y-axis label
%     xlim([0 statev1(end,1)])
%     hold off;
    text=1;
    
    
    
