%--------------------------------------------------------------------------
% Development log
% Code writer: Seong-Yong Yoon
% E-mail: theysy@postech.ac.kr
% Ver: 2.0
% Advisor: Frederic Barlat
% Affiliation: Pohang university of sicence and technology (POSTECH)
%--------------------------------------------------------------------------
% This matlab script pre-process the experimental data in terms of
% true strain and true stress.
% Moreover, it is available to process a multiple step experiment such as
% tension-compression.
clear
close all
clc
addpath 'Data'
%% Define material data
filemat='CR980XG3_shear_E20_FR';
ang=0;
filename=append(filemat,'.csv');
data0=csvread(filename);
dataNo=3000; % Selected # of data
mode=2; % 1: Eng -> True | 2: True -> True
bc0=[0.20 -0.29];
% bc0=[0.06 -0.06 0.06];
Tol=2e-3;
%% Set up boundary condition
[data_temp, indx0]=unique(data0(:,2),'stable');
data=zeros([max(size(data_temp)), 2]);
for i=1:max(size(indx0))
    data(i,1)=data0(indx0(i),1);
    data(i,2)=data0(indx0(i),2);
end
if mode==1 % Eng -> True
    data(:,2)=data(:,2).*(1+data(:,1));
    data(:,1)=log(1+data(:,1));
    bc=log(1+bc0);
    bc1=bc;
else
    bc=bc0;
    bc1=exp(bc)-1;
end
nstep=max(size(bc));
ndata=max(size(data));
pick=ceil(ndata/dataNo);
% pick=1;
%% Dividing stress-strain curves
MATdata=cell([1,nstep]); % Raw data [i,j,k]
indx2=1;
indx=ones([nstep,1]);
for k=1:nstep
    % Find index for the given boundary conditions
    for j=indx2:ndata
        if abs(bc(k)-data(j,1))<=Tol
            indx(k)=j;
            break;
        end
    end
    % Define piecewise stress-strain curve
    MATdata{k}=data(indx2:pick:indx(k),:);
    indx2=indx(k);
end
%% Plot the Data pre-processing
lgn=cell(size(nstep));
figure(51);
set(gcf, 'WindowStyle','docked')
for k=1:nstep
    plot(MATdata{k}(:,1),MATdata{k}(:,2), '-o');
    if k==1
        hold on
        grid on
    end
    lgn{k}=strcat('step:', num2str(k));
end
legend(lgn,'Location','Best');
MATdata1=cat(1,MATdata{:});
%% Save experiment data
save (append('Data\',filemat), 'ang', 'MATdata', 'bc') 
disp('#Message: Date Pre-processing is Done!!!');