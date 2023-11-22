%--------------------------------------------------------------------------
% Development log
% Code writer: Seong-Yong Yoon
% E-mail: theysy@postech.ac.kr
% Ver: 2.0
% Advisor: Frederic Barlat
% Affiliation: Pohang university of sicence and technology (POSTECH)
%--------------------------------------------------------------------------
% This matlab script pre-process the experimental data in terms of
% Effective plastic strain and true stress.
clear
close all
clc
addpath 'Data'
%% Define details of the Data pre-processing
dataNo=500;
ang(1)=90; % Direction for pre-strain (1st Loading)
ang(2)=00; % Direction for second-loading (2nd Loading)
filemat='EDDQ_TDT-RDT';

%% equivalent strain - true stress
data0=csvread('EDDQ_TDT-RDT.CSV');
[data_temp, indx0]=unique(data0(:,2),'stable');
data=zeros([max(size(data_temp)), 2]);
for i=1:max(size(indx0))
    data(i,1)=data0(indx0(i),1);
    data(i,2)=data0(indx0(i),2);
end
ndata=max(size(data));
pick=ceil(ndata/dataNo);
MATdata(:,1)=data(1:pick:end,1);
MATdata(:,2)=data(1:pick:end,2);
bc(1)=MATdata(1,1);
bc(2)=max(MATdata(:,1));
plot(MATdata(:,1), MATdata(:,2), '-o');
grid on
%% Save the experimental data
save (append('Data\',filemat), 'ang', 'MATdata', 'bc') 
disp('#Message: Date Pre-processing is Done!!!');

