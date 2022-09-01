%--------------------------------------------------------------------------
% Development log
% Code writer: Seong-Yong Yoon
% E-mail: theysy@postech.ac.kr
% Ver: 2.0
% Advisor: Frederic Barlat
% Affiliation: Pohang university of sicence and technology (POSTECH)
%--------------------------------------------------------------------------
% This matlab script process the anisotropy of a material in terms of 
% (1) Normlaized flow stress and (2) R-value
clear
close all
clc
addpath 'Data'
%% Read and write data
filemat='AA2090_DATA';
% Read material anisotropy
% The angle of Biaxial data is defined as ang=91.
data=readmatrix('AA2090_YLD.csv');
EXPdata(:,1)=data(:,1); % Angle
EXPdata(end,1)=91; % Biaxial Angle
EXPdata(:,2)=data(:,2); % Normalized flow stress
EXPdata(:,3)=data(:,3); % R-value
%% Plot the pre-processed data
figure(51);
set(gcf, 'Position',  [500, 400, 500, 400])
hold off
yyaxis left
h1=plot(EXPdata(1:end-1,1),EXPdata(1:end-1,2),'-o');
axis([0 90 min(EXPdata(1:end-1,2))*0.9 max(EXPdata(1:end-1,2))*1.1]);
% set(h1, 'markerfacecolor', get(h1, 'color'));
xticks([0:15:90]);
xlabel('Angle') % x-axis label
ylabel('Normalized flow stress (MPa)') % y-axis label
hold on
grid on

yyaxis right
hold off
h2=plot(EXPdata(1:end-1,1),EXPdata(1:end-1,3),'-o');
%     set(gca,'XColor','black','YColor','color');
axis([0 90 min(EXPdata(1:end-1,3))*0.9 max(EXPdata(1:end-1,3))*1.1]);
% set(h2, 'markerfacecolor', get(h2, 'color'));
hold on
ylabel('R-value') % y-axis label
%% Save the experimental data
save(append('Data\',filemat), 'EXPdata')
disp('#Message: Date Pre-processing is Done!!!');