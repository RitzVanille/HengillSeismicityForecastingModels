    % Model comparison 

close all
clear all
    % Path to data 
addpath('../Data/')
    % Path to models
addpath('../ETAS_models/ETAS0/Outputs')
% addpath('../ETAS_advanced/')
addpath('../SeismogenicIndex/')
    % Colour loading and appearance specifications
batlowW=[1 1 1;flip(crameri('batlow'))];
red1=[0.702 0 0.267]; orange1=[0.937 0.545 0]; purple1=[0.424 0.082 0.428]; 
green1=[0.361 0.518 0.275]; blue1=[0.29 0.471 0.808];

    % N-tests    
load('FP_times_SI')
load('SeismicityRate_recorded.mat');   
NRecorded=zeros(length(FP_times),1);
%Seismogenic Index   
load('SeismicityRate_SI.mat');
SumRateForecasted_SI=zeros(length(FP_times),1);
for i=1:length(FP_times)
   SumRateForecasted_SI(i)=sum(SeismicityRate_SI(:,i),'omitnan');
   NRecorded(i)=sum(SeismicityRate_recorded(:,i),'omitnan');
end
        % plot N-test for a forecast    
i=20; % to plot select forecast time index i. For example:
disp(['Plotting N-test for ' datestr(FP_times(i))])
figure(1)
x = 0:200;
y = poisspdf(x,SumRateForecasted_SI(i));
bar(x,y,'EdgeColor',blue1,'FaceColor',blue1)
xline(NRecorded(i),'k','LineWidth',2)
hold on
plot([poissinv(0.025,SumRateForecasted_SI(i)),poissinv(0.975,SumRateForecasted_SI(i))],[0.025 0.025],'g|-','LineWidth',2)
title(datestr(FP_times(i)))
xlabel('Number of earthquakes (M \geq M_{c})')
ylabel('PDF')


    % Information gain over time
%Seismogenic Index   
        % Log likelihood in time
load('LLe_SI.mat');
LLe_timeseries_SI=zeros(length(FP_times),1);
for i=1:length(FP_times)
    LLe_timeseries_SI(i)=sum(LLe_SI(:,i),'omitnan'); %sum of spatial LLe for every time step
end
figure(5)
p5SI=plot(FP_times,LLe_timeseries_SI,'Color',red1,'LineWidth',2)
xlabel('Date')
ylabel('Log likelihood of models')
        % Cumulative log likelihood in time
Cumul_LLe_timeseries_SI=zeros(length(FP_times),1);
Cumul_LLe_timeseries_SI(1)=LLe_timeseries_SI(1);
for i=2:length(FP_times)
    Cumul_LLe_timeseries_SI(i)=Cumul_LLe_timeseries_SI(i-1)+LLe_timeseries_SI(i); %sum of spatial LLe for every time step cumulative
end
figure(50)
p50SI=plot(FP_times,Cumul_LLe_timeseries_SI,'Color',red1,'LineWidth',2)
xlabel('Date')
ylabel('Log likelihood of models')

    %ETAS0 model
load('FP_ETAS_Times.mat');
load('LLe_ETAS_spatial.mat');
        % Log likelihood in time

LLe_timeseries_ETAS0=zeros(length(FP_ETAS_Times),1);
for i=1:length(FP_ETAS_Times)
    LLe_timeseries_ETAS0(i)=sum(LLe_ETAS_spatial(:,i),'omitnan'); %sum of spatial LLe for every time step
end
figure(5)
hold on
p5E0=plot(FP_ETAS_Times,LLe_timeseries_ETAS0,'Color',green1,'LineWidth',2)
xlabel('Date')
ylabel('Log likelihood of models')
legend([p5SI,p5E0],'Seismogenic Index Model','ETAS0 Model')
legend([],'Seismogenic Index Model')
        % Cumulative log likelihood in time
Cumul_LLe_timeseries_ETAS0=zeros(length(FP_ETAS_Times),1);
Cumul_LLe_timeseries_ETAS0(1)=LLe_timeseries_ETAS0(1);
for i=2:length(FP_ETAS_Times)
    Cumul_LLe_timeseries_ETAS0(i)=Cumul_LLe_timeseries_ETAS0(i-1)+LLe_timeseries_ETAS0(i); %sum of spatial LLe for every time step cumulative
end
figure(50)
hold on
p50E0=plot(FP_ETAS_Times,Cumul_LLe_timeseries_ETAS0,'Color',green1,'LineWidth',2)
xlabel('Date')
ylabel('Log likelihood of models')
legend([p50SI,p50E0],'Seismogenic Index Model','ETAS0 Model')

