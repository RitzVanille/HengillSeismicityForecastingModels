    % Model result visualisation 

close all
clear all
    % Path to data 
addpath('../Data/')
    % Path to models
addpath('../ETAS0/')
addpath('../ETAS_advanced/')
addpath('../ETAS_hydro/')
addpath('../SeismogenicIndex/')
    % Colour loading and appearance specifications
batlowW=[1 1 1;flip(crameri('batlow'))];
red1=[0.702 0 0.267]; orange1=[0.937 0.545 0]; purple1=[0.424 0.082 0.428]; 
green1=[0.361 0.518 0.275]; blue1=[0.29 0.471 0.808];
    %Load seismic data
load('SeismicityRate_recorded.mat');   
    %Geographic coordinates 
% Hengill region (for eqarthquakes catalogue analysis)
HengillLon=[-22 -20.8]; HengillLat=[63.8 64.2]; % Limits of the larger field
% Active region (reinjection and production) % Limits of the active region to model
ActiveRegionLon=[-21.5 -21.25]; % Longitude in degrees East (negative for degrees West)
ActiveRegionLat=[63.98 64.08]; % Latitude in degrees North (negative for degrees South)
Gridgranularity=[0.005 0.005]; % grid granularity with increment in longitude and latitude
[LonGMat,LatGMat,areamat,grid_vec_format]=GlobRectGrid_ver3(Gridgranularity,ActiveRegionLon,ActiveRegionLat); % Create grid


    %SI model
load('FP_times_SI')
        %Comparison total recorded/expected in time
load('SeismicityRate_SI.mat');
NRecorded=zeros(length(FP_times),1);
SumRateForecasted_SI=zeros(length(FP_times),1);
for i=1:length(SeismicityRate_SI)
   SumRateForecasted_SI(i)=sum(SeismicityRate_SI(:,i),'omitnan');
   NRecorded(i)=sum(SeismicityRate_recorded(:,i),'omitnan');
end
figure(4)
plot(1:length(FP_times),NRecorded)   
hold on
plot(1:length(FP_times),SumRateForecasted_SI)
xlabel('Forecast number')
ylabel('Number of events above Mc')
legend('Recorded','Simulated with Seismogenic Index model')
        % Log likelihood in time
load('LLe_SI.mat');
LLe_timeseries=zeros(length(FP_times),1);
for i=1:length(FP_times_SI)
    LLe_timeseries(i)=sum(LLe_SI(:,i),'omitnan'); %sum of spatial LLe for every time step
end
figure(5)
p5SI=plot(FP_times,LLe_timeseries,'Color',red1,'LineWidth'.2)
xlabel('Date')
ylabel('Log likelihood of models')
legend([p5SI],'Seismogenic Index Model')
        % Seismicity rate maps
forecastNum=110; %forecast number to plot (to fill manually)    
figure(6)
pcolor(LonGMat,LatGMat,log10(reshape(SeismicityRate_SI(:,forecastNum),[size(LonGMat,1),size(LonGMat,2)])))
ylim([63.98 64.08])
xlim([-21.5 -21.25])
xlabel('Longitude')
ylabel('Latitude')
colormap(batlowW)
c=colorbar;
c.Label.String='log_{10}(Seismicity rate (M \geq M_{c}))'
axis tight equal
title(datestr(FP_times(forecastNum)))
