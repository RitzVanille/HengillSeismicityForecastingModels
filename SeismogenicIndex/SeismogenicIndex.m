    % Seismogenic Index:
close all
clear all
tic
    % Path to data 
addpath('../Data/')

    % Colour loading and appearance specifications
batlowW=[1 1 1;flip(crameri('batlow'))];
red1=[0.702 0 0.267]; orange1=[0.937 0.545 0]; purple1=[0.424 0.082 0.428]; 
green1=[0.361 0.518 0.275]; blue1=[0.29 0.471 0.808];

    %Region definition

% Hengill region (for eqarthquakes catalogue analysis)
HengillLon=[-22 -20.8]; HengillLat=[63.8 64.2]; % Limits of the larger field

% Active region (reinjection and production) % Limits of the active region to model
ActiveRegionLon=[-21.5 -21.25]; % Longitude in degrees East (negative for degrees West)
ActiveRegionLat=[63.98 64.08]; % Latitude in degrees North (negative for degrees South)

Gridgranularity=[0.005 0.005]; % grid granularity with increment in longitude and latitude

[LonGMat,LatGMat,areamat,grid_vec_format]=GlobRectGrid_ver3(Gridgranularity,ActiveRegionLon,ActiveRegionLat); % Create grid
LonQuery=unique(LonGMat); %center of the cells in Longitude
LatQuery=unique(LatGMat); %center of the cells in Latitude
LonQext=[LonQuery-Gridgranularity(1)/2; LonQuery(end)+Gridgranularity(1)/2]; %edges of the cells
LatQext=[LatQuery-Gridgranularity(2)/2; LatQuery(end)+Gridgranularity(2)/2]; %edges of the cells

    % Definition of the learning (LP) and forecasting periods (FP)
LPstart=datetime(2018,12,1); % Start of the learning period
LPend=datetime(2020, 2, 1); % End of the learning period
LPDays=days(LPend-LPstart); %number of days in learning period
FPend=datetime(2021,2,1); %End of the forecasting period
NmFP=days(FPend-LPend); % Number of forecasting period (in days)
FP_times=LPend+transpose(1:NmFP); %times of forecasts produced
FP_days=days(FP_times-LPstart); % time in days since start of LP
load('../Data/TimesVolumes.mat'); % times of flow rate data (daily since start of Learning Period to end of Forecasting Period
idxFp1=find(Times==LPend); % index of end of LP in the hydraulic data stream
TimesInj=days(Times-LPstart);

    % Load hydraulic data
load('../Data/DailyRates.mat'); % daily injection/production rate in kg/s
load('../Data/DailyVolumes.mat'); % daily injected/produced volume in cubic meters
load('../Data/CumulativeDailyVolumes.mat'); %cumulative daily volume in cubic meters

    % Load catalogue Medium Quality Coseismiq Network (see Grigoli et al., 2021 (submitted))
    %format: year, month ,day, hour, minute, second, magnitude, latitude (degrees North), longitude (degrees West), depth (km), event score
RawCat=load('../Data/COSEISMIQ_catalogue'); 
RawCat(:,9)=-RawCat(:,9); %transform longitude to degrees East
Mc=0.3; %Manual input of completeness magnitude from catalogue analysis

    %restrain catalogue to Hengill limits and event score > -5 and with a defined magnitude
Cat_filtered=flip(RawCat(RawCat(:,11)>-5 & ~isnan(RawCat(:,7))...
    & RawCat(:,8)>=ActiveRegionLat(1) & RawCat(:,8)<=ActiveRegionLat(2) ...
    & RawCat(:,9)>=ActiveRegionLon(1) & RawCat(:,9)<=ActiveRegionLon(2),:));
% Extraction of times of events 
Cat_time=datetime([Cat_filtered(:,1:6)]);
% Working catalogue format: time (days since start of LP), magnitude, latitude (degrees North), longitude (degrees East), depth (km), event score
Cat(:,2:6)=Cat_filtered(:,7:11);
Cat(:,1)=days(Cat_time(:,1)-LPstart); %time in days for each event (start of LP = t0)
clear RawCat Cat_filtered 

%% Productivity (a_fb) and b-value optimisation
% Format for all internal calculations: X(number of latitude elements, number of
% longitude elements, number of time steps in the Forecasting Period)
MapRates=reshape(DailyRates,[size(LonGMat,1),size(LonGMat,2),length(Times)]);
MapVolumes=reshape(CumulativeDailyVolume,[size(LonGMat,1),size(LonGMat,2),length(Times)]);

Sigma=zeros(size(LonGMat,1),size(LonGMat,2),length(FP_times));
bvalue=zeros(size(LonGMat,1),size(LonGMat,2),length(FP_times));
RateSeismicity=zeros(size(LonGMat,1),size(LonGMat,2),length(FP_times));
VolumeAvrg=zeros(size(LonGMat,1),size(LonGMat,2),length(FP_times));
RateAvrg2=zeros(size(LonGMat,1),size(LonGMat,2),length(FP_times));

for i=1:NmFP
    for jlon=1:size(LonGMat,2)
        for klat=1:size(LonGMat,1)
            [Icc]=knnsearch(grid_vec_format(:,1:2),[LonGMat(klat,jlon) LatGMat(klat,jlon)],'K',9); %index closest cells
            
            IndexMagFor=Cat(:,1)>TimesInj(idxFp1+i-1) & Cat(:,1)<FP_days(i) & Cat(:,2)>=Mc ...
                      & Cat(:,3)>=LatGMat(klat,jlon)-Gridgranularity(2)/2 & Cat(:,3)<LatGMat(klat,jlon)+Gridgranularity(2)/2 ...
                      & Cat(:,4)>=LonGMat(klat,jlon)-Gridgranularity(1)/2 & Cat(:,4)<LonGMat(klat,jlon)+Gridgranularity(1)/2;
            
            IndexMag=Cat(:,1)<=TimesInj(idxFp1+i-1) & Cat(:,2)>=Mc ...
                      & Cat(:,3)>=LatGMat(klat,jlon)-Gridgranularity(2)/2 & Cat(:,3)<LatGMat(klat,jlon)+Gridgranularity(2)/2 ...
                      & Cat(:,4)>=LonGMat(klat,jlon)-Gridgranularity(1)/2 & Cat(:,4)<LonGMat(klat,jlon)+Gridgranularity(1)/2;

            RatesCC=reshape(MapRates(klat,jlon,1:idxFp1+i-1),1,length(1:idxFp1+i-1));
            VolumesCC=MapVolumes(klat,jlon,idxFp1+i-1);
            % Data for this grid point at this forecasting time in format
            % for optimisation algorithm
            rate=table();
            rate.N=length(Cat(IndexMag,1)); % number of recorded events above Mc 
            rate.N2=length(Cat(IndexMagFor,1)); % number of recorded events above Mc before current time stamp
            rate.m_0=Mc; % completeness magnitude
            rate.t_b_s=transpose(TimesInj(1:idxFp1+i-1)); % time before shut-in [days]
            rate.dot_V_bs=RatesCC*86400*1e-3; % injection/production rate [m^3/day]
            rate.t_sbs=Cat(IndexMag,1)'; % time stamp [days]
            rate.tot_V=VolumesCC; % cumulative volume [m^3]
            rate.m_0_m=Mc-0.05; % magnitude binning correction Dbin/2
            rate.data_magn=Cat(IndexMag,2)'; %  magnitude of events (>=Mc)
            % Optimisation of afb
            f = @(theta)log_lhood_comp(theta,rate);
            theta_0 = [0 1]; % with relaxation time tau: theta_0=[0 24 1];
            options = optimoptions('fminunc','Algorithm','quasi-newton');
            %options.Display = 'iter';
            rateTT{klat,jlon}=rate;
            if rate.dot_V_bs>0
                [theta, fval] = fminunc(f,theta_0,options);
                Sigma(klat,jlon,i)=theta(1);
                bvalue(klat,jlon,i)=theta(2);
            else
                Sigma(klat,jlon,i)=NaN;
                bvalue(klat,jlon,i)=NaN;
            end
            
            RateSeismicity(klat,jlon,i)=rateTT{klat,jlon}.N2;
            VolumeCC2=MapVolumes(klat,jlon,idxFp1+i)-MapVolumes(klat,jlon,idxFp1+i-1);
            VolumeAvrg(klat,jlon,i)=VolumeCC2;
            clear rate
        end
    end
end


%% Computation of forecasted number of event
% Expected number of events
Nexpected=zeros(size(LonGMat,1),size(LonGMat,2),length(FP_times));
SumRateForecasted=zeros(length(FP_times),1);
NRecorded=zeros(length(FP_times),1);
for i=1:length(FP_times)
   Nexpected(:,:,i)=VolumeAvrg(:,:,i).*10.^(Sigma(:,:,i)-bvalue(:,:,i).*Mc);
end
Nexpected(Nexpected==0)=NaN;
SeismicityRate_SI=reshape(Nexpected,[size(LonGMat,1)*size(LonGMat,2), length(FP_times)]);
save('SeismicityRate_SI.mat','SeismicityRate_SI');
SeismicityRate_recorded=reshape(RateSeismicity,[size(LonGMat,1)*size(LonGMat,2), length(FP_times)]);
save('SeismicityRate_recorded.mat','SeismicityRate_recorded');

% %Normalised rate:
% Nnormexp=NaN*ones(size(LonGMat,1),size(LonGMat,2),length(FP_times));
% for i=1:length(FP_times)
%     Nnormexp(:,:,i)=Nexpected(:,:,i)/max(max(Nexpected(:,:,i)));
% end

%% Computation of Log likelihood
LLe_SI=zeros(size(LonGMat,1)*size(LonGMat,2), length(FP_times));
for i=1:length(FP_times)
    for j=1:size(LonGMat,1)*size(LonGMat,2)
        LLe_SI(j,i)=-SeismicityRate_SI(j,i)+SeismicityRate_recorded(j,i)*log10(SeismicityRate_SI(j,i))-log10(factorial(SeismicityRate_recorded(j,i)));
    end    
end
save('LLe_SI.mat','LLe_SI');
save('FP_times_SI.mat','FP_times');

toc
