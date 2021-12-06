    % Model comparison 

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

    %Selection of models to compare
SIpath='../SeismogenicIndex/'; 
load('FP_times_SI')
ETAS0='../ETAS_models/ETAS0/';
ETASa='../ETAS_models/ETAS_advanced/';
ETASh='../ETAS_models/ETAS_hydro/';


    
    % N-tests
load('SeismicityRate_recorded.mat');   
NRecorded=zeros(length(FP_times),1);
%Seismogenic Index   
load('SeismicityRate_SI.mat');
SumRateForecasted_SI=zeros(length(FP_times),1);
for i=1:length(SeismicityRate_SI)
   SumRateForecasted_SI(i)=sum(SeismicityRate_SI(:,i),'omitnan');
   NRecorded(i)=sum(SeismicityRate_recorded(:,i),'omitnan');
   figure
   x = 0:200;
   y = poisspdf(x,SumRateForecasted_SI(i));
   bar(x,y,'EdgeColor',blue1,'FaceColor',blue1)
   xline(NRecorded(i),'k','LineWidth',2)
   hold on
   plot([poissinv(0.025,SumRateForecasted_SI(i)),poissinv(0.975,SumRateForecasted_SI(i))],[0.025 0.025],'g|-','LineWidth',2)
   title(datestr(FP_times(i)))
   xlabel('Number of earthquakes (M \geq M_{c})')
   ylabel('PDF')
end


    % Information gain over time
%Seismogenic Index   
load('LLe_SI.mat');
LLe_timeseries=zeros(length(FP_times),1);
for i=1:length(FP_times_SI)
    LLe_timeseries(i)=sum(LLe_SI(:,i),'omitnan'); %sum of spatial LLe for every time step
end

% 
%     %Comparison total recorded/expected in time
% figure
% plot(1:12,NumReal2)   
% hold on
% plot(1:12,SumRateForecasted)
% xlabel('forecast number')
% ylabel('number of events above Mc')
% legend('recorded','simulated')

% %% Monitoring point for each region
% load('Wells.mat');
% WellHeadCoordinates.zone=categorical(WellHeadCoordinates.zone);
% Husmuli_MP=[mean(WellHeadCoordinates.Lon(WellHeadCoordinates.zone=='Húsmúli')) mean(WellHeadCoordinates.Lat(WellHeadCoordinates.zone=='Húsmúli'))];
% Hverahlid_MP=[mean(WellHeadCoordinates.Lon(WellHeadCoordinates.zone=='Hverahlíd')) mean(WellHeadCoordinates.Lat(WellHeadCoordinates.zone=='Hverahlíd'))];
% Hellishedi_MP=[mean(WellHeadCoordinates.Lon(WellHeadCoordinates.zone=='Hellisheði')) mean(WellHeadCoordinates.Lat(WellHeadCoordinates.zone=='Hellisheði'))];
% Sleggja_MP=[mean(WellHeadCoordinates.Lon(WellHeadCoordinates.zone=='Sleggja')) mean(WellHeadCoordinates.Lat(WellHeadCoordinates.zone=='Sleggja'))];
% Grauhnukar_MP=[mean(WellHeadCoordinates.Lon(WellHeadCoordinates.zone=='Gráuhnúkar')) mean(WellHeadCoordinates.Lat(WellHeadCoordinates.zone=='Gráuhnúkar'))];
% Skardmyrarfjall_MP=[mean(WellHeadCoordinates.Lon(WellHeadCoordinates.zone=='Skardmyrarfjall')) mean(WellHeadCoordinates.Lat(WellHeadCoordinates.zone=='Skardmyrarfjall'))];
% 
% Husmuli_MP_Isnet=[mean(WellHeadCoordinates.ISNETX(WellHeadCoordinates.zone=='Húsmúli')) mean(WellHeadCoordinates.ISNETY(WellHeadCoordinates.zone=='Húsmúli'))];
% Hverahlid_MP_Isnet=[mean(WellHeadCoordinates.ISNETX(WellHeadCoordinates.zone=='Hverahlíd')) mean(WellHeadCoordinates.ISNETY(WellHeadCoordinates.zone=='Hverahlíd'))];
% Hellishedi_MP_Isnet=[mean(WellHeadCoordinates.ISNETX(WellHeadCoordinates.zone=='Hellisheði')) mean(WellHeadCoordinates.ISNETY(WellHeadCoordinates.zone=='Hellisheði'))];
% Sleggja_MP_Isnet=[mean(WellHeadCoordinates.ISNETX(WellHeadCoordinates.zone=='Sleggja')) mean(WellHeadCoordinates.ISNETY(WellHeadCoordinates.zone=='Sleggja'))];
% Grauhnukar_MP_Isnet=[mean(WellHeadCoordinates.ISNETX(WellHeadCoordinates.zone=='Gráuhnúkar')) mean(WellHeadCoordinates.ISNETY(WellHeadCoordinates.zone=='Gráuhnúkar'))];
% Skardmyrarfjall_MP_Isnet=[mean(WellHeadCoordinates.ISNETX(WellHeadCoordinates.zone=='Skardmyrarfjall')) mean(WellHeadCoordinates.ISNETY(WellHeadCoordinates.zone=='Skardmyrarfjall'))];
% 
% 
% for i=1:size(LatGMat,1)*size(LatGMat,2)
%     Lonx(i)=LonGMat(i);
%     Laty(i)=LatGMat(i);
% end
% KnotsGrid=[Lonx' Laty'];
% [Ihus,Dishus]=knnsearch(KnotsGrid,Husmuli_MP,'K',1);
% [rowhus colhus]=find(LatGMat==Laty(Ihus) & LonGMat==Lonx(Ihus));
% [Igra,Disgra]=knnsearch(KnotsGrid,Grauhnukar_MP,'K',1);
% [rowgra colgra]=find(LatGMat==Laty(Igra) & LonGMat==Lonx(Igra));
% [Iskar,Disskar]=knnsearch(KnotsGrid,Skardmyrarfjall_MP,'K',1);
% [rowskar colskar]=find(LatGMat==Laty(Iskar) & LonGMat==Lonx(Iskar));
% [Ihver,Dishver]=knnsearch(KnotsGrid,Hverahlid_MP,'K',1);
% [rowhver colhver]=find(LatGMat==Laty(Ihver) & LonGMat==Lonx(Ihver));
% [Isleg,Dissleg]=knnsearch(KnotsGrid,Sleggja_MP,'K',1);
% [rowsleg colsleg]=find(LatGMat==Laty(Isleg) & LonGMat==Lonx(Isleg));
% [Ihell,Dishell]=knnsearch(KnotsGrid,Hellishedi_MP,'K',1);
% [rowhell colhell]=find(LatGMat==Laty(Ihell) & LonGMat==Lonx(Ihell));
% 
% for i=1:length(FP_times)
%    Nrec_hus(i)= RateSeismicity(15,22,i);
%    Nrec_skar(i)= RateSeismicity(rowskar,colskar,i);
%    Nrec_hell(i)= RateSeismicity(rowskar,colskar,i);
%    Nrec_hver(i)= RateSeismicity(rowhver,colhver,i);
%    Nrec_gra(i)= RateSeismicity(rowgra,colgra,i);
%    Nrec_sleg(i)= RateSeismicity(rowsleg,colsleg,i);
%    Nexp_hus(i)= Nexpected(15,22,i);
%    Nexp_skar(i)= Nexpected(rowskar,colskar,i);
%    Nexp_hell(i)= Nexpected(rowskar,colskar,i);
%    Nexp_hver(i)= Nexpected(rowhver,colhver,i);
%    Nexp_gra(i)= Nexpected(rowgra,colgra,i);
%    Nexp_sleg(i)= Nexpected(rowsleg,colsleg,i);
% end
% 
% figure
% subplot(321)
% plot(1:length(FP_times),Nrec_sleg,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(FP_times),Nexp_sleg,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Number of events')
% title('Sleggja')
% subplot(322)
% plot(1:length(FP_times),Nrec_skar,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(FP_times),Nexp_skar,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Number of events')
% title('Skarðsmýrarfjall')
% subplot(323)
% plot(1:length(FP_times),Nrec_hus,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(FP_times),Nexp_hus,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Number of events')
% title('Húsmúli')
% subplot(324)
% plot(1:length(FP_times),Nrec_hell,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(FP_times),Nexp_hell,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Number of events')
% title('Hellisheiði')
% subplot(325)
% plot(1:length(FP_times),Nrec_gra,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(FP_times),Nexp_gra,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Number of events')
% title('Gráuhnúkar')
% subplot(326)
% plot(1:length(FP_times),Nrec_hver,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(FP_times),Nexp_hver,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Number of events')
% title('Hverahlíð')
% 
% %  Fit test
% 
%     % Recorded number of events
% Fit_times=linspace(LPstart,LPend,15);
% Nrecorded=zeros(size(LonGMat,1),size(LonGMat,2),length(Fit_times));
% for i=1:length(Fit_times)
%     for jlon=1:size(LonGMat,2)
%         for klat=1:size(LonGMat,1)
%             % Find events in the cell with mag above Mc and up to time
%             IndexMag=Cat_time(:)<Fit_times(i) & Cat(:,2)>=Mc ...
%                 & Cat(:,3)>=LatQext(klat) & Cat(:,3)<LatQext(klat+1) ...
%                 & Cat(:,4)>=LonQext(jlon) & Cat(:,4)<LonQext(jlon+1);
%             Nrecorded(klat,jlon,i)=sum(IndexMag);
%         end
%     end
% end
% 
% Nfit=zeros(size(LonGMat,1),size(LonGMat,2),length(Fit_times));
% for i=1:length(Fit_times)
%    Nfit(:,:,i)=MapCumulVolumes(:,:,i).*10.^(Sigma(:,:,1)-bvalue(:,:,1).*Mc);
% end
% for i=1:length(Fit_times)
%    fitNrec_hus(i)= Nrecorded(15,22,i);
%    fitNrec_skar(i)= Nrecorded(rowskar,colskar,i);
%    fitNrec_hell(i)= Nrecorded(rowskar,colskar,i);
%    fitNrec_hver(i)= Nrecorded(rowhver,colhver,i);
%    fitNrec_gra(i)= Nrecorded(rowgra,colgra,i);
%    fitNrec_sleg(i)= Nrecorded(rowsleg,colsleg,i);
%    fitNexp_hus(i)= Nfit(15,22,i);
%    fitNexp_skar(i)= Nfit(rowskar,colskar,i);
%    fitNexp_hell(i)= Nfit(rowskar,colskar,i);
%    fitNexp_hver(i)= Nfit(rowhver,colhver,i);
%    fitNexp_gra(i)= Nfit(rowgra,colgra,i);
%    fitNexp_sleg(i)= Nfit(rowsleg,colsleg,i);
% end
% figure
% subplot(321)
% plot(1:length(Fit_times),fitNrec_sleg,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(Fit_times),fitNexp_sleg,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Cumulative # events')
% title('Sleggja')
% subplot(322)
% plot(1:length(Fit_times),fitNrec_skar,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(Fit_times),fitNexp_skar,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Cumulative # events')
% title('Skarðsmýrarfjall')
% subplot(323)
% plot(1:length(Fit_times),fitNrec_hus,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(Fit_times),fitNexp_hus,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Cumulative # events')
% title('Húsmúli')
% subplot(324)
% plot(1:length(Fit_times),fitNrec_hell,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(Fit_times),fitNexp_hell,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Cumulative # events')
% title('Hellisheiði')
% subplot(325)
% plot(1:length(Fit_times),fitNrec_gra,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(Fit_times),fitNexp_gra,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Cumulative # events')
% title('Gráuhnúkar')
% subplot(326)
% plot(1:length(Fit_times),fitNrec_hver,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(Fit_times),fitNexp_hver,'Color',red1,'LineWidth',2)
% legend('Recorded','Expected','Location','northwest')
% ylabel('Cumulative # events')
% title('Hverahlíð')
% 
% 
% % fit check Húsmúli
% [Ihus,Dishus]=knnsearch(KnotsGrid,Husmuli_MP,'K',9);
% fit9Nrec_hus=zeros(9,length(Fit_times));
% fit9Nexp_hus=zeros(9,length(Fit_times));
% for i=1:9
%     [rowhus(i),colhus(i)]=find(LatGMat==Laty(Ihus(i)) & LonGMat==Lonx(Ihus(i)));
%     for j=1:length(Fit_times)
%         fit9Nrec_hus(i,j)= Nrecorded(rowhus(i)+1,colhus(i)+1,j);
%         fit9Nexp_hus(i,j)= Nfit(rowhus(i)+1,colhus(i)+1,j);
%     end
%      figure(9)
% %     title('Húsmúli 9 closest cells')
%     subplot(3,3,i)
%     plot(1:length(Fit_times),fit9Nrec_hus(i,:),'Color',green1,'LineWidth',2)
%     hold on
%     plot(1:length(Fit_times),fit9Nexp_hus(i,:),'Color',red1,'LineWidth',2)
%     ylabel('Cumulative # events')
%     title(['a_{fb}=',num2str(Sigma(rowhus(i),colhus(i),1)),', b=',num2str(bvalue(rowhus(i),colhus(i),1))])
% 
% end
% for i=1:length(Fit_times)
%     sumRec(i)=sum(fit9Nrec_hus(:,i));
%     sumExp(i)=sum(fit9Nexp_hus(:,i));
% end
% figure
% plot(1:length(Fit_times),sumRec,'Color',green1,'LineWidth',2)
% hold on
% plot(1:length(Fit_times),sumExp,'Color',red1,'LineWidth',2)
% ylabel('Cumulative # events')
% legend('Recorded','Expected','Location','northwest')
% title('Húsmúli sum of 9 cells')


