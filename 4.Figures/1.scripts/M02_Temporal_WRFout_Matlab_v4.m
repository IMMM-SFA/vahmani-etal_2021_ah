%Pouya, Feb 19 2021
%convert wrf outputs from .mat to .txt for EnergyPlus
clear 
close all
clc
addpath('/Volumes/GoogleDrive/My Drive/PROJECTS/2016_Ultra_High_Res_Climate_Modeling_LDRD_and_MLA/STUDY 2018 Climate Change Extreme Heat Exposure CA/2.WRF_PostProcessing/0.Pouya_matlab_fx')
%%
% folder that processed data is saved
runname1='LA_2009_WRF_AH0';plot_runname1='AH0'        ;
runname2='LA_2009_WRF_AH1';plot_runname2='AH1'        ;

%timestep for the start (hr 1) of first and last COMPLETE days
firstDay_start=7
lastDay_start=967

Foldername='/Volumes/PVahmani_G_RAID/Projects/STUDY_2020_Anthropogenic_Heat_Feedback';
Foldername_Eplus='/Volumes/GoogleDrive/My Drive/PROJECTS/2021_IM3_Phase_2/STUDY 2020 Anthropogenic Heat Feedback/2.Eplus';
Savefolder='/Volumes/GoogleDrive/My Drive/PROJECTS/2021_IM3_Phase_2/STUDY 2020 Anthropogenic Heat Feedback/1.WRF/2.Figures';

foldername1=sprintf('%s/%s',Foldername,runname1);
foldername2=sprintf('%s/%s',Foldername,runname2);
foldername1_Eplus=sprintf('%s/1.%s',Foldername_Eplus,runname1);
foldername2_Eplus=sprintf('%s/1.%s',Foldername_Eplus,runname2);
savefolder=sprintf('%s/M02_Temporal_%s-%s',Savefolder,runname1,runname2);
if ~exist(savefolder,'dir')
    mkdir(savefolder)
end

% Variables
Variables={'AH_RElEXH' 'AH_REJ' 'AH' 'CDH' 'LAI' 'HGT' 'ALBEDO' 'LU_INDEX' 'IVGTYP' 'VEGFRA' 'FRC_URB2D' 'UTYPE_URB' 'T2' 'HFX' 'LH' 'SMOIS1' 'GRDFLX' 'SWDOWN' 'GLW' 'WINDS' 'PBLH' 'RH' 'Q2' 'TSK' 'TH2' 'PSFC' 'QFX'}';% 'SST' 'ACLWDNB' 'ACLWUPB' 'ACSWDNB' 'ACSWUPB'
%%
%load lat, long, and lu_index
filename=sprintf('%s/data_wrfout_d03_%s',foldername1,'XLAT');
m=matfile(filename);data_XLAT=m.data(:,:,1);clear m filename
filename=sprintf('%s/data_wrfout_d03_%s',foldername1,'XLONG');
m=matfile(filename);data_XLONG=m.data(:,:,1);clear m filename
filename=sprintf('%s/data_wrfout_d03_%s',foldername1,'LU_INDEX');
m=matfile(filename);data_LU_INDEX=m.data(:,:,1);clear m filename
filename=sprintf('%s/data_wrfout_d03_%s',foldername1,'FRC_URB2D');
m=matfile(filename);data_FRC_URB2D=m.data(:,:,1);clear m filename
filename=sprintf('%s/data_wrfout_d03_%s',foldername1,'UTYPE_URB');
m=matfile(filename);data_UTYPE_URB=m.data(:,:,1);clear m filename
%load date time
filename=sprintf('%s/data_wrfout_d03_date_time',foldername1);
m=matfile(filename);date_time=m.date_time;clear m filename
%extract the complete days only
fprintf('only complete days are used in the analysis: timestep %d to %d\n',firstDay_start,lastDay_start+23);
date_time=date_time(firstDay_start:lastDay_start+23,:);

%% loading data
for var_i=1%1:length(Variables)
varname=Variables{var_i};

for run_i=1:2
foldername=eval(sprintf('foldername%d',run_i));
foldername_Eplus=eval(sprintf('foldername%d_Eplus',run_i));

if strcmp(varname,'AH_REJ')
   fprintf('AH_REJ is based on AH_REJ parameters calculated from E+\n');
   data_VAR=nan(size(data_LU_INDEX,1),size(data_LU_INDEX,2),size(date_time,1));
   % load AH_REJ paramters
   for k=1:24:size(date_time,1)
   % load data (all variables, all runs)
   filename4=sprintf('%s/Eplus_%s_Building_%04d-%02d-%02d.txt',foldername_Eplus,varname,date_time(k,1),date_time(k,2),date_time(k,3));
   fid  = fopen(filename4) ;
   data = textscan(fid,'AH:  %f %f %f');
   data_AH = cell2mat(data);
   data = textscan(fid,'AHDIUPRF: %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
   data_AHDIUPRF = cell2mat(data);
   fclose(fid);
   for kk=1:24
       AH_1=data_AH(1,1)*data_AHDIUPRF(1,kk);
       AH_2=data_AH(1,2)*data_AHDIUPRF(1,kk);
       AH_3=data_AH(1,3)*data_AHDIUPRF(1,kk); 
       data_VAR_temp=nan(size(data_LU_INDEX));
       data_VAR_temp(data_LU_INDEX==24)=AH_1*data_FRC_URB2D(data_LU_INDEX==24);
       data_VAR_temp(data_LU_INDEX==25)=AH_2*data_FRC_URB2D(data_LU_INDEX==25);
       data_VAR_temp(data_LU_INDEX==26)=AH_3*data_FRC_URB2D(data_LU_INDEX==26);
       data_VAR(:,:,k+kk-1)=data_VAR_temp;
   end
   clear data_AH* data
   end
   %load var desc, and unit
   varunit='W m-2';
   vardesc='AH_REJ from Eplus';

elseif strcmp(varname,'AH_RElEXH')
   fprintf('AH_RElEXH is based on AH_RElEXH parameters calculated from E+\n');
   data_VAR=nan(size(data_LU_INDEX,1),size(data_LU_INDEX,2),size(date_time,1));
   % load AH_RElEXH paramters
   for k=1:24:size(date_time,1)
   % load data (all variables, all runs)
   filename4=sprintf('%s/Eplus_%s_Building_%04d-%02d-%02d.txt',foldername_Eplus,varname,date_time(k,1),date_time(k,2),date_time(k,3));
   fid  = fopen(filename4) ;
   data = textscan(fid,'AH:  %f %f %f');
   data_AH = cell2mat(data);
   data = textscan(fid,'AHDIUPRF: %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
   data_AHDIUPRF = cell2mat(data);
   fclose(fid);
   for kk=1:24
       AH_1=data_AH(1,1)*data_AHDIUPRF(1,kk);
       AH_2=data_AH(1,2)*data_AHDIUPRF(1,kk);
       AH_3=data_AH(1,3)*data_AHDIUPRF(1,kk); 
       data_VAR_temp=nan(size(data_LU_INDEX));
       data_VAR_temp(data_LU_INDEX==24)=AH_1*data_FRC_URB2D(data_LU_INDEX==24);
       data_VAR_temp(data_LU_INDEX==25)=AH_2*data_FRC_URB2D(data_LU_INDEX==25);
       data_VAR_temp(data_LU_INDEX==26)=AH_3*data_FRC_URB2D(data_LU_INDEX==26);
       data_VAR(:,:,k+kk-1)=data_VAR_temp;
   end
   clear data_AH* data
   end
   %load var desc, and unit
   varunit='W m-2';
   vardesc='AH_RElEXH from Eplus';

elseif strcmp(varname,'AH')
   fprintf('AH is based on AH parameters calculated from E+\n');
   data_VAR=nan(size(data_LU_INDEX,1),size(data_LU_INDEX,2),size(date_time,1));
   % load AH paramters
   for k=1:24:size(date_time,1)
   % load data (all variables, all runs)
   filename4=sprintf('%s/Eplus_%s_Building_%04d-%02d-%02d.txt',foldername_Eplus,varname,date_time(k,1),date_time(k,2),date_time(k,3));
   fid  = fopen(filename4) ;
   data = textscan(fid,'AH:  %f %f %f');
   data_AH = cell2mat(data);
   data = textscan(fid,'AHDIUPRF: %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
   data_AHDIUPRF = cell2mat(data);
   fclose(fid);
   for kk=1:24
       AH_1=data_AH(1,1)*data_AHDIUPRF(1,kk);
       AH_2=data_AH(1,2)*data_AHDIUPRF(1,kk);
       AH_3=data_AH(1,3)*data_AHDIUPRF(1,kk); 
       data_VAR_temp=nan(size(data_LU_INDEX));
       data_VAR_temp(data_LU_INDEX==24)=AH_1*data_FRC_URB2D(data_LU_INDEX==24);
       data_VAR_temp(data_LU_INDEX==25)=AH_2*data_FRC_URB2D(data_LU_INDEX==25);
       data_VAR_temp(data_LU_INDEX==26)=AH_3*data_FRC_URB2D(data_LU_INDEX==26);
       data_VAR(:,:,k+kk-1)=data_VAR_temp;
   end
   clear data_AH* data
   end
   %load var desc, and unit
   varunit='W m-2';
   vardesc='AH from Eplus';

elseif strcmp(varname,'CDH')
   fprintf('CDH is calcualted using T2\n');
   % load componenets
   filename4=sprintf('%s/data_wrfout_d03_%s',foldername,'T2');
   m=matfile(filename4);data_VAR=m.data;clear m filename*
   %calculate CDH: https://www.sciencedirect.com/science/article/pii/B9780128034767000088
   t_min=26;
   data_VAR=data_VAR-273.15; %convert K to C
   data_VAR=data_VAR-26;
   data_VAR(data_VAR<0)=0;
   vardesc='CDH from T2';
   varunit='degree-hour';
    
elseif strcmp(varname,'WINDS')
   fprintf('WINDS is calcualted using u10 and v10\n');
   % load wind componenet data
   filename4=sprintf('%s/data_wrfout_d03_%s',foldername,'U10');
   m=matfile(filename4);data_U10=m.data;clear m filename*
   filename4=sprintf('%s/data_wrfout_d03_%s',foldername,'V10');
   m=matfile(filename4);data_V10=m.data;clear m filename*
   %convert to wind vector magnitudes
   data_VAR=hypot(data_U10,data_V10);clear data_U10 data_V10
   vardesc='wind speed (U10^2+V10^2)^0.5';
   filename=sprintf('%s/data_wrfout_d03_%s_varunit',foldername,'U10');
   load(filename);clear filename
elseif strcmp(varname,'WINDD')
   fprintf('WINDD is calcualted using u10 and v10\n');
   % load wind componenet data
   filename4=sprintf('%s/data_wrfout_d03_%s',foldername,'U10');
   m=matfile(filename4);data_U10=m.data;clear m filename*
   filename4=sprintf('%s/data_wrfout_d03_%s',foldername,'V10');
   m=matfile(filename4);data_V10=m.data;clear m filename*
   %convert to wind direction (north:0)
   data_VAR=-1*atan2d(data_V10,data_U10)+90;data_VAR=rem(360+data_VAR,360);%clear data_U10 data_V10
   vardesc='wind dir (N:0, E:90)';
   varunit='degree';
elseif strcmp(varname,'SMOIS1')
   fprintf('SMOSI1 is SMOIS first layer\n');
   % load componenets
   filename4=sprintf('%s/data_wrfout_d03_%s',foldername,'SMOIS');
   m=matfile(filename4);data_VAR=squeeze(m.data(:,:,1,:));clear m filename*
   %load var desc, and unit
   filename=sprintf('%s/data_wrfout_d03_%s_varunit',foldername,'SMOIS');
   load(filename);clear filename
   filename=sprintf('%s/data_wrfout_d03_%s_vardesc',foldername,'SMOIS');
   load(filename);clear filename
elseif strcmp(varname,'RH')
   fprintf('RH is calcualted using T2, Q2 and PSFC\n');
   % load componenets
   filename4=sprintf('%s/data_wrfout_d03_%s',foldername,'T2');
   m=matfile(filename4);data_T2=m.data;clear m filename*
   filename4=sprintf('%s/data_wrfout_d03_%s',foldername,'Q2');
   m=matfile(filename4);data_Q2=m.data;clear m filename*
   filename4=sprintf('%s/data_wrfout_d03_%s',foldername,'PSFC');
   m=matfile(filename4);data_PSFC=m.data;clear m filename*
   %calculate RH: https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
   pq0=379.90516;
   a2=17.2693882;
   a3=273.16;
   a4=35.86;
   data_VAR=data_Q2./((pq0./data_PSFC).*exp(a2.*(data_T2-a3)./(data_T2-a4)));clear data_Q2 data_T2 data_PSFC
   %    %approach2:calculate RH: http://mailman.ucar.edu/pipermail/wrf-users/2012/002546.html
   %    data_VAR=0.263.*data_PSFC.*data_Q2.*(exp(17.67.*(data_T2-273.16)./(data_T2-29.65))).^-1;clear data_Q2 data_T2 data_PSFC
   vardesc='RH from T2,Q2&PSFC';
   varunit='%';
else
    % load data (all variables, all runs)
    filename4=sprintf('%s/data_wrfout_d03_%s',foldername,varname);
    m=matfile(filename4);data_VAR=m.data;clear m filename*
    %load var desc, and unit
    filename=sprintf('%s/data_wrfout_d03_%s_varunit',foldername,varname);
    load(filename);clear filename
    filename=sprintf('%s/data_wrfout_d03_%s_vardesc',foldername,varname);
    load(filename);clear filename
end

%convert K to C
if strcmp(varunit,'K')
data_VAR=data_VAR-273.15; 
varunit='C';
end

%extract the complete days only
fprintf('only complete days are used in the analysis: timestep %d to %d\n',firstDay_start,lastDay_start+23);
if ~strcmp(varname,'AH') && ~strcmp(varname,'AH_RElEXH') && ~strcmp(varname,'AH_REJ')%AH is already for complete days'AH_RElEXH' 'AH_REJ' 
    data_VAR=data_VAR(:,:,firstDay_start:lastDay_start+23);
end
%%
%time series
data_VAR_timeseries=nan(size(data_VAR,3),3);
for k=1:size(data_VAR,3)
%     fprintf('processing %d-%d-%d-%d:%d:%d\n',date_time(k,1),date_time(k,2),date_time(k,3),date_time(k,4),date_time(k,5),date_time(k,6));
    data_VAR_temp=data_VAR(:,:,k);
    data_VAR_timeseries(k,1)=mean(data_VAR_temp(data_LU_INDEX==24));
    data_VAR_timeseries(k,2)=mean(data_VAR_temp(data_LU_INDEX==25));
    data_VAR_timeseries(k,3)=mean(data_VAR_temp(data_LU_INDEX==26));
end
clear data_VAR_temp k 

%daily timeseries
data_VAR_daily_mean=nan(size(data_VAR,3)/24,3);
data_VAR_daily_max =nan(size(data_VAR,3)/24,3);
data_VAR_daily_min =nan(size(data_VAR,3)/24,3);
data_VAR_diurnal   =nan(size(data_VAR,3)/24,24,3);
day_i=0;
for k=1:24:size(data_VAR,3)
    day_i=day_i+1;
    data_VAR_temp=data_VAR(:,:,k:k+23);
%     fprintf('processing %d-%d-%d-%d:%d:%d\n',date_time(k,1),date_time(k,2),date_time(k,3),date_time(k,4),date_time(k,5),date_time(k,6));
    %daily mean
    data_VAR_temp_mean=mean(data_VAR_temp,3);
    data_VAR_daily_mean(day_i,1)=nanmean(data_VAR_temp_mean(data_LU_INDEX==24));
    data_VAR_daily_mean(day_i,2)=nanmean(data_VAR_temp_mean(data_LU_INDEX==25));
    data_VAR_daily_mean(day_i,3)=nanmean(data_VAR_temp_mean(data_LU_INDEX==26));
    %daily max
    data_VAR_temp_max=max(data_VAR_temp,[],3);
    data_VAR_daily_max(day_i,1)=nanmean(data_VAR_temp_max(data_LU_INDEX==24));
    data_VAR_daily_max(day_i,2)=nanmean(data_VAR_temp_max(data_LU_INDEX==25));
    data_VAR_daily_max(day_i,3)=nanmean(data_VAR_temp_max(data_LU_INDEX==26));
    %daily min
    data_VAR_temp_min=min(data_VAR_temp,[],3);
    data_VAR_daily_min(day_i,1)=nanmean(data_VAR_temp_min(data_LU_INDEX==24));
    data_VAR_daily_min(day_i,2)=nanmean(data_VAR_temp_min(data_LU_INDEX==25));
    data_VAR_daily_min(day_i,3)=nanmean(data_VAR_temp_min(data_LU_INDEX==26));
    %diurnal profile
    for kk=1:24
        data_VAR_diurnal_temp=data_VAR_temp(:,:,kk);
        data_VAR_diurnal(day_i,kk,1)=nanmean(data_VAR_diurnal_temp(data_LU_INDEX==24));
        data_VAR_diurnal(day_i,kk,2)=nanmean(data_VAR_diurnal_temp(data_LU_INDEX==25));
        data_VAR_diurnal(day_i,kk,3)=nanmean(data_VAR_diurnal_temp(data_LU_INDEX==26));        
    end
end
clear *temp* kk k day_i

eval(sprintf('data_VAR%d=data_VAR;clear data_VAR',run_i));
eval(sprintf('data_VAR_timeseries%d=data_VAR_timeseries;clear data_VAR_timeseries',run_i));
eval(sprintf('data_VAR_daily_mean%d=data_VAR_daily_mean;clear data_VAR_daily_mean',run_i));
eval(sprintf('data_VAR_daily_max%d =data_VAR_daily_max ;clear data_VAR_daily_max' ,run_i));
eval(sprintf('data_VAR_daily_min%d =data_VAR_daily_min ;clear data_VAR_daily_min' ,run_i));
eval(sprintf('data_VAR_diurnal%d   =data_VAR_diurnal   ;clear data_VAR_diurnal'   ,run_i));

end
%% stats
%land cover land use of interest
hiUrb=26;meUrb=25;loUrb=24;
agric=38;
water=17;
land=100;land_nonurb=101;land_nonag=102;land_nonurb_nonag=103;
all=1000;
%set the land cover land use of interest for each domain
LULC ={hiUrb;meUrb;loUrb;land_nonurb};clear hiUrb meUrb loUrb agric water land land_nonurb land_nonag land_nonurb_nonag all
LULCl={'hiUrb';'meUrb';'loUrb';'land_nonurb'};

% statistics
stat_base= pouya_time_series2( data_VAR1, data_LU_INDEX, LULC, LULCl );
stat_est = pouya_time_series2( data_VAR2, data_LU_INDEX, LULC, LULCl );

%calc changes
delta=(data_VAR2-data_VAR1);
% statistics
stat_delta= pouya_time_series2( delta, data_LU_INDEX, LULC, LULCl );

%calc percent changes
delta_perc=(data_VAR2-data_VAR1)./data_VAR1*100;
%avoid unrealistic percentage changes over close to zero baselines (this mostly affect error bars)
delta_perc(data_VAR1<=prctile(data_VAR1(:),0.1))=NaN;
fprintf('note: %s < %f (%s) were not condidered in perc change calc\n',varname,prctile(data_VAR1(:),0.1),varunit);
% statistics
stat_delta_perc= pouya_time_series2( delta_perc,data_LU_INDEX, LULC, LULCl );

%% plot timeseries
%hourly
figure('position',[500 500 1200 420])
x=1:size(data_VAR_timeseries1,1);
y11=data_VAR_timeseries1(:,1); y12=data_VAR_timeseries1(:,2); y13=data_VAR_timeseries1(:,3);
y21=data_VAR_timeseries2(:,1); y22=data_VAR_timeseries2(:,2); y23=data_VAR_timeseries2(:,3);
yyaxis left
plot(x,y11,':','Color',[0      0.4470 0.7410],'LineWidth',2);hold on
plot(x,y12,':','Color',[0.4660 0.6740 0.1880],'LineWidth',2);hold on
plot(x,y13,':','Color',[0.8500 0.3250 0.0980],'LineWidth',2);hold on
plot(x,y21,'-','Color',[0      0.4470 0.7410],'LineWidth',2);hold on
plot(x,y22,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',2);hold on
plot(x,y23,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',2);
ylabel(sprintf('%s (%s)',varname,varunit))
yl = get(gca,'ylim');
ylim([yl(1,1)-(yl(1,2)-yl(1,1))/4 yl(1,2)])
hold on
y30=zeros(size(y13));
y31=data_VAR_timeseries2(:,1)-data_VAR_timeseries1(:,1);
y32=data_VAR_timeseries2(:,2)-data_VAR_timeseries1(:,2);
y33=data_VAR_timeseries2(:,3)-data_VAR_timeseries1(:,3);
yyaxis right
plot(x,y31,'-','Color',[0      0.4470 0.7410],'LineWidth',1);hold on
plot(x,y32,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',1);hold on
plot(x,y33,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1);hold on
plot(x,y30,'-','Color',[0      0      0     ],'LineWidth',0.5);
ylabel(sprintf('%s difference (%s)',varname,varunit))
yl = get(gca,'ylim');
ylim([yl(1,1) yl(1,2)+(yl(1,2)-yl(1,1))*2])
title(sprintf('%s (%s): hourly',varname,varunit))
legend('Utype 1 AH0','Utype 2 AH0','Utype 3 AH0','Utype 1 AH1','Utype 2 AH1','Utype 3 AH1','Utype 1 diff','Utype 2 diff','Utype 3 diff');
legend off
xlim([1 size(data_VAR_timeseries1,1)]);
set(gca,'FontSize',18);
savename=sprintf('%s/temporal_hourly_%s_%s-%s',savefolder,varname,plot_runname2,plot_runname1);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');
%% plot timeseries
%diurnal profile
figure('position',[500 500 600 420])
x=1:size(data_VAR_diurnal1,2);
y11=mean(data_VAR_diurnal1(:,:,1),1); y12=mean(data_VAR_diurnal1(:,:,2),1); y13=mean(data_VAR_diurnal1(:,:,3),1);
y21=mean(data_VAR_diurnal2(:,:,1),1); y22=mean(data_VAR_diurnal2(:,:,2),1); y23=mean(data_VAR_diurnal2(:,:,3),1);
yyaxis left
plot(x,y11,':','Color',[0      0.4470 0.7410],'LineWidth',2);hold on
plot(x,y12,':','Color',[0.4660 0.6740 0.1880],'LineWidth',2);hold on
plot(x,y13,':','Color',[0.8500 0.3250 0.0980],'LineWidth',2);hold on
plot(x,y21,'-','Color',[0      0.4470 0.7410],'LineWidth',2);hold on
plot(x,y22,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',2);hold on
plot(x,y23,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',2);
ylabel(sprintf('%s (%s)',varname,varunit))
yl = get(gca,'ylim');
ylim([yl(1,1)-(yl(1,2)-yl(1,1))/4 yl(1,2)])
hold on
y30=zeros(size(y13));
y31=mean(data_VAR_diurnal2(:,:,1),1)-mean(data_VAR_diurnal1(:,:,1),1);
y32=mean(data_VAR_diurnal2(:,:,2),1)-mean(data_VAR_diurnal1(:,:,2),1);
y33=mean(data_VAR_diurnal2(:,:,3),1)-mean(data_VAR_diurnal1(:,:,3),1);
yyaxis right
plot(x,y31,'-','Color',[0      0.4470 0.7410],'LineWidth',1);hold on
plot(x,y32,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',1);hold on
plot(x,y33,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1);hold on
plot(x,y30,'-','Color',[0      0      0     ],'LineWidth',0.5);
ylabel(sprintf('%s difference (%s)',varname,varunit))
yl = get(gca,'ylim');
ylim([yl(1,1) yl(1,2)+(yl(1,2)-yl(1,1))*2])
title(sprintf('%s (%s): diurnal profile',varname,varunit))
legend('Utype 1 AH0','Utype 2 AH0','Utype 3 AH0','Utype 1 AH1','Utype 2 AH1','Utype 3 AH1','Utype 1 diff','Utype 2 diff','Utype 3 diff');
legend off
xlim([1 size(data_VAR_diurnal1,2)]);
set(gca,'FontSize',18);
savename=sprintf('%s/temporal_diurnal_%s_%s-%s',savefolder,varname,plot_runname2,plot_runname1);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');
%% plot timeseries
%daily mean
figure('position',[500 500 1200 420])
x=1:size(data_VAR_daily_mean1,1);
y11=data_VAR_daily_mean1(:,1); y12=data_VAR_daily_mean1(:,2); y13=data_VAR_daily_mean1(:,3);
y21=data_VAR_daily_mean2(:,1); y22=data_VAR_daily_mean2(:,2); y23=data_VAR_daily_mean2(:,3);
yyaxis left
plot(x,y11,':','Color',[0      0.4470 0.7410],'LineWidth',2);hold on
plot(x,y12,':','Color',[0.4660 0.6740 0.1880],'LineWidth',2);hold on
plot(x,y13,':','Color',[0.8500 0.3250 0.0980],'LineWidth',2);hold on
plot(x,y21,'-','Color',[0      0.4470 0.7410],'LineWidth',2);hold on
plot(x,y22,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',2);hold on
plot(x,y23,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',2);
ylabel(sprintf('%s (%s)',varname,varunit))
yl = get(gca,'ylim');
ylim([yl(1,1)-(yl(1,2)-yl(1,1))/4 yl(1,2)])
hold on
y30=zeros(size(y13));
y31=data_VAR_daily_mean2(:,1)-data_VAR_daily_mean1(:,1);
y32=data_VAR_daily_mean2(:,2)-data_VAR_daily_mean1(:,2);
y33=data_VAR_daily_mean2(:,3)-data_VAR_daily_mean1(:,3);
yyaxis right
plot(x,y31,'-','Color',[0      0.4470 0.7410],'LineWidth',1);hold on
plot(x,y32,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',1);hold on
plot(x,y33,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1);hold on
plot(x,y30,'-','Color',[0      0      0     ],'LineWidth',0.5);
ylabel(sprintf('%s difference (%s)',varname,varunit))
yl = get(gca,'ylim');
ylim([yl(1,1) yl(1,2)+(yl(1,2)-yl(1,1))*2])
title(sprintf('%s (%s): daily mean',varname,varunit))
legend('Utype 1 AH0','Utype 2 AH0','Utype 3 AH0','Utype 1 AH1','Utype 2 AH1','Utype 3 AH1','Utype 1 diff','Utype 2 diff','Utype 3 diff');
legend off
xlim([1 size(data_VAR_daily_mean1,1)]);
set(gca,'FontSize',18);
savename=sprintf('%s/temporal_dailyMean_%s_%s-%s',savefolder,varname,plot_runname2,plot_runname1);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');
%% plot timeseries
%daily max
figure('position',[500 500 1200 420])
x=1:size(data_VAR_daily_max1,1);
y11=data_VAR_daily_max1(:,1); y12=data_VAR_daily_max1(:,2); y13=data_VAR_daily_max1(:,3);
y21=data_VAR_daily_max2(:,1); y22=data_VAR_daily_max2(:,2); y23=data_VAR_daily_max2(:,3);
yyaxis left
plot(x,y11,':','Color',[0      0.4470 0.7410],'LineWidth',2);hold on
plot(x,y12,':','Color',[0.4660 0.6740 0.1880],'LineWidth',2);hold on
plot(x,y13,':','Color',[0.8500 0.3250 0.0980],'LineWidth',2);hold on
plot(x,y21,'-','Color',[0      0.4470 0.7410],'LineWidth',2);hold on
plot(x,y22,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',2);hold on
plot(x,y23,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',2);
ylabel(sprintf('%s (%s)',varname,varunit))
yl = get(gca,'ylim');
ylim([yl(1,1)-(yl(1,2)-yl(1,1))/4 yl(1,2)])
hold on
y30=zeros(size(y13));
y31=data_VAR_daily_max2(:,1)-data_VAR_daily_max1(:,1);
y32=data_VAR_daily_max2(:,2)-data_VAR_daily_max1(:,2);
y33=data_VAR_daily_max2(:,3)-data_VAR_daily_max1(:,3);
yyaxis right
plot(x,y31,'-','Color',[0      0.4470 0.7410],'LineWidth',1);hold on
plot(x,y32,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',1);hold on
plot(x,y33,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1);hold on
plot(x,y30,'-','Color',[0      0      0     ],'LineWidth',0.5);
ylabel(sprintf('%s difference (%s)',varname,varunit))
yl = get(gca,'ylim');
ylim([yl(1,1) yl(1,2)+(yl(1,2)-yl(1,1))*2])
title(sprintf('%s (%s): daily max',varname,varunit))
legend('Utype 1 AH0','Utype 2 AH0','Utype 3 AH0','Utype 1 AH1','Utype 2 AH1','Utype 3 AH1','Utype 1 diff','Utype 2 diff','Utype 3 diff');
legend off
xlim([1 size(data_VAR_daily_max1,1)]);
set(gca,'FontSize',18);
savename=sprintf('%s/temporal_dailyMax_%s_%s-%s',savefolder,varname,plot_runname2,plot_runname1);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');
%% plot timeseries
%daily min
figure('position',[500 500 1200 420])
x=1:size(data_VAR_daily_min1,1);
y11=data_VAR_daily_min1(:,1); y12=data_VAR_daily_min1(:,2); y13=data_VAR_daily_min1(:,3);
y21=data_VAR_daily_min2(:,1); y22=data_VAR_daily_min2(:,2); y23=data_VAR_daily_min2(:,3);
yyaxis left
plot(x,y11,':','Color',[0      0.4470 0.7410],'LineWidth',2);hold on
plot(x,y12,':','Color',[0.4660 0.6740 0.1880],'LineWidth',2);hold on
plot(x,y13,':','Color',[0.8500 0.3250 0.0980],'LineWidth',2);hold on
plot(x,y21,'-','Color',[0      0.4470 0.7410],'LineWidth',2);hold on
plot(x,y22,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',2);hold on
plot(x,y23,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',2);
ylabel(sprintf('%s (%s)',varname,varunit))
yl = get(gca,'ylim');
ylim([yl(1,1)-(yl(1,2)-yl(1,1))/4 yl(1,2)])
hold on
y30=zeros(size(y13));
y31=data_VAR_daily_min2(:,1)-data_VAR_daily_min1(:,1);
y32=data_VAR_daily_min2(:,2)-data_VAR_daily_min1(:,2);
y33=data_VAR_daily_min2(:,3)-data_VAR_daily_min1(:,3);
yyaxis right
plot(x,y31,'-','Color',[0      0.4470 0.7410],'LineWidth',1);hold on
plot(x,y32,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',1);hold on
plot(x,y33,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1);hold on
plot(x,y30,'-','Color',[0      0      0     ],'LineWidth',0.5);
ylabel(sprintf('%s difference (%s)',varname,varunit))
yl = get(gca,'ylim');
ylim([yl(1,1) yl(1,2)+(yl(1,2)-yl(1,1))*2])
title(sprintf('%s (%s): daily min',varname,varunit))
legend('Utype 1 AH0','Utype 2 AH0','Utype 3 AH0','Utype 1 AH1','Utype 2 AH1','Utype 3 AH1','Utype 1 diff','Utype 2 diff','Utype 3 diff');
legend off
xlim([1 size(data_VAR_daily_min1,1)]);
set(gca,'FontSize',18);
savename=sprintf('%s/temporal_dailyMin_%s_%s-%s',savefolder,varname,plot_runname2,plot_runname1);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');
%%
% %save data
% savename=sprintf('%s/plot_data_%s.mat',savefolder,varname);
% save(savename,'varname','varunit','vardesc','data*','date_time','runname*','plot_runname*','delta*','stat_base','stat_est','stat_delta','stat_delta_perc','LULC*','savefolder','-v7.3')
end %var_i
