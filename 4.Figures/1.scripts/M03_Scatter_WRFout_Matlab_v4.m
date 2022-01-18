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
savefolder=sprintf('%s/M03_Scatter_%s-%s',Savefolder,runname1,runname2);
if ~exist(savefolder,'dir')
    mkdir(savefolder)
end

% Variables
% Variables={'AH_RElEXH' 'AH_REJ' 'AHinEB' 'AH' 'CDH' 'LAI' 'HGT' 'ALBEDO' 'LU_INDEX' 'IVGTYP' 'VEGFRA' 'FRC_URB2D' 'UTYPE_URB' 'T2' 'HFX' 'LH' 'SMOIS1' 'GRDFLX' 'SWDOWN' 'GLW' 'WINDS' 'PBLH' 'RH' 'Q2' 'TSK' 'TH2' 'PSFC' 'QFX'}';% 'SST' 'ACLWDNB' 'ACLWUPB' 'ACSWDNB' 'ACSWUPB'
variable_x='T2';
var_att1_x='max';%daily mean max min, (or - for CDH)
var_att2_x='base';%base, est, or change 

variable_y='PBLH';
var_att1_y='mean';%daily mean max min, (or - for CDH)
var_att2_y='base';%base, est, or change 

variable_sz='WINDS';
var_att1_sz='mean';%daily mean max min, (or - for CDH)
var_att2_sz='base';%base, est, or change 

Variables={variable_x variable_y variable_sz};
Var_att1s={var_att1_x var_att1_y var_att1_sz};
Var_att2s={var_att2_x var_att2_y var_att2_sz};
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
%daily date_time
date_time_daily=nan(size(date_time,1)/24,3);
day_i=0;
for k=1:24:size(date_time,1)
    day_i=day_i+1;
    date_time_daily(day_i,:)=date_time(k,1:3);
end;clear *temp* k day_i 

%% loading data
for var_i=1:length(Variables)
varname=Variables{var_i};
var_att1=Var_att1s{var_i};
var_att2=Var_att2s{var_i};

for run_i=1:2
    
if strcmp(var_att2,'base')
    foldername=eval(sprintf('foldername%d',1));
    foldername_Eplus=eval(sprintf('foldername%d_Eplus',1));
else
    foldername=eval(sprintf('foldername%d',run_i));
    foldername_Eplus=eval(sprintf('foldername%d_Eplus',run_i));
end

if strcmp(varname,'AHinEB')%AH percentage of surface Energy Budget
   fprintf('AHinEB: AH/LH+HFX where AH is based on AH parameters calculated from E+\n');
   data_VAR=nan(size(data_LU_INDEX,1),size(data_LU_INDEX,2),size(date_time,1));
   % load AH paramters
   for k=1:24:size(date_time,1)
   % load data (all variables, all runs)
   filename4=sprintf('%s/Eplus_%s_Building_%04d-%02d-%02d.txt',foldername_Eplus,'AH',date_time(k,1),date_time(k,2),date_time(k,3));
   fid  = fopen(filename4);
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
   %daily mean/max/min
    data_VAR_daily=nan(size(data_VAR,3)/24,3);
    day_i=0;
    for k=1:24:size(data_VAR,3)
        day_i=day_i+1;
        data_VAR_temp=data_VAR(:,:,k:k+23);
        if     strcmp(var_att1,'mean')
            data_VAR_temp=nanmean(data_VAR_temp,3);
        else
            fprintf('ERROR: AHinEB max or min doesnt make sense\n');
            break
        end
        data_VAR_daily(day_i,1)=nanmean(data_VAR_temp(data_LU_INDEX==24));
        data_VAR_daily(day_i,2)=nanmean(data_VAR_temp(data_LU_INDEX==25));
        data_VAR_daily(day_i,3)=nanmean(data_VAR_temp(data_LU_INDEX==26));
    end;clear *temp* k day_i data_VAR
    data_VAR_daily_AH=data_VAR_daily;clear data_VAR_daily
    
      % load LH, HFX componenet data
   filename4=sprintf('%s/data_wrfout_d03_%s',foldername,'LH');
   m=matfile(filename4);data_LH=m.data;clear m filename*
   filename4=sprintf('%s/data_wrfout_d03_%s',foldername,'HFX');
   m=matfile(filename4);data_HFX=m.data;clear m filename*
    %daily mean/max/min
    %extract the complete days only
    fprintf('only complete days are used in the analysis: timestep %d to %d\n',firstDay_start,lastDay_start+23);
    data_LH    =data_LH    (:,:,firstDay_start:lastDay_start+23);
    data_HFX   =data_HFX   (:,:,firstDay_start:lastDay_start+23);
    data_VAR_daily_LH    =nan(size(data_LH    ,3)/24,3);
    data_VAR_daily_HFX   =nan(size(data_HFX   ,3)/24,3);
    day_i=0;
    for k=1:24:size(data_LH,3)
        day_i=day_i+1;
        data_VAR_temp_LH    =data_LH    (:,:,k:k+23);
        data_VAR_temp_HFX   =data_HFX   (:,:,k:k+23);
        if     strcmp(var_att1,'mean')
            data_VAR_temp_LH    =mean(data_VAR_temp_LH    ,3);
            data_VAR_temp_HFX   =mean(data_VAR_temp_HFX   ,3);
        else
            fprintf('ERROR: AHinEB max or min doesnt make sense\n');
            break
        end
        data_VAR_daily_LH    (day_i,1)=nanmean(data_VAR_temp_LH    (data_LU_INDEX==24));
        data_VAR_daily_LH    (day_i,2)=nanmean(data_VAR_temp_LH    (data_LU_INDEX==25));
        data_VAR_daily_LH    (day_i,3)=nanmean(data_VAR_temp_LH    (data_LU_INDEX==26));
        data_VAR_daily_HFX   (day_i,1)=nanmean(data_VAR_temp_HFX   (data_LU_INDEX==24));
        data_VAR_daily_HFX   (day_i,2)=nanmean(data_VAR_temp_HFX   (data_LU_INDEX==25));
        data_VAR_daily_HFX   (day_i,3)=nanmean(data_VAR_temp_HFX   (data_LU_INDEX==26));
    end;clear *temp* k day_i data_LH data_HFX 

   %total surface avaialbe energy: LH+HFX+GRDFLX (note that GRDFLX is
   %negative during the day)
   data_VAR_daily=data_VAR_daily_AH./(data_VAR_daily_LH+data_VAR_daily_HFX)*100;
   clear data_VAR_daily_AH data_VAR_daily_LH data_VAR_daily_HFX
   vardesc='AH/(LH+HFX)*100';
   varunit='%';
    
   
elseif strcmp(varname,'AH_REJ')
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
   %daily mean/max/min
    data_VAR_daily=nan(size(data_VAR,3)/24,3);
    day_i=0;
    for k=1:24:size(data_VAR,3)
        day_i=day_i+1;
        data_VAR_temp=data_VAR(:,:,k:k+23);
        if     strcmp(var_att1,'mean')
            data_VAR_temp=nanmean(data_VAR_temp,3);
        elseif strcmp(var_att1,'max')
            data_VAR_temp=nanmax (data_VAR_temp,[],3);
        elseif strcmp(var_att1,'min')
            data_VAR_temp=nanmin (data_VAR_temp,[],3);
        end
        data_VAR_daily(day_i,1)=nanmean(data_VAR_temp(data_LU_INDEX==24));
        data_VAR_daily(day_i,2)=nanmean(data_VAR_temp(data_LU_INDEX==25));
        data_VAR_daily(day_i,3)=nanmean(data_VAR_temp(data_LU_INDEX==26));
    end;clear *temp* k day_i data_VAR

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
   %daily mean/max/min
    data_VAR_daily=nan(size(data_VAR,3)/24,3);
    day_i=0;
    for k=1:24:size(data_VAR,3)
        day_i=day_i+1;
        data_VAR_temp=data_VAR(:,:,k:k+23);
        if     strcmp(var_att1,'mean')
            data_VAR_temp=nanmean(data_VAR_temp,3);
        elseif strcmp(var_att1,'max')
            data_VAR_temp=nanmax (data_VAR_temp,[],3);
        elseif strcmp(var_att1,'min')
            data_VAR_temp=nanmin (data_VAR_temp,[],3);
        end
        data_VAR_daily(day_i,1)=nanmean(data_VAR_temp(data_LU_INDEX==24));
        data_VAR_daily(day_i,2)=nanmean(data_VAR_temp(data_LU_INDEX==25));
        data_VAR_daily(day_i,3)=nanmean(data_VAR_temp(data_LU_INDEX==26));
    end;clear *temp* k day_i data_VAR

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
   %daily mean/max/min
    data_VAR_daily=nan(size(data_VAR,3)/24,3);
    day_i=0;
    for k=1:24:size(data_VAR,3)
        day_i=day_i+1;
        data_VAR_temp=data_VAR(:,:,k:k+23);
        if     strcmp(var_att1,'mean')
            data_VAR_temp=nanmean(data_VAR_temp,3);
        elseif strcmp(var_att1,'max')
            data_VAR_temp=nanmax (data_VAR_temp,[],3);
        elseif strcmp(var_att1,'min')
            data_VAR_temp=nanmin (data_VAR_temp,[],3);
        end
        data_VAR_daily(day_i,1)=nanmean(data_VAR_temp(data_LU_INDEX==24));
        data_VAR_daily(day_i,2)=nanmean(data_VAR_temp(data_LU_INDEX==25));
        data_VAR_daily(day_i,3)=nanmean(data_VAR_temp(data_LU_INDEX==26));
    end;clear *temp* k day_i data_VAR

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
   %daily values
    %extract the complete days only
    fprintf('only complete days are used in the analysis: timestep %d to %d\n',firstDay_start,lastDay_start+23);
    data_VAR=data_VAR(:,:,firstDay_start:lastDay_start+23);
    data_VAR_daily=nan(size(data_VAR,3)/24,3);
    day_i=0;
    for k=1:24:size(data_VAR,3)
        day_i=day_i+1;
        data_VAR_temp=data_VAR(:,:,k:k+23);
        data_VAR_temp=sum(data_VAR_temp,3);
        data_VAR_daily(day_i,1)=nanmean(data_VAR_temp(data_LU_INDEX==24));
        data_VAR_daily(day_i,2)=nanmean(data_VAR_temp(data_LU_INDEX==25));
        data_VAR_daily(day_i,3)=nanmean(data_VAR_temp(data_LU_INDEX==26));
    end;clear *temp* kk k day_i data_VAR
    
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
   %daily mean/max/min
    %extract the complete days only
    fprintf('only complete days are used in the analysis: timestep %d to %d\n',firstDay_start,lastDay_start+23);
    data_VAR=data_VAR(:,:,firstDay_start:lastDay_start+23);
    data_VAR_daily=nan(size(data_VAR,3)/24,3);
    day_i=0;
    for k=1:24:size(data_VAR,3)
        day_i=day_i+1;
        data_VAR_temp=data_VAR(:,:,k:k+23);
        if     strcmp(var_att1,'mean')
            data_VAR_temp=mean(data_VAR_temp,3);
        elseif strcmp(var_att1,'max')
            data_VAR_temp=max (data_VAR_temp,[],3);
        elseif strcmp(var_att1,'min')
            data_VAR_temp=min (data_VAR_temp,[],3);
        end
        data_VAR_daily(day_i,1)=nanmean(data_VAR_temp(data_LU_INDEX==24));
        data_VAR_daily(day_i,2)=nanmean(data_VAR_temp(data_LU_INDEX==25));
        data_VAR_daily(day_i,3)=nanmean(data_VAR_temp(data_LU_INDEX==26));
    end;clear *temp* k day_i data_VAR
    
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
   %daily mean/max/min
    %extract the complete days only
    fprintf('only complete days are used in the analysis: timestep %d to %d\n',firstDay_start,lastDay_start+23);
    data_VAR=data_VAR(:,:,firstDay_start:lastDay_start+23);
    data_VAR_daily=nan(size(data_VAR,3)/24,3);
    day_i=0;
    for k=1:24:size(data_VAR,3)
        day_i=day_i+1;
        data_VAR_temp=data_VAR(:,:,k:k+23);
        if     strcmp(var_att1,'mean')
            data_VAR_temp=mean(data_VAR_temp,3);
        elseif strcmp(var_att1,'max')
            data_VAR_temp=max (data_VAR_temp,[],3);
        elseif strcmp(var_att1,'min')
            data_VAR_temp=min (data_VAR_temp,[],3);
        end
        data_VAR_daily(day_i,1)=nanmean(data_VAR_temp(data_LU_INDEX==24));
        data_VAR_daily(day_i,2)=nanmean(data_VAR_temp(data_LU_INDEX==25));
        data_VAR_daily(day_i,3)=nanmean(data_VAR_temp(data_LU_INDEX==26));
    end;clear *temp* k day_i data_VAR
    
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
   %daily mean/max/min
    %extract the complete days only
    fprintf('only complete days are used in the analysis: timestep %d to %d\n',firstDay_start,lastDay_start+23);
    data_VAR=data_VAR(:,:,firstDay_start:lastDay_start+23);
    data_VAR_daily=nan(size(data_VAR,3)/24,3);
    day_i=0;
    for k=1:24:size(data_VAR,3)
        day_i=day_i+1;
        data_VAR_temp=data_VAR(:,:,k:k+23);
        if     strcmp(var_att1,'mean')
            data_VAR_temp=mean(data_VAR_temp,3);
        elseif strcmp(var_att1,'max')
            data_VAR_temp=max (data_VAR_temp,[],3);
        elseif strcmp(var_att1,'min')
            data_VAR_temp=min (data_VAR_temp,[],3);
        end
        data_VAR_daily(day_i,1)=nanmean(data_VAR_temp(data_LU_INDEX==24));
        data_VAR_daily(day_i,2)=nanmean(data_VAR_temp(data_LU_INDEX==25));
        data_VAR_daily(day_i,3)=nanmean(data_VAR_temp(data_LU_INDEX==26));
    end;clear *temp* k day_i data_VAR
    
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
   %daily mean/max/min
    %extract the complete days only
    fprintf('only complete days are used in the analysis: timestep %d to %d\n',firstDay_start,lastDay_start+23);
    data_VAR=data_VAR(:,:,firstDay_start:lastDay_start+23);
    data_VAR_daily=nan(size(data_VAR,3)/24,3);
    day_i=0;
    for k=1:24:size(data_VAR,3)
        day_i=day_i+1;
        data_VAR_temp=data_VAR(:,:,k:k+23);
        if     strcmp(var_att1,'mean')
            data_VAR_temp=mean(data_VAR_temp,3);
        elseif strcmp(var_att1,'max')
            data_VAR_temp=max (data_VAR_temp,[],3);
        elseif strcmp(var_att1,'min')
            data_VAR_temp=min (data_VAR_temp,[],3);
        end
        data_VAR_daily(day_i,1)=nanmean(data_VAR_temp(data_LU_INDEX==24));
        data_VAR_daily(day_i,2)=nanmean(data_VAR_temp(data_LU_INDEX==25));
        data_VAR_daily(day_i,3)=nanmean(data_VAR_temp(data_LU_INDEX==26));
    end;clear *temp* k day_i data_VAR
    
else
    % load data (all variables, all runs)
    filename4=sprintf('%s/data_wrfout_d03_%s',foldername,varname);
    m=matfile(filename4);data_VAR=m.data;clear m filename*
    %load var desc, and unit
    filename=sprintf('%s/data_wrfout_d03_%s_varunit',foldername,varname);
    load(filename);clear filename
    filename=sprintf('%s/data_wrfout_d03_%s_vardesc',foldername,varname);
    load(filename);clear filename
    %daily mean/max/min
    %extract the complete days only
    fprintf('only complete days are used in the analysis: timestep %d to %d\n',firstDay_start,lastDay_start+23);
    data_VAR=data_VAR(:,:,firstDay_start:lastDay_start+23);
    data_VAR_daily=nan(size(data_VAR,3)/24,3);
    day_i=0;
    for k=1:24:size(data_VAR,3)
        day_i=day_i+1;
        data_VAR_temp=data_VAR(:,:,k:k+23);
        if     strcmp(var_att1,'mean')
            data_VAR_temp=mean(data_VAR_temp,3);
        elseif strcmp(var_att1,'max')
            data_VAR_temp=max (data_VAR_temp,[],3);
        elseif strcmp(var_att1,'min')
            data_VAR_temp=min (data_VAR_temp,[],3);
        end
        data_VAR_daily(day_i,1)=nanmean(data_VAR_temp(data_LU_INDEX==24));
        data_VAR_daily(day_i,2)=nanmean(data_VAR_temp(data_LU_INDEX==25));
        data_VAR_daily(day_i,3)=nanmean(data_VAR_temp(data_LU_INDEX==26));
    end;clear *temp* k day_i data_VAR
    
end

%convert K to C
if strcmp(varunit,'K')
data_VAR_daily=data_VAR_daily-273.15; 
varunit='C';
end

%keep data for 2 variables and runs
if var_i==1
eval(sprintf('data_VAR_daily_x_%d=data_VAR_daily;clear data_VAR_daily',run_i));
eval(sprintf('varunit_x=varunit;clear varunit'));
eval(sprintf('vardesc_x=vardesc;clear vardesc'));
elseif var_i==2
eval(sprintf('data_VAR_daily_y_%d=data_VAR_daily;clear data_VAR_daily',run_i));
eval(sprintf('varunit_y=varunit;clear varunit'));
eval(sprintf('vardesc_y=vardesc;clear vardesc'));
elseif var_i==3
eval(sprintf('data_VAR_daily_sz_%d=data_VAR_daily;clear data_VAR_daily',run_i));
eval(sprintf('varunit_sz=varunit;clear varunit'));
eval(sprintf('vardesc_sz=vardesc;clear vardesc'));
end

end

end
%% plot Scatter
%plot data
close all
if strcmp(var_att2_x,'change')
    plot_data_x=data_VAR_daily_x_2-data_VAR_daily_x_1;
elseif strcmp(var_att2_x,'base')
    plot_data_x=data_VAR_daily_x_1;
elseif strcmp(var_att2_x,'est')
    plot_data_x=data_VAR_daily_x_2;
end
if strcmp(var_att2_y,'change')
    plot_data_y=data_VAR_daily_y_2-data_VAR_daily_y_1;
elseif strcmp(var_att2_y,'base')
    plot_data_y=data_VAR_daily_y_1;
elseif strcmp(var_att2_y,'est')
    plot_data_y=data_VAR_daily_y_2;
end
if strcmp(var_att2_sz,'change')
    plot_data_sz=data_VAR_daily_sz_2-data_VAR_daily_sz_1;
elseif strcmp(var_att2_sz,'base')
    plot_data_sz=data_VAR_daily_sz_1;
elseif strcmp(var_att2_sz,'est')
    plot_data_sz=data_VAR_daily_sz_2;
end
%rescale sizes to range 100-200
a=50;b=500;
plot_data_sz=a+(plot_data_sz-min(plot_data_sz(:))).*(b-a)./(max(plot_data_sz(:))-min(plot_data_sz(:)));
%% plot
% figure('position',[500 500 2000 420])
figure
scatter(plot_data_x(:,1),plot_data_y(:,1),plot_data_sz(:,1),'MarkerFaceColor',[0      0.4470 0.7410],'MarkerEdgeColor',[0      0.4470 0.7410]);
hold on
scatter(plot_data_x(:,2),plot_data_y(:,2),plot_data_sz(:,2),'MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880]);
hold on
scatter(plot_data_x(:,3),plot_data_y(:,3),plot_data_sz(:,3),'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980]);

%trendlines
coef = polyfit(plot_data_x(:,1),plot_data_y(:,1), 1);
hline = refline(coef(1), coef(2));
hline.Color = [0      0.4470 0.7410];
hline.LineStyle='--';
coef = polyfit(plot_data_x(:,2),plot_data_y(:,2), 1);
hline = refline(coef(1), coef(2));
hline.Color = [0.4660 0.6740 0.1880];
hline.LineStyle='--';
coef = polyfit(plot_data_x(:,3),plot_data_y(:,3), 1);
hline = refline(coef(1), coef(2));
hline.Color = [0.8500 0.3250 0.0980];
hline.LineStyle='--';

%zero line
if strcmp(var_att2_y,'change')
mu = 0;
hline = refline([0 mu]);
hline.Color = 'k';
hline.LineStyle='-';
end

ylabel(sprintf('%s, %s, daily %s (%s)',variable_y,var_att2_y,var_att1_y,varunit_y))
xlabel(sprintf('%s, %s, daily %s (%s)',variable_x,var_att2_x,var_att1_x,varunit_x))
title(sprintf('%s-%s',runname2,runname1));
legend('Utype 1','Utype 2','Utype 3');
legend off
set(gca,'FontSize',14);

%save figure
savename=sprintf('%s/scatter_%s%sDaily%s-%s%sDaily%s_%s%sDaily%s_%s-%s',savefolder,variable_x,var_att2_x,var_att1_x,variable_y,var_att2_y,var_att1_y,variable_sz,var_att2_sz,var_att1_sz,plot_runname2,plot_runname1);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');


% %save data
% savename=sprintf('%s/plot_data_%s.mat',savefolder,varname);
% save(savename,'varname','varunit','vardesc','data*','date_time','runname*','plot_runname*','delta*','stat_base','stat_est','stat_delta','stat_delta_perc','LULC*','savefolder','-v7.3')
