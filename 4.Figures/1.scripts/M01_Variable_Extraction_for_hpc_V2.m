%pouya 04/12/17
%extracts all variables from wrfout (and combines restarted ones) and saves
%them as matlab (.mat) files
clear 
close all
clc
%% display availble varialbes (of interst)
for domainname={'wrfout_d03'}

%[~,date_time,~,~,~,~,~,~,~,variables_to_display_name]=pouya_netcdf_read_restart_V2(domainname{1},'XLAT');
variables_to_display_name= {'SMOIS' 'ACGRDFLX' 'ACHFX' 'ACLHF' 'ACLWDNB' 'ACLWUPB' 'ACSWDNB' 'ACSWUPB' 'XLAT' 'XLONG' 'V10' 'U10' 'LU_INDEX' 'Q2' 'T2' 'TH2' 'PSFC' 'IVGTYP' 'VEGFRA' 'GRDFLX' 'FRC_URB2D' 'UTYPE_URB' 'LAI' 'HGT' 'TSK' 'SWDOWN' 'GLW' 'ALBEDO' 'PBLH' 'HFX' 'QFX' 'LH' 'SST'};
% all variables
varnumberS=[1:1:length(variables_to_display_name)];

%% load data (all variables, all runs) and convert to timeseries    

for iii=1%9:10%:length(varnumberS)
var_i=varnumberS(iii);
fprintf('%s\n',variables_to_display_name{var_i});
%load data
[data,date_time,~,~,~,vardesc,varunit,~,~]=pouya_netcdf_read_restart_V2(domainname{1},variables_to_display_name{var_i});
%save data
savename=sprintf('data_%s_%s',domainname{1},variables_to_display_name{var_i});
if size(data,4)==1
save     (savename,'data');
else
save     (savename,'data','-v7.3');
end
%save unit and description
savename=sprintf('data_%s_%s_vardesc',domainname{1},variables_to_display_name{var_i});
save     (savename,'vardesc');
savename=sprintf('data_%s_%s_varunit',domainname{1},variables_to_display_name{var_i});
save     (savename,'varunit');
%save date_time once
if iii==1
savename=sprintf('data_%s_%s',domainname{1},'date_time');
save     (savename,'date_time');
end

end

%eval(sprintf('delete %s_*',domainname{1}));
%eval(sprintf('zip(''data_%s'',{''data_%s_*.mat''});',domainname{1},domainname{1}));
%eval(sprintf('delete data_%s_*.mat',domainname{1}));
end
%dirname=pwd;
%mkdir(dirname(length(dirname)-25:length(dirname)))
%movefile('*.mat',(dirname(length(dirname)-25:length(dirname))))

%make an empthy text file to signal matlab succeded
fopen('matlabSucceded','wt+')
