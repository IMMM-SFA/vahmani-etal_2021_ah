function [data,date_time,dimlen,dimname,MemoryOrder,description,units,stagger,coordinates,Varname]=pouya_netcdf_read_restart(filename_in,varname_in)
% pouya_function: read netcdf and and for a varname returns:
%  data, data_time, dimension lengths, dimension names, MemoryOrder, description, units, stagger, coordinates


listing = dir(sprintf('%s*',filename_in));
for restart_i=1:length(listing)
    filename=sprintf('%s%s',filename_in(1:length(filename_in)-10),listing(restart_i).name);
%open the netcdf
ncid = netcdf.open(filename,'NOWRITE');
% ncdisp(filename)
%
if restart_i==1
%get number of variables
[ndims,nvars,natts,unlimdimID] = netcdf.inq(ncid);clear ndims natts unlimdimID
%
%get information: varnames, var IDs, num of attributes, dimensions (# and name)
for j=1:nvars-1
    
    %useful infor on each varable
    [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,j); 
    Varname  {j,1}=varname;
    VarAtts  (j,1)=varAtts;
    
    % ID of each variable
    varid = netcdf.inqVarID(ncid,varname); 
    Varid(j,1) = varid;
    
    %info on attributes
    for i=1:varAtts
    attname = netcdf.inqAttName(ncid,varid,i-1);%get variable att namies
    Attname  {j,i}=attname;
    attval = netcdf.getAtt(ncid,varid,attname);
    Attval   {j,i}=attval;
    end;
    
    %info on diminstions
    for i=1:size(varDimIDs,2);
    [dimname,dimlen]=netcdf.inqDim(ncid,varDimIDs(1,i)); %use diminstion ID to get dimtion name and length
    Dimname  {j,i}=dimname;
    Dimlen   (j,i)=dimlen;
    end; 
        
end

clear j i dimname dimlen varDimIDs varAtts xtype varname attval attname varid
% load data and info for a specific varname
for j=1:size(Varname,1)
if strcmpi(Varname{j,1},varname_in)
    j_var=j;
    break
end
end
clear j
%get info (data, dimension lengths, dimension name, MemoryOrder, description, units, stagger, coordinates for one variable
for i=1:size(Dimlen,2)
    if Dimlen(j_var,i)~=0
dimlen (1,i)=Dimlen (j_var,i);
dimname(1,i)=Dimname(j_var,i);
    end
end
attname=Attname(j_var,1:VarAtts(j_var));
attval =Attval (j_var,1:VarAtts(j_var));
MemoryOrder='';
description='';
units='';
stagger='';
coordinates='';
for i=2:size(attname,2)
eval(sprintf('%s=''%s'';',attname{1,i},attval{1,i}));
end
clear i
end
%%

varid  =Varid  (j_var,1);
%load variable
data = netcdf.getVar(ncid,varid);
data=double(data);


%%
% time arraye
% number of timesteps
for j=1:size(Dimname,1);for i=1:size(Dimname,2)
if strcmp(Dimname{j,i},'Time');
Time(j,1)=Dimlen(j,i);
end
end;                    end
n_timestep=max(Time);

% start time (character string)
data_starttime = netcdf.getVar(ncid,0);

% the date-time columns
DateString = data_starttime';
% %sometimes the last timestep is not properly saved so this will fix that
if sum(DateString(size(DateString,1),:))==0
   DateString(size(DateString,1),:)=datestr(datevec(addtodate(datenum(datevec(DateString(size(DateString,1)-1,:))),+1,'hour')),'yyyy-mm-dd_HH:MM:SS'); 
end
DateVector = datevec(DateString);
date_time=DateVector; 
clear DateString DateVector data_starttime
% UTC to PST
for j=1:size(date_time,1)
    date_time(j,1:6)=datevec(addtodate(datenum(date_time(j,1:6)),-8,'hour'));
end
clear j
netcdf.close(ncid)
clearvars -except j_var Varid Dimlen Dimname Time data date_time data_* date_time_* dimlen dimname MemoryOrder description units stagger coordinates Varname varname_in restart_i listing filename*

eval(sprintf('data_%d=data;clear data',restart_i));
eval(sprintf('date_time_%d=date_time;clear date_time',restart_i));
end
%combine the restarted files
if length(listing)>1
    for combine_i=1:length(listing)-1
   if combine_i==1
       data_first=data_1;
       date_time_first=date_time_1;
   end

eval(sprintf('data_second=data_%d;',combine_i+1));
eval(sprintf('date_time_second=date_time_%d;',combine_i+1));

clear timestep_start
for search_i=1:length(date_time_first)
    time_to_match=date_time_second(1,:);
        if time_to_match(1)==date_time_first(search_i,1)&&time_to_match(2)==date_time_first(search_i,2)&&time_to_match(3)==date_time_first(search_i,3)&&time_to_match(4)==date_time_first(search_i,4)
        timestep_start=search_i;
    end
end
while ~exist('timestep_start')&&size(date_time_first,1)<8761
%add one timestep to the date_time_first for cases there is no overlap
date_time_first(size(date_time_first,1)+1,:)=datevec(addtodate(datenum(date_time_first(size(date_time_first,1),:)),+1,'hour'));
if size(data_first,3)==1
data_first (size(data_first,1)+1,:)=NaN;
elseif size(data_first,4)==1
data_first (:,:,size(data_first,3)+1)=NaN;
else
data_first (:,:,:,size(data_first,4)+1)=NaN;
end
for search_i=1:length(date_time_first)
    time_to_match=date_time_second(1,:);
        if time_to_match(1)==date_time_first(search_i,1)&&time_to_match(2)==date_time_first(search_i,2)&&time_to_match(3)==date_time_first(search_i,3)&&time_to_match(4)==date_time_first(search_i,4)
        timestep_start=search_i;
    end
end
end

clear timestep_end
for search_i=1:length(date_time_first)
    time_to_match=date_time_second(size(date_time_second,1),:);
        if time_to_match(1)==date_time_first(search_i,1)&&time_to_match(2)==date_time_first(search_i,2)&&time_to_match(3)==date_time_first(search_i,3)&&time_to_match(4)==date_time_first(search_i,4)
        timestep_end=search_i;
    end
end
if ~exist('timestep_end')
date_time_first(timestep_start:timestep_start+size(date_time_second,1)-1,:)=date_time_second(1:size(date_time_second,1),:);
if size(data_first,3)==1
data_first (timestep_start:timestep_start+size(data_second     ,1)-1,:)  =data_second (1:size(data_second     ,1),:);
elseif size(data_first,4)==1
    %irrigation water saved as noahres starts at zero after restart
    if strcmp(varname_in,'NOAHRES')&&sum(sum(data_second(:,:,1)))==0
        last_nonNaN_i=timestep_start;
        while isnan(sum(sum(data_first(:,:,last_nonNaN_i))))
            last_nonNaN_i=last_nonNaN_i-1;
        end
        for noahres_i=1:size(data_second,3)
            data_second(:,:,noahres_i)=data_second(:,:,noahres_i)+data_first(:,:,last_nonNaN_i);
        end
    end
data_first (:,:,timestep_start:timestep_start+size(data_second     ,3)-1)  =data_second (:,:,1:size(data_second     ,3));
else
%4-D data
for  forthLVL_current=1:size(data_first,3)
data_temp_first=zeros(size(data_first,1),size(data_first,2),size(data_first,4));
for jjjj=1:size(data_first,4)
data_temp_first(:,:,jjjj)=data_temp_first(:,:,jjjj)+data_first(:,:,forthLVL_current,jjjj);
end
data_temp_second=zeros(size(data_second,1),size(data_second,2),size(data_second,4));
for jjjj=1:size(data_second,4)
data_temp_second(:,:,jjjj)=data_temp_second(:,:,jjjj)+data_second(:,:,forthLVL_current,jjjj);
end
data_temp_first(:,:,timestep_start:timestep_start+size(data_temp_second     ,3)-1)=data_temp_second(:,:,1:size(data_temp_second     ,3));
data_first_combined(:,:,forthLVL_current,:)=data_temp_first;clear data_temp_first
end
data_first=data_first_combined;clear data_first_combined
end
end
end
data=data_first;clear data_*
date_time=date_time_first; clear date_time_*
else
data=data_1;clear data_*
date_time=date_time_1; clear date_time_*
end


end

