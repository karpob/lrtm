function [tcme,tcmp]= get_adams_data(oblateness_factor,data_directory,suffix)

%oblateness factor==ratio of polar to equatorial radius.
%data directory==relative path name where the data is stored
%suffix==filename.* where * is the suffix ie (mixrat.jup_test suffix=jup_test)

name=strcat(data_directory,'/mixrat.',suffix);
file_handle=fopen(name,'r');
oblateness_factor
%skip over the header
buffer_line1=fgetl(file_handle);
buffer_line2=fgetl(file_handle);
buffer_line3=fgetl(file_handle);
buffer_line4=fgetl(file_handle);
buffer_line5=fgetl(file_handle);
buffer_line6=fgetl(file_handle);
buffer_line7=fgetl(file_handle);

%read the data file 7c used to get rid of garbage fortran output

raw_data=fscanf(file_handle,'%f %7c %f %f %f %f %f %f %f %f %f %f %f %f %f ',[21 inf]);
raw_data=raw_data';
fclose(file_handle);

%filter out garbage Temperature output
data=[raw_data(:,1),zeros(length(raw_data),1),raw_data(:,9:21)];
z=data(:,1);
%pressure=data(:,3);
XH2O=2.0*data(:,4);
DH2O=data(:,5);
XNH3=data(:,6);
DNH3=data(:,7);
XH2S=data(:,8);
DH2S=data(:,9);
XCH4=data(:,10);
DCH4=data(:,11);
XAr=data(:,12);
DAr=data(:,13);
XPH3=data(:,14);
DPH3=data(:,15);

%Now get pressure and Temperature from initial model input file (need to
%for precision, and garbage fortran output in Temperature.
name2=strcat(data_directory,'/model.',suffix);
file_handle=fopen(name2,'r');
%Skip header
buffer_line1=fgetl(file_handle);
buffer_line2=fgetl(file_handle);
%read file
raw_data2=fscanf(file_handle,'%f %f %f %f %f %f %f %f %f %f',[10 inf]);
raw_data2=raw_data2';
fclose(file_handle);

temperature=raw_data2(:,3);
pressure=raw_data2(:,4);
XH2=raw_data2(:,5);
XHe=raw_data2(:,6);

zpp=oblateness_factor*z;

%flag to indicate "1 bar" pressure level has been found
located=0;
for ii=1:length(z)
    if(located==0)
        if((pressure(ii)-1.0)< 0.01)
            ze_offset=z(ii);
            zp_offset=zpp(ii);
            pressure(ii)-1.0
            located=1;        
        end 
    end
end

for ii=1:length(z)
    ze(ii)=z(ii)-ze_offset;
    zp(ii)=zpp(ii)-zp_offset;
end

blank_col=zeros(length(raw_data),1);
tcme=[pressure,temperature,ze',XH2,XHe,XH2S,XNH3,XH2O,XCH4,XPH3,blank_col,blank_col,DH2S,DNH3,DH2S,DNH3,DH2O,DCH4,DPH3,blank_col,blank_col,blank_col,blank_col,blank_col];
tcmp=[pressure,temperature,zp',XH2,XHe,XH2S,XNH3,XH2O,XCH4,XPH3,blank_col,blank_col,DH2S,DNH3,DH2S,DNH3,DH2O,DCH4,DPH3,blank_col,blank_col,blank_col,blank_col,blank_col];
size(tcme)
size(tcmp)
%tcme=tcme(75:486,:)
%tcmp=tcmp(75:486,:)