% QBLADE file CL/CD file cleaner %
% Batch process Aerodyn formatted Lift and Drag curves %
% from Qblade/XFLR 5

clear all
close all

st  = 200000; %Start Reynolds number
en  = 4000000; %End Reynolds number
inc = 100000; %Increment of Reynolds number
% 
Re = st:inc:en;
Re = Re./1E6;
% Re = 3.230;

%INPUT file info%
folder = 'Tripped Foils\N63214\';
fname1 = 'N63214_T1_Re'; %Beginning of file name
fname2 = '_M0.00_N9.0_360_M';  %End of file name, leaving out Reynolds #
% fname2 = '';

% Middle points are ones actually solved for in XFLR5, other points are
% extrapolated to make curve fit 360degrees
aoarange = [-180:1:-11 -10:0.5:30 31:1:180];
%   aoarange = -10:0.5:25;

%OUTPUT file info%
wfolder = folder;
idx = strfind(fname1,'_T1'); %Strip off airfoil name
wfname = fname1(1:idx-1);
%Variable Declaration
% Cldum = zeros(size(aoarange)); Cddum = Cldum; %dummy variables
% Cl = zeros(numel(Cldum),numel(Re)); Cd = Cl;





for ii = 1:numel(Re)
    %Read in all data from XFOIL style file
    data = dlmread([folder fname1 num2str(Re(ii),'%2.3f') fname2 '.dat'],'',14,0);
    data = data(1:end-1,:); %Disclude last repeated 180deg point

    Cldum = interp1(data(:,1),data(:,2),aoarange);
    Cddum = interp1(data(:,1),data(:,3),aoarange);

    Cl(:,ii) = Cldum';
    Cd(:,ii) = Cddum';

    clear data Cddum Cldum
end

  hdr = {'First column is AoA (deg) and first row is Reynolds Number'};
  
  Clw = [00  Re.*1E6 ; aoarange' Cl];
  Cdw = [00  Re.*1E6 ; aoarange' Cd];
  
   %Write Lift Data
   fmt = repmat('%s\t', 1, length(hdr));
   fmt(end-1:end+2) = '\r\n';
   fid = fopen([wfolder wfname '_CL.txt'],'w');      
   fprintf(fid, fmt, hdr{:});
   fclose(fid); 
   dlmwrite([wfolder wfname '_CL.txt'],Clw,'-append','delimiter','\t','precision','%0.6f'); 
    
   %Write Drag Data
%    fmt = repmat('%s\t', 1, length(hdr));
%    fmt(end-1:end+2) = '\r\n';
%    fid = fopen([wfolder wfname '_CD.txt'],'w');      
%    fprintf(fid, fmt, hdr{:});
%    fclose(fid); 
%    dlmwrite([wfolder wfname '_CD.txt'],Cdw,'-append','delimiter','\t','precision','%0.6f'); 
