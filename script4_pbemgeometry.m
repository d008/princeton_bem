%  Script for making PBEM V27 Rotor Geometry Files  %
% Outputs new rotor file in same directory as input %
%   M.Miller 1-18-18 ;]                             %

%Get the original rotor geometry file%
rfolder = 'BEM_Models';
rotorname = 'NREL_5MW_Model_BEM_fixRe.txt';
fid = fopen([rfolder rotorname],'r'); %Format input file
rotorgem = textscan(fid,'%f %f %f %s','Delimiter',' ','MultipleDelimsAsOne',1,'Headerlines',1);
fclose(fid);
rotor = [rotorgem{1} rotorgem{2} rotorgem{3}];
foil = rotorgem{4}; 

%% Options for rotor creation %%
rchop = 1; %Controls where last element is located in r/R
NE = 25  ; %Number of blade elements
% spacing = 'linear'; %'linear' or 'sine', ADD THIS FUNCTIONALITY LATER

% Discretize Rotor into Blade Elements %
    % Find end blade element %
    rend = rotor(end,1).* rchop;
    % Find length of blade elements
    db = (rend - rotor(1,1)) / NE;
    % create array of radial locations for blade elements  
    abr = linspace(rotor(1,1) + db/2, rend - db/2, NE)'; 
    % Interpolate array of chord elements
    abc = interp1(rotor(:,1),rotor(:,2),abr);
    % Interpolate twist for elements, in degrees
    abt = interp1(rotor(:,1),rotor(:,3),abr);
    % Find airfoils for each blade element
    Afidx = round(interp1(rotor(:,1),1:1:length(rotor(:,1)-1),abr));
    for mm = 1:numel(Afidx)
        %populate list with airfoils
       aft(mm) = foil(Afidx(mm)); 
    end 
    
 %% Write new rotor to file %%
 [~,fname,~] = fileparts(rotorname);
 wname = [rfolder fname '_PBEM-geom_NE-' num2str(NE) '.txt'];
    fileID = fopen(wname,'w');
    fprintf(fileID,'%s %s %s %s\n','R(m)	|',' c(m)	|',' Twist (deg)|',' Airfoil');
    
    for mm = 1:numel(Afidx)
      fprintf(fileID,'%4.6f \t %4.6f \t %4.6f \t',[abr(mm) abc(mm) abt(mm)]);  
      fprintf(fileID,'%s\n',aft{mm});
    end
    
    fclose(fileID);
    
    
    
    
    