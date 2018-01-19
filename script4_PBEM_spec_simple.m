%%% Run PBEM Specify Inputs %%%
clc
close all 
clear all

% Folder to write results %
wfolder = '\Results\';
% Save File Identifier %
ident = 'V27_Model_BEM_3xchord_PBEM-geom_NE-25';

dwrite = 0 ;         % Set to 1 to write data to file
pitch  = 5 ;         % Rotor Pitch (degrees)
pstat  = 1528 ;      % Tunnel Static Pressure (PSI)
T      = 20  ;       % Tunnel Temperature (degC)
Ufs    = 5   ;       % Tunnel Velocity (m/s), single or array input
TSR    = 3:1:8 ;     % Tip Speed Ratio, single or array input
nb     = 3;          % Number of Blades
%Rotor geometry file location
    %change this to location of rotor file
    [cpath , ~, ~] = fileparts(mfilename('fullpath'));
    rloc = [cpath '\BEM_Models\'];
    rname = 'V27_Model_BEM_3xchord_PBEM-geom_NE-25.txt';
    P       = pstat * 6894.75728;
    rotorspeed
%% Run BEM Code %% 
    nf = numel(Ufs);
for m = 1:nf
    disp([ num2str(Ufs(m)) ' rot/min (' num2str(m) ' of ' num2str(nf) ')'])
    % Assume tunnel conditions constant for run %
    tic;
    rotorfile = [rloc rname];
[ bemd(m), bld(m)] = pbem_3(rotorfile, pitch , nb, T, P, U, TSR );
    t(m) = toc;
end

%% Plotting Section %%
set(0,'DefaultTextInterpreter','tex')
disp(['Total Run Time: ' num2str(sum(t))])
fsize = 18; %Font Size
xlimits = [min([bemd.TSR]) - 0.5 max([bemd.TSR]) + 0.5];

       % Power Coefficient %
        fmain = figure(1);
        subplot(2,3,1)
        plot([bemd.TSR],[bemd.Cp],'k--','Linewidth',2);
            ylabel('C_P')
            set(gca,'Fontunits','points','Fontsize',fsize)
            ylim([0 0.6])
            xlim(xlimits)
            grid on
        subplot(2,3,4)
        plot([bemd.TSR],[bemd.Power],'k--','Linewidth',2);
            xlabel('TSR')
            ylabel('Power (W)')
            xlim(xlimits)
            set(gca,'Fontunits','points','Fontsize',fsize)
            grid on
   % Thrust Coefficient %
       subplot(2,3,2)
       plot([bemd.TSR],[bemd.Ct],'k--','Linewidth',2);
            ylabel('C_T')
            set(gca,'Fontunits','points','Fontsize',fsize)
            ylim([0 1.0])
            xlim(xlimits)
            grid on
        subplot(2,3,5)
        plot([bemd.TSR],[bemd.fx],'k--','Linewidth',2);
            xlabel('TSR')
            ylabel('Thrust (N)') 
            xlim(xlimits)
            set(gca,'Fontunits','points','Fontsize',fsize)
            grid on
   % Torque %
        subplot(2,3,6)
        plot([bemd.TSR],[bemd.Power]./([bemd.speed].*pi./30),'k--','LineWidth',2)
        grid on
        xlabel('TSR')
        ylabel('Torque (Nm)')
        xlim(xlimits)
        set(gca,'FontUnits','points','FontSize',fsize)

        scrsz = get(0,'screensize');
        set(gcf,'Position',[5 50+scrsz(2) 0.55*scrsz(3) 0.85*scrsz(4)]) 
        
   if dwrite == 1
      figure(fmain)
      savefig([wfolder 'PBEMspec-Plots_' ident])
   end  

   
 %% Data Saving %%
 
 % Write Data to File %
   if dwrite == 1
       save([wfolder 'PBEMspec_Results_' ident],'bemd')
       save([wfolder 'PBEMspec_Results-blade_' ident],'bld')
   end

    
    
    
    
     
