%%% Run PBEM Specify Inputs %%%
clc
close all 
clear all
% Folder to write results %
wfolder = 'C:\Users\mamil\Documents\GitHub\princeton_bem\Results\';
% Save File Identifier %
ident = 'NREL_5MW_BEM_test';

dwrite = 1 ;         % Set to 1 to write data to file
pitch  = 5 ;         % Rotor Pitch (degrees)
pstat  = 3100 ;      % Tunnel Static Pressure (PSI)
T      = 23  ;       % Tunnel Temperature (degC)
Ufs    = 3:1:8   ;   % Tunnel Velocity (m/s)
rotorspeed = 12.1;   % Rotor Speed (RPM)
%Rotor geometry file location
    %change this to location of rotor file
    rloc = '\BEM_Models\';
    rname = 'NREL_5MW_BEM_fixRe.txt';
    
    fid = fopen([rloc rname],'r'); %Format input file
    rotorgem = textscan(fid,'%f %f %f %s','Delimiter',' ','MultipleDelimsAsOne',1,'Headerlines',1);
    fclose(fid);
    rotor = [rotorgem{1} rotorgem{2} rotorgem{3}];
    foil = rotorgem{4}; 

    P       = pstat * 6894.75728;
    Rmax    = max(rotor(:,1));
    speed   = rotorspeed;
    TSR     = speed .* Rmax .*pi ./ (30 .* Ufs);
    
%% Run BEM Code %% 
for m = 1:numel(TSR)
    disp(['TSR: ' num2str(TSR(m)) ' (' num2str(m) ' of ' num2str(numel(TSR)) ')'])
    % Assume tunnel conditions constant for run %
    U = Ufs(m);
    tic;
[ bemd(m), bld(m)] = pbem_2(rotor, foil, T, P, speed, pitch , U , 3 );
    t(m) = toc;
end

%% Plotting %%
close all
dwrite = 1 ;         % Set to 1 to write data to file
lw = 1.5; %LineWidth, radial plots
set(0,'DefaultTextInterpreter','tex')
disp(['Total Run Time: ' num2str(sum(t))])
fsize = 15; %Font Size
xlimits = [min(TSR) - 0.5 max(TSR) + 0.5];


       % Power Coefficient %
        fmain = figure(1);
        subplot(2,3,1)
        plot([bemd.TSR],[bemd.Cp],'x-',...
                'Linewidth',2);
            ylabel('C_P')
            set(gca,'Fontunits','points','Fontsize',fsize)
            ylim([0 0.6])
            xlim(xlimits)
            grid on

        subplot(2,3,4)
        plot([bemd.TSR],[bemd.Power],'x-',...
            'Linewidth',2);
            xlabel('TSR')
            ylabel('Power (W)')
            xlim(xlimits)
            set(gca,'Fontunits','points','Fontsize',fsize)
            grid on

   % Thrust Coefficient %
       subplot(2,3,2)
       plot([bemd.TSR],[bemd.Ct],'x-',...
           'Linewidth',2);
            ylabel('C_T')
            set(gca,'Fontunits','points','Fontsize',fsize)
            ylim([0 1.0])
            xlim(xlimits)
            grid on

        subplot(2,3,5)
        plot([bemd.TSR],[bemd.fx],'x-',...
            'Linewidth',2);
            xlabel('TSR')
            ylabel('Thrust (N)') 
            xlim(xlimits)
            set(gca,'Fontunits','points','Fontsize',fsize)
            grid on

   % Torque %
        subplot(2,3,6)
        a6(m) = plot([bemd.TSR], [bemd.Power]./([bemd.speed].*pi./30),...
            'x-','LineWidth',2);
            grid on
            hold on
            xlabel('TSR')
            ylabel('Torque (Nm)')
            xlim(xlimits)
            set(gca,'FontUnits','points','FontSize',fsize)
    % Legend %
        subplot(2,3,3)
        aleg = plot(1,1,'x-','Linewidth',2);
        
%----- Blade-Level Data -----%        
for m = 1:numel(TSR)

    X = bld(m).abr ./ max(bld(m).abr);
        frad = figure(2);
        %-- AoA --%
        subplot(2,4,1)
       b1(m) = plot(X,bld(m).alpha);
            if m == 1
                hold on
                ylabel('\alpha (degrees)') 
                xlim([0 1])
                ylim([0 40])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
                grid minor
            end
        %-- Cl --%
        subplot(2,4,2)
        b2(m) = plot(X,bld(m).Cl);
             if m == 1
                hold on
                ylabel('C_L') 
                xlim([0 1])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
             end 
        %-- Cd --%
        subplot(2,4,6)
        b3(m) = plot(X,bld(m).Cd);
             if m == 1
                hold on
                ylabel('C_D') 
                xlim([0 1])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
             end 
        %-- Re_c --%
        subplot(2,4,5)
        b4(m) = plot(X,bld(m).Rec);
             if m == 1
                hold on
                ylabel('Re_c') 
                xlabel('r/R')
                xlim([0 1])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
             end
        %-- a --%
        subplot(2,4,3)
        b5(m) = plot(X,bld(m).a);
             if m == 1
                hold on
                ylabel('a')
                xlim([0 1])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
             end 
        %-- a' --%
        subplot(2,4,7)
        b6(m) = plot(X,bld(m).ap);
             if m == 1
                hold on
                ylabel('a''')
                xlabel('r/R')
                xlim([0 1])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
             end 
        %-- Pn --%
        subplot(2,4,4)
        b7(m) = plot(X,bld(m).Pn);
             if m == 1
                hold on
                ylabel('P_n')
                xlim([0 1])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
             end 
        %-- Pt --%
        subplot(2,4,8)
        b8(m) = plot(X,bld(m).Pt);
             if m == 1
                hold on
                ylabel('P_t')
                xlabel('r/R')
                xlim([0 1])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
             end 
     %Fix legend entries%
%         dum = char(TSR{m});
%         ind = strfind(dum,'_PBEM');
%         dum = dum(1:ind-1);
%         if length(dum) >= 25
%         dum = [dum(1:25) char(10) dum(26:end)];        
%         end
%         legentry{m} = dum;
    legentry{m} = [ident '_TSR-' num2str(TSR(m),'%2.0f') ];
        
    if m == numel(TSR)
        set(b1,'LineWidth',lw) 
        set(b2,'LineWidth',lw) 
        set(b3,'LineWidth',lw) 
        set(b4,'LineWidth',lw) 
        set(b5,'LineWidth',lw) 
        set(b6,'LineWidth',lw) 
        set(b7,'LineWidth',lw)
        set(b8,'LineWidth',lw)
        
        legh = legend(b1,legentry,'Location','NorthWest');
        set(legh,'Fontunits','points','Fontsize',fsize-6,'Interpreter','none') 
        
        figure(fmain)
        legh = legend(aleg,legentry,'Location','East');
        set(legh,'Fontunits','points','Fontsize',fsize-6,'Interpreter','none') 
        set(gca,'Visible','off')
        
        %Text Box%
%             set(gcf,'Units','Normalized');
        subplot(2,3,3)
            tb = get(gca,'Position');
            fstr = {['Re_D = ' num2str(mean([bemd.ReD])./1E6,'%0.2f\n')...
                     ' x 10^6;' '  \lambda = ' num2str(TSR)] ...
                     ['U_\infty = ' num2str(Ufs)]...
                     ['P = ' num2str(pstat) ' (PSI);'  ...
                      '  T = ' num2str(T) ' degC'] ...
                     ['Pitch = ' num2str(pitch)]};
       annotation('textbox',[tb(1:2) 0.2 0.1],'String',fstr,...
            'FitBoxToText','on','Interpreter','tex');
     end
     
       
end
 scrsz = get(0,'screensize');
 set(frad,'Position',[5 50+scrsz(2) 0.8*scrsz(3) 0.85*scrsz(4)])  
 set(fmain,'Position',[5 50+scrsz(2) 0.8*scrsz(3) 0.85*scrsz(4)]) 
  % Data Saving %
    if dwrite == 1
    figure(fmain)
        savefig([wfolder 'PBEMspec-Plots_' ident])
    figure(frad)
        savefig([wfolder 'PBEM-Rad-Plots_' ident])
    save([wfolder 'PBEMspec_Results_' ident],'bemd')
    save([wfolder 'PBEMspec_Results-blade_' ident],'bld')
   end
