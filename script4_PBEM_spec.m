%%% Run PBEM Specify Inputs %%%
% More advanced plotting options%

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
%% Run BEM Code %% 
    nf = numel(Ufs);
for m = 1:nf
    disp(['Running: ' num2str(Ufs(m)) ' m/s (' num2str(m) ' of ' num2str(nf) ')'])
    % Assume tunnel conditions constant for run %
    tic;
    rotorfile = [rloc rname];
[ bemd(m), bld(m)] = pbem_3(rotorfile, pitch , nb, T, P, Ufs(m), TSR );
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

for k = 1:numel(Ufs)
       % Power Coefficient %
        fmain = figure(1);
        subplot(2,3,1)
        plot([bemd(k).TSR],[bemd(k).Cp],'x-',...
                'Linewidth',2);
            ylabel('C_P')
            set(gca,'Fontunits','points','Fontsize',fsize)
            ylim([0 0.6])
            xlim(xlimits)
            grid on

        subplot(2,3,4)
        plot([bemd(k).TSR],[bemd(k).Power],'x-',...
            'Linewidth',2);
            xlabel('TSR')
            ylabel('Power (W)')
            xlim(xlimits)
            set(gca,'Fontunits','points','Fontsize',fsize)
            grid on

   % Thrust Coefficient %
       subplot(2,3,2)
       plot([bemd(k).TSR],[bemd(k).Ct],'x-',...
           'Linewidth',2);
            ylabel('C_T')
            set(gca,'Fontunits','points','Fontsize',fsize)
            ylim([0 1.0])
            xlim(xlimits)
            grid on

        subplot(2,3,5)
        plot([bemd(k).TSR],[bemd(k).fx],'x-',...
            'Linewidth',2);
            xlabel('TSR')
            ylabel('Thrust (N)') 
            xlim(xlimits)
            set(gca,'Fontunits','points','Fontsize',fsize)
            grid on

   % Torque %
        subplot(2,3,6)
        a6(m) = plot([bemd(k).TSR], [bemd(k).Power]./([bemd(k).speed].*pi./30),...
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

    X = bld(k).abr ./ max(bld(k).abr);
        frad = figure(2);
        %-- AoA --%
        subplot(2,4,1)
       b1(m) = plot(X,bld(k).alpha(:,m));
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
        b2(m) = plot(X,bld(k).Cl(:,m));
             if m == 1
                hold on
                ylabel('C_L') 
                xlim([0 1])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
             end 
        %-- Cd --%
        subplot(2,4,6)
        b3(m) = plot(X,bld(k).Cd(:,m));
             if m == 1
                hold on
                ylabel('C_D') 
                xlim([0 1])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
             end 
        %-- Re_c --%
        subplot(2,4,5)
        b4(m) = plot(X,bld(k).Rec(:,m));
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
        b5(m) = plot(X,bld(k).a(:,m));
             if m == 1
                hold on
                ylabel('a')
                xlim([0 1])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
             end 
        %-- a' --%
        subplot(2,4,7)
        b6(m) = plot(X,bld(k).ap(:,m));
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
        b7(m) = plot(X,bld(k).Pn(:,m));
             if m == 1
                hold on
                ylabel('P_n')
                xlim([0 1])
                set(gca,'Fontunits','points','Fontsize',fsize)
                grid on
             end 
        %-- Pt --%
        subplot(2,4,8)
        b8(m) = plot(X,bld(k).Pt(:,m));
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
    legentry{m} = ['TSR-' num2str(TSR(m),'%2.0f') ...
                '_U' num2str(Ufs(k),'%2.0f')];
        
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
        uident = ['_U' num2str(Ufs(k),'%2.0f')];
    figure(fmain)
        savefig([cpath wfolder 'PBEM-Plots_' ident uident])
    figure(frad)
        savefig([cpath wfolder 'PBEM-Rad-Plots_' ident uident])
    save([cpath wfolder 'PBEM_Results_' ident uident],'bemd')
    save([cpath wfolder 'PBEM_Results-blade_' ident uident],'bld')
    end

   close all
end