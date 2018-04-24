classdef pbem_cls
    %pbemcls is a Blade-Element-Momentum class specifically developed
    %for analyzing the performance of wind turbine models operating in a
    %pressurized environment
    %   Detailed explanation goes here
    
    properties
        %These define rotor geometry for the object
        name    = {};
        radius  = [];
        chord   = [];
        twist   = [];
        foil    = [];
        %These define run conditions%
        pitch   = [];
        num_bld = [];
        Ufs     = [];
        T       = [];
        P       = [];
        TSR     = [];
        rotor_folder = '';
        %plot options struct%
        popt = struct('lw',1.5,'fsize',15);
%         lw = 1.5; %LineWidth, radial plots
%         fsize = 15; %Font Size
    end
    
    methods
        %Object builder for the pbem_cls class%
        %List and Load Rotor Geometry%
        function rt = pbem_cls
            %-- Location of BEM Rotor Models --%
            [cfold, ~, ~] = fileparts(mfilename('fullpath'));
            rt.rotor_folder = [cfold '\BEM_Models\'];
            avail_rot = dir([rt.rotor_folder '*.txt']);
            names = {avail_rot.name}';
            index = (1:1:numel(names))';
            t = table(index,names);
            disp('Available Rotors: ')
            disp(t);
            x = input('Index of rotor:');
            disp(['Rotor chosen is: ' names{x}])
            rt.name = names{x};
            %-Load rotor geometry-%
            fid = fopen([rt.rotor_folder rt.name],'r');
            rotorgem = textscan(fid,'%f %f %f %s','Delimiter',' ','MultipleDelimsAsOne',1,'Headerlines',1);
            fclose(fid);
            rt.radius = rotorgem{1};
            rt.chord  = rotorgem{2};
            rt.twist  = rotorgem{3};
            rt.foil   = rotorgem{4}; 
        end
        %fcn for interpolating rotor geometry
        
        %fcn for graphical plot of rotor geometry
        
        %fcn for loading experimental run conditions, lots options
        
        %Run single power curve (at 1 or several TSR values)%
        %Re-format for the current version of pbem
        function [bemd, bld] = run_pbem(rt, pitch, nb, T, P, U, TSR)
            %bemd is a global-results struct.
            %bld is the blade-level results for each run in bemd
            rotor = [rt.radius rt.chord rt.twist];
            [bemd, bld ] = pbem_solver(rotor, rt.foil, pitch , nb, T, P, U, TSR );
        end
        %Simple Plotting Function%
        function hs = plot_s(bemd)
            
            
        end
        
        %plotting results, global
        function hp = plot_g(bemd)
            tsr = [bemd.TSR]; Cp =  [bemd.Cp]; Ct =  [bemd.Ct];
            torq = [bemd.torq]; fx = [bemd.fx]; Power = [bemd.Power];
           %--- Power Coefficient ---%
            hp = figure(1);
            set(0,'DefaultTextInterpreter','Latex')
            subplot(2,3,1)
            plot(tsr,Cp,'x-','Linewidth',rt.lw);
                ylabel('$C_P$')
                set(gca,'Fontunits','points','Fontsize',rt.fsize)
                ylim([0 0.6])
                xlim(xlimits)
                grid on

            subplot(2,3,4)
            plot(tsr,Power,'x-','Linewidth',rt.lw);
                xlabel('$\lambda$')
                ylabel('Power (W)')
                xlim(xlimits)
                set(gca,'Fontunits','points','Fontsize',rt.fsize)
                grid on

       %--- Thrust Coefficient ---%
           subplot(2,3,2)
           plot(tsr,Ct,'x-','Linewidth',rt.lw);
                ylabel('$C_t$')
                set(gca,'Fontunits','points','Fontsize',rt.lw)
                ylim([0 1.0])
                xlim(xlimits)
                grid on

            subplot(2,3,5)
            plot(tsr,fx,'x-','Linewidth',rt.lw);
                xlabel('$\lambda$')
                ylabel('Thrust (N)') 
                xlim(xlimits)
                set(gca,'Fontunits','points','Fontsize',rt.fsize)
                grid on

       %--- Torque ---%
            subplot(2,3,6)
            plot(tsr,torq,'x-','LineWidth',rt.lw);
                grid on
                hold on
                xlabel('$\lambda$')
                ylabel('Torque (Nm)')
                xlim(xlimits)
                set(gca,'FontUnits','points','FontSize',rt.fsize)
        %--- Legend ---%
            subplot(2,3,3)
            aleg = plot(1,1,'x-','Linewidth',2);
            legentry{m} = ['$\lambda = $' num2str(bemd.TSR,'%2.0f') ...
            ' $U =$' num2str(bemd.U,'%2.0f')];
            legh = legend(aleg, legentry,'Location','East');
            set(legh,'Fontunits','points','Fontsize',rt.fsize-6,'Interpreter','none') 
            set(gca,'Visible','off')

            %Text Box%
    %             set(gcf,'Units','Normalized');
%             subplot(2,3,3)
%                 tb = get(gca,'Position');
%                 fstr = {['$Re_D =$ ' num2str(mean([bemd.ReD])./1E6,'%0.2f\n')...
%                          ' $\times 10^6$;' '  $\lambda =$ ' num2str(TSR)] ...
%                          ['$U_\infty =$ ' num2str(bemd.U)]...
%                          ['P = ' num2str(pstat) ' (PSI);'  ...
%                           '  T = ' num2str(T) ' degC'] ...
%                          ['Pitch = ' num2str(pitch)]};
%            annotation('textbox',[tb(1:2) 0.2 0.1],'String',fstr,...
%                 'FitBoxToText','on','Interpreter','tex');
            
            end
            %Plotting results, blade level

        
        %--EXPERIMENTAL--%
        %Fancy functions for changing out rotor geometry??
        %running pbem as a fortran or C program instead?
        %deflection code and analysis
    end
    
end

