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
        pitch   = []; %Collective rotor pitch, defined as (+) into the flow
        nb      = []; %Number of blades
        U       = []; %Free-stream velocity
        T       = []; %Dry-Air temperature (degrees Celsius)
        p       = []; %Static pressure (Pascals Gage) i.e. 0 = atmospheric
        TSR     = []; %Tip Speed Ratio, can be single value or array
        rotor_folder = '';
        %plot options struct%
        popt = struct('lw',1.5,'fsize',15);

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
        %Function to interactively set run conditions%
        function rt = runcon(rt)
            disp('Set the run conditions')
            rt.pitch = input('Rotor Pitch: ');
            rt.nb    = input('Number of blades: ');
            rt.T     = input('Tunnel Temperature (degrees C): ');
            rt.p     = input('Tunnel Static Pressure (Pa, gage): ');
            rt.U     = input('Free-stream Velocity (m/s): ');
            rt.TSR   = input('Tip Speed Ratio(s): ');
        end
        %fcn for graphical plot of rotor geometry
        function hp = rotor_plot(rt)
           hp = figure;
            scrsz = get(0,'screensize');
            set(hp,'Position',[5 50+scrsz(2) 0.8*scrsz(3) 0.85*scrsz(4)]) 
            
           subplot(2,3,[1 4])
            plot(rt.radius,rt.chord,'b-','LineWidth',rt.popt.lw)
            xlabel('Radius (m)')
            ylabel('Chord (m)')
            set(gca,'Fontunits','points','Fontsize',rt.popt.fsize)
            xlim([0 max(rt.radius)])
            
           subplot(2,3,[2 5])
            plot(rt.radius,rt.twist,'r-','LineWidth',rt.popt.lw)
            xlabel('Radius (m)')
            ylabel('Twist (^\circ)')
            set(gca,'Fontunits','points','Fontsize',rt.popt.fsize)
            xlim([0 max(rt.radius)])
            
           subplot(2,3,[3 6])
            tpos = get(gca,'OuterPosition');
            colnames = {'R (m)','c (m)','Twist (deg)','Airfoil'};
            numdat = num2cell([rt.radius rt.chord rt.twist]);
            cfm = {'short','short','short','char'};
            hui = uitable(hp,'Data',[numdat rt.foil],...
                'ColumnName',colnames,'ColumnFormat',cfm);
            set(hui,'Units','normalized','Position',tpos,...
                'ColumnWidth',{'auto','auto','auto',100});
 
        end
        %fcn for interpolating rotor geometry
         
        %fcn for loading experimental run conditions, lots options
        
        %Run single power curve (at 1 or several TSR values)%
        %Re-format for the current version of pbem
        function [bemd, bld] = run_pbem(rt)
            %bemd is a global-results struct.
            %bld is the blade-level results for each run in bemd
            rotor = [rt.radius rt.chord rt.twist];
            [bemd, bld ] = pbem_solver(rotor, rt.foil, rt.pitch ,...
                rt.nb, rt.T, rt.p, rt.U, rt.TSR );
        end
        %Simple Plotting Function%
        function hs = plot_s(rt,bemd)
            hs = 'none';
            
        end
        
        %plotting results, global
        function hp = plot_g(rt,bemd)
            tsr = [bemd.TSR]; Cp =  [bemd.Cp]; Ct =  [bemd.Ct];
            torq = [bemd.torq]; fx = [bemd.fx]; Power = [bemd.Power];
           %--- Power Coefficient ---%
            hp = figure;
            scrsz = get(0,'screensize');
            set(hp,'Position',[5 50+scrsz(2) 0.8*scrsz(3) 0.85*scrsz(4)])
            set(0,'DefaultTextInterpreter','Latex')
            subplot(2,3,1)
            plot(tsr,Cp,'x-','Linewidth',rt.popt.lw);
                ylabel('$C_P$')
                set(gca,'Fontunits','points','Fontsize',rt.popt.fsize)
                ylim([0 0.6])
                grid on

            subplot(2,3,4)
            plot(tsr,Power,'x-','Linewidth',rt.popt.lw);
                xlabel('$\lambda$')
                ylabel('Power (W)')
                set(gca,'Fontunits','points','Fontsize',rt.popt.fsize)
                grid on

       %--- Thrust Coefficient ---%
           subplot(2,3,2)
           plot(tsr,Ct,'x-','Linewidth',rt.popt.lw);
                ylabel('$C_t$')
                set(gca,'Fontunits','points','Fontsize',rt.popt.fsize)
                ylim([0 1.0])
                grid on

            subplot(2,3,5)
            plot(tsr,fx,'x-','Linewidth',rt.popt.lw);
                xlabel('$\lambda$')
                ylabel('Thrust (N)') 
                set(gca,'Fontunits','points','Fontsize',rt.popt.fsize)
                grid on

       %--- Torque ---%
            subplot(2,3,6)
            plot(tsr,torq,'x-','LineWidth',rt.popt.lw);
                grid on
                hold on
                xlabel('$\lambda$')
                ylabel('Torque (Nm)')
                set(gca,'FontUnits','points','FontSize',rt.popt.fsize)
        %--- Legend ---%
            subplot(2,3,3)
            aleg = plot(1,1,'x-','Linewidth',2);
            [~,nme,~] = fileparts(rt.name);
            legentry = [nme '; U = ' num2str(rt.U,'%2.0f')...
            ', T = ' num2str(rt.U,'%2.0f') ', P = ' num2str(rt.p,'%2.2E')];

            legh = legend(aleg, legentry,'Location','East');
            set(legh,'Fontunits','points','Fontsize',rt.popt.fsize-6,...
                'Interpreter','none') 
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

