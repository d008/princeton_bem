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
        rotor_folder = '';
        %Simulation Inputs%
        runcon = struct([]);
        %Simulation Outputs%
        bemd = struct([]);
        bld  = struct([]);
        %plot options struct%
        popt = struct('lw',1.5,'fsize',15);

    end
    
    methods

        function rt = pbem_cls
            %Object builder for the pbem_cls class%
             %List and Load Rotor Geometry%
             %Location of BEM Rotor Models%
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
            %---Struct for solver inputs---%
            fn = {'pitch','nb','U','T','p','TSR'};
            fv = cell(1,numel(fn));
            fo = [fn;fv];
            rt.runcon = struct(fo{:});
            %pitch      %Collective rotor pitch, defined as (+) into the flow
            %nb         %Number of blades
            %U          %Free-stream velocity
            %T          %Dry-Air temperature (degrees Celsius)
            %p          %Static pressure (Pascals Gage) i.e. 0 = atmospheric
            %TSR        %Tip Speed Ratio, can be single value or array
            %---Stucts for solver output---%
           fn = {'Ct','Cp','ReD','Retip','TSR','fx','torq','Power',...
                 'rho','mu','speed','omega'};
           fv = cell(1,numel(fn));
           fo = [fn;fv];
           rt.bemd = struct(fo{:});
           fn = {'Pt','Pn','Rec','a','ap','alpha',...
                 'Cl','Cd','phi','Urtan','Urnorm'};
           fv = cell(1,numel(fn));
           fo = [fn;fv];
           rt.bld = struct(fo{:});
           
        end
        
        function rt = runsetup(rt)
        %Interactively set run conditions%
            %Pick the index for setting run conditions%
            if isempty(rt.runcon(1).pitch)
                ind = 1;
            else
                disp('Stored Runs: ')
                rt.disp_runs
                dum = size(rt.runcon,2) + 1;
                disp(['--> Or select index ' num2str(dum,'%d') ' for new run'])
                ind = input('Select column index: ');
            end
            disp('Set the run conditions')
            rt.runcon(ind).pitch = input('Rotor Pitch: ');
            rt.runcon(ind).nb    = input('Number of blades: ');
            rt.runcon(ind).T     = input('Tunnel Temperature (degrees C): ');
            rt.runcon(ind).p     = input('Tunnel Static Pressure (Pa, gage): ');
            rt.runcon(ind).U     = input('Free-stream Velocity (m/s): ');
            tdum   = input('Tip Speed Ratio(s): ');
            rt.runcon(ind).TSR(1:numel(tdum)) = tdum;
            %Display runs to user
            disp_runs(rt)
        end
        
        function hp = rotor_plot(rt)
        % Graphical plot of rotor geometry %
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
            ylabel('Twist ($^\circ$)')
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
        
        function disp_runs(rt,varargin)
        %Display the run conditions currently available by index
            ic = size(rt.runcon,2);
            for j = 1:ic
                ntsr = numel( rt.runcon(j).TSR(:) );
                disp(['---------- Index ' num2str(j,'%d') ' ----------'])
                if ntsr == 0
                    disp('No run conditions set for this index!')
                else
                    dispstr = ['U: %2.1f P: %2.1E T: %2.1f TSR: ' ...
                        repmat('%2.1f; ',1,ntsr) ' \n' ];
                    disparray = {rt.runcon(j).U, rt.runcon(j).p,...
                        rt.runcon(j).T, rt.runcon(j).TSR(1:ntsr)};
                    fprintf(dispstr,disparray{:});
                    %Display the estimated ReD and Rec values
                    disp('Estimated Run conditions: ')
                    [rho, mu] = ZSI(rt.runcon(j).T, rt.runcon(j).p);
                    D = rt.radius(end) *2; tipc = rt.chord(end);
                    ReD = rho .* D .* rt.runcon(j).U ./ mu;
                    Rec = rho .* tipc .* (rt.runcon(j).U).* ...
                        sqrt(1 + rt.runcon(j).TSR(1:ntsr).^2) ./ mu;
                    disp(['Re_D =  ' num2str(ReD,'%2.2E')])
                    dispstr2 = ['Re_c = ' repmat('%2.2E; ',1,ntsr) '\n'];
                    fprintf(dispstr2,Rec) 
                end
                disp('-----------------------------')
            end
        end
        
       
        function rt = run_pbem(rt)
        %Run single power curve (at 1 or several TSR values)% 
            if isempty(rt.runcon(1).pitch)
                disp('No run conditions set, use rt = rt.runsetup')
            else
                disp('Choose index for run: ')
                disp('Run results will be overwritten!')
                rt.disp_runs
                    
                indrun = input('Chosen Run Index: ');
            end
            rotor = [rt.radius rt.chord rt.twist];
            [bemdo, bldo] = pbem_solver(rotor, rt.foil, rt.runcon(indrun));
            rt.bemd(indrun) = bemdo;
            rt.bld(indrun)  = bldo;
            
        end
        %Simple Plotting Function%
        function hs = plot_s(rt,bemd)
            hs = 'none';
            
        end
        
        %plotting results by index, global
        function hp = plot_g(rt,ind)
            tsr = [rt.bemd(ind).TSR];  Cp    =  [rt.bemd(ind).Cp]; 
            Ct  = [rt.bemd(ind).Ct];   torq  = [rt.bemd(ind).torq];
            fx  = [rt.bemd(ind).fx];   Power = [rt.bemd(ind).Power];
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
            legentry = [nme '; U = ' num2str(rt.runcon(ind).U,'%2.0f')...
            ', T = ' num2str(rt.runcon(ind).U,'%2.0f') ...
            ', P = ' num2str(rt.runcon(ind).p,'%2.2E')];

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

