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
        %-- Location of BEM Rotor Models --%
        rotor_folder = [mfilename('fullpath') '\BEM_Models\'];
    end
    
    methods
        %Object builder for the pbem_cls class%
        %Name specifies the rotor name
        function rt = pbem_cls(name)
           rt.name = name;
        end
        %Function for loading rotor geometry from file%
        function rt = rotor_load(rt)
            fid = fopen([rt.rotor_folder rt.name '.txt'],'r');
            rotorgem = textscan(fid,'%f %f %f %s','Delimiter',' ','MultipleDelimsAsOne',1,'Headerlines',1);
            fclose(fid);
            rt.radius = rotorgem{1};
            rt.chord  = rotorgem{2};
            rt.twist  = rotorgem{3};
            rt.foil   = rotorgem{4}; 
        end
        %fcn for interpolating loaded rotor geometry
        
        %fcn for graphical plot of rotor geometry
        
        %fcn for loading experimental run conditions, lots options
        
        %Run single power curve (1+ TSR values)%
        %Re-format for the current version of pbem
        function [bemd, bld] = run_pbem(rt, pitch, nb, T, P, U, TSR)
            rotor = [rt.radius rt.chord rt.twist];
            [bemd, bld ] = pbem_solver(rotor, rt.foil, pitch , nb, T, P, U, TSR );
        end
        %plotting results, global
        
        %Plotting results, blade level
        
        
        %--EXPERIMENTAL--%
        %Fancy functions for changing out rotor geometry??
        %running pbem as a fortran or C program instead?
        %deflection code and analysis
    end
    
end

