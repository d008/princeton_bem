function [bemd, bld ] = pbem_solver(rotor, foil, pitch , nb, T, P, U, TSR )
%Blade Element Momentum code for HRTF
%INPUTS:
% rotor =  struct created with load_rotor
% nb = # of blades; T = Tunnel Temp (degC); 
% P = Ambient (Tunnel) Pressure (Pa) gage; U = Tunnel Velocity (m/s)
% TSR = Tip Speed Ratio (array input accepted);
%OUTPUTS:
% bemd = global properties in struct;
% bld = blade level results in struct array, run conditions are index values;
%%% Mark Miller 4-20-18 %%% ;]
%  

%Location of the airfoil data%
rfolder = [mfilename('fullpath') '\Airfoil_Data\Tripped Foils\'];
nmax = 1000;    % Maximum number of iterations allowed
tol = 1E-6;     % Convergence tolerance
relax = 0.5;    % Relaxation factor

     % Check to make sure input values are all the correct size  %
      if numel(T) + numel(P) + numel(U) ~= 3
          error('Input Conditons not equal size (T,P,U)!!')
      elseif numel(foil) ~= numel(rotor(:,1))
          disp(['Number of airfoils: ' num2str(numel(foil))])
          disp(['Number of radial locations: ' num2str(numel(rotor(:,1)))])
          error('Airfoils not assigned or too few radius locations given')
      end

% Determine tunnel Conditions %
      [rho, mu ] = ZSI(T,P);
      % Find local pitch angle in radians %
      theta = (rotor(:,3) + pitch).*pi./180;
      % Find rotor Solidity %
      sigma = rotor(:,2) * nb ./ (2*pi .* rotor(:,1));
      % Find rotor rotational frequency %
      omega = TSR .* U ./ rotor(end,1);
% Load airfoil data  
    for m = 1:numel(rotor(:,1))
        %Load data into 3-dimensional matrix, 
        % ( AoA , Re , Foil Type )
        dum = dlmread([rfolder char(foil(m)) '\' char(foil(m)) '_CL.txt'],'\t',1,0);
        Cl(:,:,m) = dum(2:end,2:end);
        Re(:,m) = dum(1,2:end);
        aoa(:,m) = dum(2:end,1).*pi./180; %Convert to radians
        Cd(:,:,m) = dlmread([rfolder char(foil(m)) '\' char(foil(m)) '_CD.txt'],'\t',2,1);
        clear dum
    end

%-------------------%    
%  Enter BEM Code   %
%-------------------%

%Loop through all input conditions%
for k = 1:numel(omega)
    disp(['Currently on TSR value ' num2str(k) ' of ' num2str(numel(omega))])
      %  Reserve variables  %
      a = zeros(length(rotor(:,1)),1); ap = a; alpha = a; phi = a;
      Clc = a; Cdc = a; Cn = a; Ct = a; Ulocal = a; Rec = a;
    for jj = 1:length(rotor(:,1))-1 %Loads at end of rotor are 0
        m = 1; %re-initialize while loop counter
        while m <= nmax + 1
            if m == nmax %Error message if maximum iteration count is reached
                disp('Maximum # iterations reached for airfoil at ')
                disp(['R = ' num2str(rotor(jj,1)) ' and foil: ' char(foil(jj))])
                break
            end
         % Find local blade velocity   
             Ulocal(jj) = ((U .* (1-a(jj))).^2 + ...
                 (omega(k) .* rotor(jj,1) .* (1 + ap(jj))).^2).^(0.5);  
             Rec(jj) = rho .* rotor(jj,2) .* Ulocal(jj) ./ mu;
%              Rec(jj) = 0.2E6;  %Forces a specific Re
         % Find inflow angle for each position  
             phi(jj) = atan((1 - a(jj)) .* U ./ ...
                    ((1 + ap(jj)) .* omega(k) .* rotor(jj,1) ));
         % Find local aoa
             alpha(jj) = phi(jj) - theta(jj);

         % Find closest Re for each section %
            dmin = abs(Re(:,jj) - repmat(Rec(jj)',[size(Re,1) 1])); 
            [~, idx] = min(dmin);
            clear dmin   
         % Interpolate between available AoA
             Clc(jj) = interp1(aoa(:,jj),Cl(:,idx,jj),alpha(jj));
             Cdc(jj) = interp1(aoa(:,jj),Cd(:,idx,jj),alpha(jj));
        %Transform to Normal and Tangential Directions
            Cn(jj) = Clc(jj)*cos(phi(jj)) + Cdc(jj)*sin(phi(jj));
            Ct(jj) = Clc(jj)*sin(phi(jj)) - Cdc(jj)*cos(phi(jj));
        % Prandtl's Tip Loss Factor, set F=1 to turn off
            if sin(phi(jj)) < 0.02 %Avoid sigularities in correction
                F = 1;
            else
                eff = nb/2*(rotor(end,1) - rotor(jj))/(rotor(jj)*sin(phi(jj)));
                F = 2*acos(exp(-eff)) / pi;
            end 
%             F = 1;
       % Re-calculate a and a' %
        %Use Glauert Correction
            if a(jj) > 0.2 %From Spera (1994), see Hansen CH 6
                K = 4.*F.*sin(phi(jj))^2. / (sigma(jj) * Cn(jj));
                anew = 0.5*(2 + K*(1-2*0.2) - sqrt((K*(1-2*0.2)+2)^2 + ...
                    4*(K*0.2^2 - 1)));
            else
             anew = 1. / (4*F*sin(phi(jj))^2 / (sigma(jj)*Cn(jj)) + 1);
            end
             apnew = 1. / (4*F*sin(phi(jj))*cos(phi(jj)) / (sigma(jj)*Ct(jj)) - 1 );
        %Relaxation constant set to 1 for first few iterations as in Qblade docs    
            if m < 2 
                arelax = 1;
            elseif m == 2
                %Use 3 point eqn to find center of initial oscillations
                %See Qblade documentation
                arelax = 1;
                anew = 0.25*anew + 0.5*a(jj) +0.25*aprev; 
            else
                arelax = relax;
            end
        % Calculate error for each foil %
            if abs(anew - a(jj)) < tol && abs(apnew - ap(jj)) < tol
                a(jj) =  anew;
                ap(jj) = apnew;
                break
            else
                aprev = a(jj); %Save previous iteration value
                a(jj) = arelax*anew +(1-arelax)*a(jj); 
                ap(jj) = arelax*apnew + (1-arelax)*ap(jj);
            end

        clear anew apnew pl pr idx Y K
        
        m = m + 1; %Iteration Counter  
        end
    end

        %Find normal and tangential forces per unit length
            Pn = 0.5*rho.*Ulocal'.^2.*rotor(:,2)'.*Cn';
            Pt = 0.5*rho.*Ulocal'.^2.*rotor(:,2)'.*Ct';
        %Calculate hub moment, linear variation btwn points
            A = zeros(length(rotor(:,1)) - 1,1); B = A; M = A;
        for m = 1:length(rotor(:,1))-1
            A(m) = (Pt(m+1) - Pt(m))/(rotor(m+1,1) - rotor(m,1));
            B(m) = (Pt(m)*rotor(m+1,1) - Pt(m+1)*rotor(m,1)) / ...
                (rotor(m+1,1) - rotor(m,1));
            M(m) = (1/3)*A(m)*(rotor(m+1,1)^3 - rotor(m,1)^3) + ...
                (1/2)*B(m)*(rotor(m+1,1)^2 - rotor(m,1)^2);
        end
        %Calculate Thrust force loading
            Tn = trapz(rotor(:,1),Pn);
            Tt = trapz(rotor(:,1),Pt);
            Mtot = nb * sum(M); %Total hub moment
            Ttot = nb * Tn;     %Total axial thrust
            
        %Calculate global properties:
        bemd.Ct(k) = Ttot ./ (0.5 .* rho .* U.^2 .* pi .* rotor(end,1).^2);
        bemd.Cp(k) = (Mtot * omega(k)) ./ (0.5*rho*U.^3*pi*rotor(end,1)^2);
        bemd.fx(k) = Ttot; bemd.M(k) = Mtot; bemd.ReD(k) = rho*U*rotor(end,1)*2./mu;
        bemd.Power(k) = Mtot * omega(k); 
        bemd.Retip(k) = rotor(end,2).*rho.*sqrt(U.^2 + (omega(k).*0.1).^2)./mu;
        %Blade level variables
        bld.Pt(:,k) = Pt';  bld.Pn(:,k) = Pn';       bld.Rec(:,k)   = Rec;
        bld.a(:,k)  = a;    bld.ap(:,k) = ap;   bld.Cl(:,k)    = Clc;
        bld.Cd(:,k) = Cdc;   bld.alpha(:,k) = alpha .* 180 ./ pi;
        bld.phi(:,k) = phi .* 180 ./ pi; 
        bld.Urtan(:,k) = omega(k).*rotor(:,1).*(1+ap); 
        bld.Urnorm(:,k) = U.*(1-a);
end
        bld.abr = rotor(:,1); bld.chord = rotor(:,2); bld.twist = rotor(:,3);
        bld.foil = foil;
        bemd.rho = rho; bemd.mu = mu; bemd.U = U; 
        bemd.TSR = omega .* rotor(end,1) ./ U ;   
        bemd.speed = omega .* 30 ./ pi;

end


