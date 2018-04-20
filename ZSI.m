function [rho, mu] = ZSI(TempC,P_Pa_gage)
%--------------------------------------------------------------------------
%   [ Rho, mu] = ZSI(TempC,P_Pa)
%   THIS PROGRAM CONTAINS EQUATIONS FOR CALCULATING THE DENSITY AND VISCOSITY 
%   See Zagarola (1996) Thesis, p294
%   -INPUT Pressure in Pascals gage pressure
%   -INPUT Temperature in degC 
%   -Input BOTH as vectors -OR- both as scalars
%--------------------------------------------------------------------------
oneatm = 101325.01;

TempK = TempC+273.15;
P_Pa_abs = P_Pa_gage + oneatm;
P_atm  = P_Pa_abs ./ oneatm;    %Convert to atmospheres
%
if size(TempC)~=size(P_Pa_abs)
    warning('Input dimensions do not agree!')
end
[l, w] = size(TempC);
%
Z1 = (-9.5378*10^-3) + ( 5.1986*10^-5).*TempK + (-7.0621*10^-8 ).*TempK.^2;
Z2 = ( 3.1753*10^-5) + (-1.7155*10^-7).*TempK + ( 2.4630*10^-10).*TempK.^2;
Z3 = ( 6.3764*10^-7) + (-6.4678*10^-9).*TempK + ( 2.1880*10^-11).*TempK.^2 + (-2.4691*10^-14).*TempK.^3;
%
for i = 1:l
    for j =1:w
    Z(i,j) = 1.0 + Z1(i,j)*(P_atm(i,j) - 1) + Z2(i)*(P_atm(i,j) - 1)^2 + Z3(i,j)*(P_atm(i,j) - 1)^3;
    end
end
%
for i = 1:l
    for j =1:w
    rho(i,j) = P_Pa_abs(i,j)/(TempK(i,j)*Z(i,j)*287.1);
    if (rho(i,j) ==0)
        rho(i,j) =1;
    end
    end
end
%
mu_0 = (1.458*10^-6)*((TempK.^1.5)./(110.4 + TempK));
mu_1 = (0) + (1.021*10^-8)*rho + (5.969*10^-11)*rho.^2 ;
%
for i = 1:l
    for j =1:w
    mu(i,j)   = mu_0(i,j) + mu_1(i,j);
    end
end
%
