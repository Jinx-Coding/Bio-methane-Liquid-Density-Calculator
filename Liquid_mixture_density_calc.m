
% ----------------  Liquid Mixture Density Calculator ----------------- %
% Giovanni Correra - 09/2024 %

clc
close all
clear variables

% ------------------------------ Data ---------------------------------- %

T = -146.3 + 273.15; % Temperature in (K) %
P = ( 3 + 1.01325); % Pressure in (bar) %

x = [0.9984, 0.0016, 0.00]; % Molar fractions of CH4, N2, CO2 %

Tc = [190.4, 126.2, 304.1]; % Critical temperatures in (K) %
Pc = [46.0, 33.9, 73.8]; % Critical pressures in (bar) %
om = [0.011, 0.039, 0.239]; % Omega coefficient (-) %
MW = [16.043, 28.013, 44.010]; % Molar masses in (g/mol) %

% -------------------------- Data check -------------------------------- %

if sum(x)~=1
    fprintf('ERROR : feed molar fractions do not add up to 1\n')
    return
end

% -------------------------- Main script ------------------------------- %

[Z,~,~,~,~] = SRK(T,P,Tc,Pc,om,x,2);
[rho,~] = density(Z,T,P,x,MW);

% ------------------------ Post - Processing --------------------------- %

fprintf('Liquid density = %.3f [kg/m3]\n',rho)
fprintf('\n')
fprintf('T = %.3f [°C]    ',T-273.15)
fprintf('P, abs = %.5f [bar]    ', P)
fprintf('P, rel = %.5f [bar]\n', P-1.01325)
fprintf('\n')
fprintf('Liquid composition [mol/mol]\n')
fprintf('x, CH4 = %.4f %%\n',x(1)*100)
fprintf('x, N2  = %.4f  %%\n',x(2)*100)
fprintf('x, CO2 = %.4f  %%\n',x(3)*100)
fprintf('\n')

% --------------------------- Functions -------------------------------- %

function [Z,AS,BS,A,B] = SRK(T,P,Tc,Pc,om,y,phase)

% SRK solving function with basic mixture laws %
% a(i,j) = (a(i)*a(j))^0.5 %
% b(i,j) = (b(i)+b(j))/2 %

% Phase = 1 (vapour), Phase = 2 (liquid) %

R = 8.3145;
RT = R*T;
RTc = R*Tc;

S = 0.48 + 1.574*om - 0.176*om.^2;
k = (1 + S.*(1-sqrt(T./Tc))).^2;
a = (0.42748*k.*RTc.^2)./Pc;
b = 0.08664*RTc./Pc;
AS = a.*P/(RT)^2;
BS = b.*P/(RT);
aM = sqrt(a'*a);
bM = zeros(length(Tc),length(Tc));
for i = 1:length(Tc)
    for j = 1:length(Tc)
        bM(i,j) = (b(i)+b(j))/2;
    end
end
am = y*aM*y';
bm = y*bM*y';
A = am*P/(RT)^2;
B = bm*P/(RT);
alfa = -1;
beta = A-B-B^2;
gamma = -A*B;

% Analytic solution %

p = beta - (alfa^2)/3;
q = 2*(alfa^3)/27 - alfa*beta/3 + gamma;
q2 = q/2;
a3 = alfa/3;
D = (q^2) / 4 + p^3 / 27;

if D>0
    Z1 = nthroot((-q2+sqrt(D)),3) + nthroot((-q2-sqrt(D)),3) - a3;
    Z = [Z1,Z1,Z1];
elseif D == 0
    Z1 = -2*nthroot(q2,3) - a3;
    Z2 = nthroot(q2,3) - a3;
    Z = [Z1,Z2,Z2];
elseif D<0
    r = sqrt(-p^3 / 27);
    teta = acos(-q2*sqrt(-27/p^3));
    Z1 = 2*nthroot(r,3)*cos(teta/3)-a3;
    Z2 = 2*nthroot(r,3)*cos((2*pi+teta)/3)-a3;
    Z3 = 2*nthroot(r,3)*cos((4*pi+teta)/3)-a3;
    Z = [Z1,Z2,Z3];
end

Z = sort(Z);

if phase == 1
    Z = max(Z);
elseif phase == 2
    Z = min(Z);
end

end

function [rho,v] = density(Z,T,P,x,MW)

v = Z*T*8.3451/(P*1e5); % Molar volume in (m3/mol) %
MW_tot = sum(x.*MW); % Molecular weight of mixture (g/mol) %
rho = (1/v)*MW_tot/1e3; % Mixture density (kg/m3) %

end

