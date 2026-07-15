% Standard Entropy Calculator for Air
clear
close all
clc

%% Given / Reference Information

Rbar = 8.314462618; % [J/mol]
M = 28.95866;          % Molar mass of air [kg/kmol]
R = Rbar/M;         % Gas constant for air [kJ/kg-K]

% Temperature range
TK = 298:1300; %273:1800;   % Temperature K
t = TK./1000;    % Shomate temperature variable, t = T/1000

% Reference value from Lemmon et al. J. Phys. Chem. Ref. Data, Vol. 29, No. 3, 2000
TK_ref = 298.15;     % Reference temperature K
t_ref = TK_ref/1000; % Reference Shomate temperature variable

s_ref = 194.0/M;     % s^o at 298.15 K [J/mol] -> [kJ/kg-K]
h_ref = 8649.34/M;   % h^o at 298.15 K [J/mol] -> [kJ/kg-K]
u_ref = 6170.38/M;   % u^o at 298.15 K [J/mol] -> [kJ/kg-K]

%% Shomate Coefficients for Air

% Shomate form:
% C_{P0} = A + B*t + C*t^2 + D*t^3 + E/t^2
% where C_{P0} is in [J/mol-K] and t = TK/1000
% Dividing constants by M gives C_{P0} in [kJ/kg-K]
A = 28.11/M;
B = 1.967/M;
C = 4.802/M;
D = -1.966/M;
E = 0/M;

%% Define C_{P0} and Standard Entropy Function

% Standard-state constant-pressure specific heat
cp0 = @(t) A + B.*t + C.*t.^2 + D.*t.^3 + E./t.^2;
cp0 = @(TK) A + B.*(TK./1000) + C.*(TK./1000).^2 + D.*(TK./1000).^3 + E./(TK./1000).^2;

% Entropy antiderivative:
% integral(C_{P0}/TK) dTK
% Since t = TK/1000, this becomes:
% integral(C_{P0}/t) dt
s_int = @(t) A.*log(t) + B.*t + (C/2).*t.^2 + (D/3).*t.^3 - (E/2)./t.^2;

% Standard entropy function:
s0 = @(TK) s_ref + s_int(TK./1000) - s_int(t_ref);

%% Enthalpy Function
% dh = C_{P0} dTK -> h_2 - h_1 = int(C_{P0}, dTK)
% Since t = TK/1000 and dTK = 1000 dt:
% integral(C_{P0} dTK) = 1000*integral(C_{P0} dt)
h_int = @(t) 1000.*( A*t + (B/2).*t.^2 + (C/3).*t.^3 + (D/4).*t.^4 - E./t );

% Mass-specific standard enthalpy [kJ/kg]
h0 = @(TK) h_ref + h_int(TK./1000) - h_int(t_ref);

%% Internal Energy Function
% For an ideal gas:
% h0 = u0 + Rbar*T
% Therefore:
% u0 = h0 - Rbar*T

% Mass-specific standard internal energy [kJ/kg]
u0 = @(TK) h0(TK) - R.*TK;

%% Calculate Properties Over Temperature Range
% 
% cp0_air = cp0(TK);       % [kJ/kg-K]
[cp0_T,cp0_vals,~,~] = dT_Optimizer(298,1300,cp0,0.2);
% h0_air = h0(TK);         % [kJ/kg]
[h0_T,h0_vals,~,~] = dT_Optimizer(298,1300,h0,0.2);
% u0_air = u0(TK);         % [kJ/kg]
[u0_T,u0_vals,~,~] = dT_Optimizer(298,1300,u0,0.2);
% s0_air = s0(TK);         % [kJ/kg-K]
[s0_T,s0_vals,~,~] = dT_Optimizer(298,1300,s0,0.2);

%% Create Output Table
% air_table = table( ...
%     TK', ...
%     round(cp0_air',4), ...
%     round(h0_air',2), ...
%     round(u0_air',2), ...
%     round(s0_air',5), ...
%     'VariableNames', {'T_K','cp0_kJ_kgK','h0_kJ_kg','u0_kJ_kg','s0_kJ_kgK'} );
% disp(air_table)
% 

cp0_table = table( cp0_T',round(cp0_vals',4),'VariableNames',{'T [K]','cp0 [kJ/kg-K]'});