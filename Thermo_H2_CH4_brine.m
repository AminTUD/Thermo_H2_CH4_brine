function [H2, CH4, H2O, Z] = Thermo_H2_CH4_brine(T, Salinity)
    
% This MATLAB function calculates thermodynamic properties and phase 
% equilibria for a system containing hydrogen (H2), methane (CH4), and brine
% (H2O with salinity). The code is based on the Peng-Robinson Equation of 
% State (PR EOS) and follows the methodology described by 
% Ziabakhsh-Ganji and Kooi (2012).
%
% References:
%
% 1. Ziabakhsh-Ganji, Z., & Kooi, H. (2012). An Equation of State for 
% thermodynamic equilibrium of gas mixtures and brines to allow simulation 
% of the effects of impurities in subsurface CO2 storage. International 
% Journal of Greenhouse Gas Control, 11, S21-S34.
%
% 2. Shabani, B., & Vilcáez, J. (2019). TOUGHREACT-CO2Bio–A new module to 
% simulate geological carbon storage under biotic conditions (Part 1): 
% The multiphase flow of CO2-CH4-H2-H2S gas mixtures. Journal of Natural 
% Gas Science and Engineering, 63, 85-94.
%
% 3. Nielsen, M. H., Arekhov, V., Whitson, C. H., Clemens, T., Zhainakov, 
% T., & Wegner, J. (2023, June). Fluid Modeling of Underground Hydrogen 
% Storage in a Depleted Natural Gas Field. In SPE Europec featured at EAGE 
% Conference and Exhibition (p. D021S001R002). SPE.
    
    % === NOTES:
    % 1. check the units
    % 2. volume shift method
    % 3. where did you get kij? H2-CH4: 0.362 H2-H2O: 0.1 CH4-H2O: 0.47893
    % 4. Check the Z equation solution
    % 5. Why P_app?
    % 6. what typo in Ganji2012? (A5) and (20)
    % 7. fugacity coefficient of H2O
    % 8. relation between molar volme and V
    
    R = 83.1447; % Gas constant [bar cm^3 mol^(-1) K^(-1)]
    
    % Defining 
    P_vec = [1:0.5:10, 11:50, 60:10:500, 550:50:1000]; % Pressure range [bar]
    y1_vec = 10.^(-5:0.01:0); % Mole fraction of H2 in vapor phase
    y2_vec = 1 - y1_vec;   % Mole fraction of CH4 in vapor phase
    
    P = zeros(length(P_vec), length(y1_vec));
    y1 = zeros(length(P_vec), length(y1_vec));
    
    for i = 1:length(y1_vec)
        y1(:, i) = y1_vec(i);
    end
    for i = 1:length(P_vec)
        P(i, :) = P_vec(i);
    end
    y2 = 1 - y1;
    T = ones(length(P_vec), length(y1_vec)) * T;
    
    % Temperature
    T_C = T - 273.15; % Temperature in [C]
    P_app = P - 1.01325 ; % Applied pressure at P_app = 0 bars, the absolute pressure is 1 atm or 1.01325 bars
    
    % brine
    % Salinity is in [ppm] (1000 [ppm] means 1 [g] solute in 1 [kg] water)
    m_Na = Salinity / 1e3 / 58.44; % [mol/kg]
    m_Cl = m_Na; % [mol/kg]
    
    % H2
    H2.Pc = 12.964; % Critical pressure [bar]
    H2.Tc = 33.145; % Critical temperature [K]
    H2.w = -0.219;  % Acentric factor
    H2.Mw = 2.0157; % Molecular weight [g/mol]
    H2.c = -2.959;  % Volume shift factor (from Shabani (2019))
    H2.X = 4.8070; % Constants fitting parameter for correct density
    H2.Vc = 64.4828; % the molar volume in the critical point. cm3/mol
    
    % CH4
    CH4.Pc = 45.992;  % Critical pressure [bar]
    CH4.Tc = 190.564; % Critical temperature [K]
    CH4.w = 0.01142;  % Acentric factor
    CH4.Mw = 16.0425; % Molecular weight [g/mol]
    CH4.c = -0.1557;  % Volume shift factor (from Shabani (2019))
    CH4.X = 6.404378; % Constants fitting parameter for correct density calculation for PR EOS.
    CH4.Vc = 99.3; % the molar volume in the critical point. cm3/mol
    
    % H2O
    H2O.Pc = 220.640; % Critical pressure [bar]
    H2O.Tc = 647.096; % Critical temperature [K]
    H2O.w = 0.3443;   % Acentric factor
    H2O.Mw = 18.0152; % Molecular weight [g/mol]
    H2O.X = 6.404378; % Constants fitting parameter for correct density calculation for PR EOS.
    H2O.Vc = 55.9478; % the molar volume in the critical point. cm3/mol
    
    % Peng-Robinson EOS
    alpha = @(T, Tc, w) (1 + (0.37464 + 1.54226*w - 0.26992*w^2) * (1 - sqrt(T / Tc))).^2;
    a = @(Tc, Pc, alpha) (0.45724 * R^2 * Tc^2 * alpha) / Pc;
    b = @(Tc, Pc) 0.07780 * R * Tc / Pc;
    
    H2.alpha = alpha(T, H2.Tc, H2.w);
    H2.a = a(H2.Tc, H2.Pc, H2.alpha);
    H2.b = b(H2.Tc, H2.Pc);
    H2.B = H2.b .* P ./ (R * T);
    
    CH4.alpha = alpha(T, CH4.Tc, CH4.w);
    CH4.a = a(CH4.Tc, CH4.Pc, CH4.alpha);
    CH4.b = b(CH4.Tc, CH4.Pc);
    CH4.B = CH4.b .* P ./ (R * T);
    
    
    H2O.alpha = alpha(T, H2O.Tc, H2O.w);
    H2O.a = a(H2O.Tc, H2O.Pc, H2O.alpha);
    H2O.b = b(H2O.Tc, H2O.Pc);
    H2O.B = H2O.b .* P ./ (R * T);
    
    % Mixture rule
    mixture.a11 = H2.a;
    mixture.a22 = CH4.a;
    mixture.a12 = sqrt(H2.a .* CH4.a) * (1 - 0.362);
    mixture.b1 = H2.b;
    mixture.b2 = CH4.b;
    
    mixture.a = y1.^2 .* mixture.a11 + 2 * y1 .* y2 .* mixture.a12 + y2.^2 .* mixture.a22;
    mixture.b = y1 * mixture.b1 + y2 * mixture.b2;
    
    mixture.A = mixture.a .* P ./ (R*T).^2;
    mixture.B = mixture.b .* P ./ (R*T);
    
    % Compressibility factor
    Z_guess = 1; % Initial guess
    Z = arrayfun(@(A, B)  fzero(@(Z) Z^3 - (1 - B) * Z^2 + (A - 2*B - 3*B^2) * Z - (A*B - B^2 - B^3), Z_guess),  mixture.A, mixture.B);
    if any(Z <= 0, 'all') || any(Z <= mixture.B, 'all')
        error('There is a problem in calculation of compressibility factor (Z).')
    end
    
    
    % Fugacity coefficient
    H2.phi = exp( ...
                  H2.B ./ mixture.B .* (Z - 1) ...
                - log(Z - mixture.B) ...
                + (mixture.A ./ (2.828 * mixture.B)) .* ...
                  (H2.B ./ mixture.B - (2 * (y1 .* mixture.a11 + y2 .* mixture.a12)) ./ mixture.a) .* ...
                  (log((Z + 2.414 * mixture.B) ./ (Z - 0.414 * mixture.B))) ...
                                                              );
    
    
    CH4.phi = exp( ...
                  CH4.B ./ mixture.B .* (Z - 1) ...
                - log(Z - mixture.B) ...
                + (mixture.A ./ (2.828 * mixture.B)) .* ...
                  (CH4.B ./ mixture.B - (2 * (y1 .* mixture.a12 + y2 .* mixture.a22)) ./ mixture.a) .* ...
                  (log((Z + 2.414 * mixture.B) ./ (Z - 0.414 * mixture.B))) ...
                                                              );
    
    % Properties of pure water
    V0 = (1 + 18.159725E-3 * T_C) ./ (0.9998396 + 18.224944E-3 * T_C - 7.922210E-6 * T_C.^2 - 55.44846E-9 * T_C.^3 + 149.7562E-12 * T_C.^4 - 393.2952E-15 * T_C.^5);
    B  = 19654.320 + 147.037   * T_C - 2.21554   * T_C.^2 + 1.0478E-2  * T_C.^3 - 2.2789E-5 * T_C.^4;
    A1 = 3.2891    - 2.3910E-3 * T_C + 2.8446E-4 * T_C.^2 - 2.8200E-6  * T_C.^3 + 8.477E-9  * T_C.^4;
    A2 = 6.2450E-5 - 3.9130E-6 * T_C + 3.4990E-8 * T_C.^2 - 7.9420E-10 * T_C.^3 + 3.299E-12 * T_C.^4;
    
    V = V0 - (V0 .* P_app) ./ (B + A1 .* P_app + A2 .* P_app.^2); % [cm3/g], typo in Ganji2012
    H2O.rho = 1 ./ V; % [g/cm3]
    
    % Fugucity of pure water (Appendix A)
    V_molar = 18.0152 * V;
    tau = 1 - T/H2O.Tc; 
    
    % vapor pressure of pure water, only applicable up to critical point
    P_s = H2O.Pc .* exp(H2O.Tc ./ T .*(-7.85951783 * tau + 1.84408259 * tau.^1.5 - 11.7866497 * tau.^3 + 22.6807411 * tau.^3.5 - 15.9618719 * tau.^4 + 1.80122502 * tau.^7.5)); 
    H2O.f0 = P_s .* exp((P - P_s) .* V_molar ./ (R*T)); % Poynting equation
    % H2O.phi = H2O.f ./ P;
    H2O.phi = exp( ...
                  H2O.B ./ mixture.B .* (Z - 1) ...
                - log(Z - mixture.B) ...
                + (mixture.A ./ (2.828 * mixture.B)) .* ...
                  (H2O.B ./ mixture.B - (2 * (y1 .* sqrt(H2.a .* H2O.a) .* (1 - 0.1) + y2 .* sqrt(CH4.a .* H2O.a) .* (1 - 0.47893))) ./ mixture.a) .* ...
                  (log((Z + 2.414 * mixture.B) ./ (Z - 0.414 * mixture.B))) ...
                                                              );
    
    
    % Henry's constant
    eta = 0.294994;
    tau = -4.23407;
    beta = 5.41551;
    H2.deltaB = tau + beta * sqrt(1E3 ./ T);
    H2.KH = exp((1 - eta) * log(H2O.f0) + eta * log(R * T / H2O.Mw .* H2O.rho) + 2 * H2O.rho .* H2.deltaB);
     
    eta = -0.092248;
    tau = -5.779280;
    beta = 7.262730;
    CH4.deltaB = tau + beta * sqrt(1E3 ./ T);
    CH4.KH = exp((1 - eta) * log(H2O.f0) + eta * log(R * T / H2O.Mw .* H2O.rho) + 2 * H2O.rho .* CH4.deltaB);
    
    % Activity coefficient
    Par = @(c, T, P) c(1) + c(2) * T + c(3) ./ T + c(4) * P + c(5) ./ P + c(6) * P ./ T + c(7) * T ./ P.^2 + c(8) * P ./ (630 - T) + c(9) * T .* log(P) + c(10) * P ./ T.^2;
    
    H2.c_Na =     [-2.1432831E+00 +3.1411257E-03 +3.9220546E+02 -2.8601200E-05 +0.0000000E+00 +2.3527160E-03 +0.0000000E+00 -2.4422000E-03 +2.9806000E-06 -3.8900000E-01];
    H2.c_Na_Cl =  [-4.0631000E-03 +0.0000000E+00 +0.0000000E+00 -3.6650000E-06 +1.7004000E-01 -4.1800000E-04 +5.8800000E-04 -3.9100000E-04 +0.0000000E+00 +1.9913000E-01];
    H2.gamma = exp(2 * m_Na * Par(H2.c_Na, T, P) + m_Na * m_Cl * Par(H2.c_Na_Cl, T, P));
    
    CH4.c_Na =    [-5.7066455E-01 +7.2997588E-04 +1.5176903E+02 +3.1927112E-05 +0.0000000E+00 -1.6426510E-05 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00];
    CH4.c_Na_Cl = [-2.9990084E-03 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00 +0.0000000E+00];
    CH4.gamma = exp(2 * m_Na * Par(CH4.c_Na, T, P) + m_Na * m_Cl * Par(CH4.c_Na_Cl, T, P));
    
    % K-value
    H2.K = H2.KH .* H2.gamma ./ (H2.phi .* P);
    
    CH4.K = CH4.KH .* CH4.gamma ./ (CH4.phi.* P);
    
    H2O.K0 = 10 .^ (-2.209 + 3.097E-2 * T_C - 1.098E-4 * T_C.^2 + 2.048E-7 * T_C.^3);
    p_ref = 1; % [bar];
    H2O.Vm_avg = 18.1; % [cm^3 mol^(-1)];
    H2O.K = H2O.K0 .* exp((P - p_ref) * H2O.Vm_avg ./ (R * T))./ (H2O.phi .* P); %(eq 11)
    
    H2.K = griddedInterpolant(P, y1, H2.K);
    CH4.K = griddedInterpolant(P, y1, CH4.K);
    H2O.K = griddedInterpolant(P, y1, H2O.K);

    Z = griddedInterpolant(P, y1, Z);

end