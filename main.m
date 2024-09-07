%% SCRIPT FOR TESTING LIGHT PROPAGATION THROUGH BRAGG WITH GLASS
% -------------------------------------------------------------------------
%% GENERAL INPUT DATA
% -------------------------------------------------------------------------
lam = linspace(380, 1080, 851);

Bragg1 = ["MgF2"];  % First Bragg stack
Bragg2 = ["SiO2"];  % Second Bragg stack
matPSC = [Bragg1, "gls", Bragg2];  % Combining Bragg1 and Bragg2 with glass in between
dgls=1000
dBragg1 = [110.18];  % Thicknesses for Bragg1 (SiO2/TiO2 x2)
dBragg2 = [100];  % Thicknesses for Bragg2 (TiO2/SiO2 x2)
dPSC = [dBragg1, dgls, dBragg2];
incoh = 1e3; % incoherent layer is 1000nm thickness
theta = 0;

% -------------------------------------------------------------------------
% PSC with AZO and EVA on PV
% -------------------------------------------------------------------------
materials = ["air",matPSC,"air","air","air","air"];
d = [dPSC,dAZO,dEVA,dARC];
stack = set_stack( materials, d, lam, theta, incoh );
[R_PSC_AZO_EVA_PV, T_PSC_AZO_EVA_PV, ~] = ATR1D( stack );
A_PSC_AZO_EVA_PV = 1 - T_PSC_AZO_EVA_PV.sp - R_PSC_AZO_EVA_PV.sp;  % Absorptance

Jsc_PSC_AZO_EVA_PV = get_electric(lam, T_PSC_AZO_EVA_PV.sp).jsc;

% -------------------------------------------------------------------------
% PLOTTING RESULTS
% -------------------------------------------------------------------------
figure();
plot(lam, T_PSC_AZO_EVA_PV.sp, 'b', 'LineWidth', 1.5); % Transmittance
hold on;
plot(lam, R_PSC_AZO_EVA_PV.sp, 'r', 'LineWidth', 1.5); % Reflectance
plot(lam, A_PSC_AZO_EVA_PV, 'g', 'LineWidth', 1.5);    % Absorptance
hold off;

xlabel('Wavelength [nm]');
ylabel('Fraction');
xlim([lam(1) lam(end)]);
ylim([0 1]);
legend('Transmittance (T)', 'Reflectance (R)', 'Absorptance (A)');
title('TRA Plot for PSC Structure with Separate Bragg Stacks');

disp(['Jsc for classical PSC: ' num2str(Jsc_PSC_AZO_EVA_PV)]);

