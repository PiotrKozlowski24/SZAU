clear; clc; close all;
%% Stałe
alpha2 = 27;
FD_wart = [0 5 10 20 30];     % cm^3/s, wartości zakłócenia FD
F1 = 0:0.5:100;               % cm^3/s, zakres F1_in
%% Charakterystyki statyczne
% W stanie ustalonym: F3 = F2 = F1 + FD, a F3 = alpha2*sqrt(h2)
% => h2 = ((F1 + FD)/alpha2)^2
figure; hold on;
for i = 1:length(FD_wart)
    h2 = ((F1 + FD_wart(i))/alpha2).^2;
    plot(F1, h2, 'LineWidth', 2, ...
        'DisplayName', sprintf('F_D = %g cm^3/s', FD_wart(i)));
end
%% Punkt pracy (F1 = 60, FD = 10)
h2_0 = ((60 + 10)/alpha2)^2;
plot(60, h2_0, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'punkt pracy');
grid on;
xlabel('F_{1,in} [cm^3/s]'); ylabel('h_2 [cm]');
title('Charakterystyki statyczne h_2(F_{1,in}) dla różnych F_D');
legend('Location', 'northwest');