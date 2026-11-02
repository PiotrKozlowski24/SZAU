clear; clc; close all;
%% Stałe
C1 = 0.6;
C2 = 0.9;
alpha1 = 28;
alpha2 = 27;
FD0 = 10;    % cm^3/s, wartość zakłócenia w punkcie pracy
tau = 50;    % s, opóźnienie wejścia F1_in

%% Plan skoków sterowania F1_in
t_skok  = [0   50  150  250  350 550 700];
F1_wart = [60  70   10   60  90   60   30];

%% Plan skoków zakłócenia FD
tFD_skok = [850 900 950];    % s, moment skoku zakłócenia
FD_wart  = [20 30 5];     % cm^3/s, wartość po skoku

%% Punkt pracy
F1_0 = 60;
F2_0 = F1_0 + FD0;
h1_0 = (F2_0/alpha1)^2;
h2_0 = (F2_0/alpha2)^2;
V1_0 = C1*h1_0^3;
V2_0 = C2*h2_0^2;
x0   = [V1_0; V2_0];

%% Symulacja
dt = 0.01;
t  = 0:dt:1000;
N  = numel(t);
x_nl  = zeros(2,N);              % model nieliniowy, x = [V1; V2]
x_lin = zeros(2,N);              % model liniowy
x_nl(:,1)  = x0;
x_lin(:,1) = x0;

%% Sygnał sterujący F1_in i opóźnione wejście F1
F1_in = F1_0*ones(1,N);
for i = 1:length(t_skok)
    F1_in(t >= t_skok(i)) = F1_wart(i);
end

% Opóźnienie wejścia o tau
F1 = F1_0*ones(1,N);
for i = 1:length(t_skok)
    F1(t >= t_skok(i) + tau) = F1_wart(i);
end

%% Zakłócenie FD (bez opóźnienia)
FD = FD0*ones(1,N);
for i = 1:length(tFD_skok)
    FD(t >= tFD_skok(i)) = FD_wart(i);
end

%% RK4
for n = 1:N-1
    % model nieliniowy
    k1 = model_procesu(x_nl(:,n),          F1(n), FD(n), C1, C2, alpha1, alpha2);
    k2 = model_procesu(x_nl(:,n)+dt*k1/2,  F1(n), FD(n), C1, C2, alpha1, alpha2);
    k3 = model_procesu(x_nl(:,n)+dt*k2/2,  F1(n), FD(n), C1, C2, alpha1, alpha2);
    k4 = model_procesu(x_nl(:,n)+dt*k3,    F1(n), FD(n), C1, C2, alpha1, alpha2);
    x_nl(:,n+1) = x_nl(:,n) + dt/6*(k1 + 2*k2 + 2*k3 + k4);

    % model liniowy
    k1 = model_liniowy(x_lin(:,n),          F1(n), FD(n), F1_0, FD0, C1, C2, alpha1, alpha2);
    k2 = model_liniowy(x_lin(:,n)+dt*k1/2,  F1(n), FD(n), F1_0, FD0, C1, C2, alpha1, alpha2);
    k3 = model_liniowy(x_lin(:,n)+dt*k2/2,  F1(n), FD(n), F1_0, FD0, C1, C2, alpha1, alpha2);
    k4 = model_liniowy(x_lin(:,n)+dt*k3,    F1(n), FD(n), F1_0, FD0, C1, C2, alpha1, alpha2);
    x_lin(:,n+1) = x_lin(:,n) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

%% Wyniki
h2_nl  = (x_nl(2,:)/C2).^(1/2);                          % zmienna regulowana
h2_lin = h2_0 + (x_lin(2,:) - V2_0)/(2*C2*h2_0);         % wzór zlinearyzowany

figure;
subplot(3,1,1);
plot(t, F1_in, 'Color', [0 0.447 0.741], 'LineWidth', 2); hold on;
plot(t, F1, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 0.8);
grid on;
ylabel('F [cm^3/s]');
title('Sterowanie F_{1,in} i przepływ F_1 (z opóźnieniem)');
legend('F_{1,in}', 'F_1', 'Location', 'best');

subplot(3,1,2);
plot(t, FD, 'Color', [0.929 0.694 0.125], 'LineWidth', 2); grid on;
ylabel('F_D [cm^3/s]'); title('Zakłócenie F_D');

subplot(3,1,3);
plot(t, h2_nl, 'Color', [0.466 0.674 0.188], 'LineWidth', 2); hold on;
plot(t, h2_lin, '--', 'Color', [0.494 0.184 0.556], 'LineWidth', 1.5);
grid on;
xlabel('Czas [s]'); ylabel('h_2 [cm]'); title('Zmienna regulowana h_2');
legend('nieliniowy', 'liniowy', 'Location', 'best');

%% Modele
function dx = model_liniowy(x, F1, FD, F1_0, FD0, C1, C2, a1, a2)
% Model zlinearyzowany w punkcie pracy F1_0 i FD0, x = [V1; V2]
    F2_0 = F1_0 + FD0;
    h1_0 = (F2_0/a1)^2;
    h2_0 = (F2_0/a2)^2;
    V1_0 = C1*h1_0^3;
    V2_0 = C2*h2_0^2;
    k2 = a1/6 * C1^(-1/6) * V1_0^(-5/6);
    k3 = a2/4 * C2^(-1/4) * V2_0^(-3/4);
    dV1 = x(1) - V1_0;
    dV2 = x(2) - V2_0;
    dx = [(F1 - F1_0) + (FD - FD0) - k2*dV1;
          k2*dV1 - k3*dV2];
end

function dx = model_procesu(x, F1, FD, C1, C2, a1, a2)
% Model nieliniowy, x = [V1; V2]
    h1 = (x(1)/C1)^(1/3);
    h2 = (x(2)/C2)^(1/2);
    F2 = a1*sqrt(h1);
    F3 = a2*sqrt(h2);
    dx = [F1 + FD - F2;
          F2 - F3];
end