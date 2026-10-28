clear; clc; close all;
%% Stałe
C1 = 0.6;
C2 = 0.9;
alpha1 = 28;
alpha2 = 27;
FD  = 10;    % cm^3/s, stały dodatkowy dopływ
tau = 50;    % s, opóźnienie wejścia F1_in
%% Plan skoków sterowania F1_in
t_skok  = [0   50  150  250  350 450 700];   % s, momenty zadania skoków
F1_wart = [60  65   55   70  50  90   30];    % cm^3/s, wartości po skoku
%% Punkt pracy
F1_0 = 60;
F2_0 = F1_0 + FD;
h1_0 = (F2_0/alpha1)^2;
h2_0 = (F2_0/alpha2)^2;
x0   = [C1*h1_0^3; C2*h2_0^2];
%% Symulacja
dt = 0.1;
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

%% RK4
for n = 1:N-1
    % model nieliniowy
    k1 = model_procesu(x_nl(:,n),          F1(n), FD, C1, C2, alpha1, alpha2);
    k2 = model_procesu(x_nl(:,n)+dt*k1/2,  F1(n), FD, C1, C2, alpha1, alpha2);
    k3 = model_procesu(x_nl(:,n)+dt*k2/2,  F1(n), FD, C1, C2, alpha1, alpha2);
    k4 = model_procesu(x_nl(:,n)+dt*k3,    F1(n), FD, C1, C2, alpha1, alpha2);
    x_nl(:,n+1) = x_nl(:,n) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    % model liniowy
    k1 = model_liniowy(x_lin(:,n),          F1(n), F1_0, FD, C1, C2, alpha1, alpha2);
    k2 = model_liniowy(x_lin(:,n)+dt*k1/2,  F1(n), F1_0, FD, C1, C2, alpha1, alpha2);
    k3 = model_liniowy(x_lin(:,n)+dt*k2/2,  F1(n), F1_0, FD, C1, C2, alpha1, alpha2);
    k4 = model_liniowy(x_lin(:,n)+dt*k3,    F1(n), F1_0, FD, C1, C2, alpha1, alpha2);
    x_lin(:,n+1) = x_lin(:,n) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end
%% Wyniki
h2_nl  = (x_nl(2,:)/C2).^(1/2);                  % zmienna regulowana
h2_lin = (max(x_lin(2,:),0)/C2).^(1/2);
figure;
subplot(2,1,1);
plot(t, F1_in, 'LineWidth', 1.5);
grid on;
ylabel('F_{1,in} [cm^3/s]');
title('Sygnał sterujący F_{1,in}');
subplot(2,1,2);
plot(t, h2_nl, t, h2_lin, '--', 'LineWidth', 1.5); grid on;
xlabel('Czas [s]'); ylabel('h_2'); title('Zmienna regulowana h_2');
legend('nieliniowy', 'liniowy');

%% Modele
function dx = model_liniowy(x, F1, F1_0, FD, C1, C2, a1, a2)
    F2_0 = F1_0 + FD;
    h1_0 = (F2_0/a1)^2;
    h2_0 = (F2_0/a2)^2;
    V1_0 = C1*h1_0^3;
    V2_0 = C2*h2_0^2;
    k2 = a1/6 * C1^(-1/6) * V1_0^(-5/6);
    k3 = a2/4 * C2^(-1/4) * V2_0^(-3/4);
    dV1 = x(1) - V1_0;
    dV2 = x(2) - V2_0;
    dx = [(F1 - F1_0) - k2*dV1;
          k2*dV1 - k3*dV2];
end

function dx = model_procesu(x, F1, FD, C1, C2, a1, a2)
    h1 = (max(x(1),0)/C1)^(1/3);
    h2 = (max(x(2),0)/C2)^(1/2);
    F2 = a1*sqrt(h1);
    F3 = a2*sqrt(h2);
    dx = [F1 + FD - F2;
          F2 - F3];
end