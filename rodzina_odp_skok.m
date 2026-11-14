clear; clc; close all;
%% Stałe
C1 = 0.6;
C2 = 0.9;
alpha1 = 28;
alpha2 = 27;
FD0 = 10;    % cm^3/s, wartość zakłócenia (stała)
tau = 50;    % s, opóźnienie wejścia F1_in
%% Wielkości skoków F1_in
t_skok = 50;                          % s, moment skoku
dF1    = [-50 -30 -15 -5 5 15 30 50]; % cm^3/s, wielkości skoków
kolory = lines(length(dF1));
%% Punkt pracy
F1_0 = 60;
F2_0 = F1_0 + FD0;
h1_0 = (F2_0/alpha1)^2;
h2_0 = (F2_0/alpha2)^2;
V1_0 = C1*h1_0^3;
V2_0 = C2*h2_0^2;
x0   = [V1_0; V2_0];
%% Czas
dt = 0.05;
t  = 0:dt:600;
N  = numel(t);
figure;
%% Pętla po wielkościach skoków
for j = 1:length(dF1)
    x_nl  = zeros(2,N);
    x_lin = zeros(2,N);
    x_nl(:,1)  = x0;
    x_lin(:,1) = x0;
    % wejście F1 (skok z opóźnieniem tau) i stałe zakłócenie FD
    F1 = F1_0*ones(1,N);
    F1(t >= t_skok + tau) = F1_0 + dF1(j);
    FD = FD0*ones(1,N);
    % RK4
    for n = 1:N-1
        k1 = model_procesu(x_nl(:,n),          F1(n), FD(n), C1, C2, alpha1, alpha2);
        k2 = model_procesu(x_nl(:,n)+dt*k1/2,  F1(n), FD(n), C1, C2, alpha1, alpha2);
        k3 = model_procesu(x_nl(:,n)+dt*k2/2,  F1(n), FD(n), C1, C2, alpha1, alpha2);
        k4 = model_procesu(x_nl(:,n)+dt*k3,    F1(n), FD(n), C1, C2, alpha1, alpha2);
        x_nl(:,n+1) = x_nl(:,n) + dt/6*(k1 + 2*k2 + 2*k3 + k4);

        k1 = model_liniowy(x_lin(:,n),          F1(n), FD(n), F1_0, FD0, C1, C2, alpha1, alpha2);
        k2 = model_liniowy(x_lin(:,n)+dt*k1/2,  F1(n), FD(n), F1_0, FD0, C1, C2, alpha1, alpha2);
        k3 = model_liniowy(x_lin(:,n)+dt*k2/2,  F1(n), FD(n), F1_0, FD0, C1, C2, alpha1, alpha2);
        k4 = model_liniowy(x_lin(:,n)+dt*k3,    F1(n), FD(n), F1_0, FD0, C1, C2, alpha1, alpha2);
        x_lin(:,n+1) = x_lin(:,n) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    end
    h2_nl  = (x_nl(2,:)/C2).^(1/2);
    h2_lin = h2_0 + (x_lin(2,:) - V2_0)/(2*C2*h2_0);
    % skoki ujemne na górnym wykresie, dodatnie na dolnym
    if dF1(j) < 0
        subplot(2,1,1);
    else
        subplot(2,1,2);
    end
    plot(t, h2_nl, 'Color', kolory(j,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('nieliniowy, \\DeltaF_1 = %+g', dF1(j))); hold on;
    plot(t, h2_lin, '--', 'Color', kolory(j,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('liniowy, \\DeltaF_1 = %+g', dF1(j)));
end
%% Wykresy
subplot(2,1,1);
grid on;
ylabel('h_2 [cm]');
title('Skoki ujemne: model nieliniowy (linia ciągła) i liniowy (linia przerywana)');
legend('Location', 'eastoutside');
subplot(2,1,2);
grid on;
xlabel('Czas [s]'); ylabel('h_2 [cm]');
title('Skoki dodatnie: model nieliniowy (linia ciągła) i liniowy (linia przerywana)');
legend('Location', 'eastoutside');
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