clear; clc; close all;
%% Stale
C1 = 0.6;
C2 = 0.9;
alpha1 = 28;
alpha2 = 27;
FD  = 10;    % cm^3/s, staly dodatkowy doplyw
tau = 50;    % s, opoznienie wejscia F1_in

%% Plan skokow sterowania F1_in
t_skok  = [0   50  100  200  700  900  1100  1200];   % s
F1_wart = [60  65  75   100  50   30  20    60];    % cm^3/s

%% Punkt pracy
F1_0 = 60;
F2_0 = F1_0 + FD;
h1_0 = (F2_0/alpha1)^2;
h2_0 = (F2_0/alpha2)^2;
x0   = [C1*h1_0^3; C2*h2_0^2];

%% Symulacja
dt = 0.1;
t  = 0:dt:1400;
N  = numel(t);

F1 = F1_0*ones(1,N);
for i = 1:length(t_skok)
    F1(t >= t_skok(i) + tau) = F1_wart(i);
end

x  = zeros(2,N);  x(:,1)  = x0;   % nieliniowy
xl = zeros(2,N);  xl(:,1) = x0;   % liniowy
for n = 1:N-1
    k1 = model_procesu(x(:,n),          F1(n), FD, C1, C2, alpha1, alpha2);
    k2 = model_procesu(x(:,n)+dt*k1/2,  F1(n), FD, C1, C2, alpha1, alpha2);
    k3 = model_procesu(x(:,n)+dt*k2/2,  F1(n), FD, C1, C2, alpha1, alpha2);
    k4 = model_procesu(x(:,n)+dt*k3,    F1(n), FD, C1, C2, alpha1, alpha2);
    x(:,n+1) = max(x(:,n) + dt/6*(k1 + 2*k2 + 2*k3 + k4), 0);

    l1 = model_liniowy(xl(:,n),          F1(n), F1_0, FD, C1, C2, alpha1, alpha2);
    l2 = model_liniowy(xl(:,n)+dt*l1/2,  F1(n), F1_0, FD, C1, C2, alpha1, alpha2);
    l3 = model_liniowy(xl(:,n)+dt*l2/2,  F1(n), F1_0, FD, C1, C2, alpha1, alpha2);
    l4 = model_liniowy(xl(:,n)+dt*l3,    F1(n), F1_0, FD, C1, C2, alpha1, alpha2);
    xl(:,n+1) = xl(:,n) + dt/6*(l1 + 2*l2 + 2*l3 + l4);
end

%% Wyjscie h2
h2_nl  = (x(2,:)/C2).^(1/2);
h2_lin = h2_0 + (xl(2,:) - x0(2))/(2*C2*h2_0);   % dh2/dV2 = 1/(2*C2*h2_0)

err = h2_lin - h2_nl;

%% Wykresy
figure;
subplot(3,1,1);
plot(t, F1, 'LineWidth', 1.5); grid on;
ylabel('F_1 [cm^3/s]'); title('Sterowanie (z opoznieniem)');
subplot(3,1,2);
plot(t, h2_nl, 'LineWidth', 1.5); hold on;
plot(t, h2_lin, '--', 'LineWidth', 1.5); grid on;
ylabel('h_2'); legend('nieliniowy','liniowy','Location','best');
title('Zmienna regulowana h_2');
subplot(3,1,3);
yline(0, '--k', 'LineWidth', 1);
hold on;
plot(t, err, 'LineWidth', 1.5); hold on;
grid on;
xlabel('Czas [s]');
ylabel('h_{2,lin} - h_{2,nl}');
title('Blad linearyzacji');