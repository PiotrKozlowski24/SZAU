clear; clc; close all;
%% Stałe
C1 = 0.6;
C2 = 0.9;
alpha1 = 28;
alpha2 = 27;
FD0 = 10;    % cm^3/s, początkowa wartość zakłócenia
tau = 50;    % s, opóźnienie wejścia F1_in

%% Plan skoków sterowania F1_in
t_skok  = [0   50  150  250  350 550 700];
F1_wart = [60  70   10   60  90   60   30];

%% Plan skoków zakłócenia FD
tFD_skok = [850 900 950];    % s, moment skoku zakłócenia
FD_wart  = [20 30 5];     % cm^3/s, wartość po skoku

%% Punkt pracy (warunek początkowy)
F1_0 = 60;
F2_0 = F1_0 + FD0;
h1_0 = (F2_0/alpha1)^2;
h2_0 = (F2_0/alpha2)^2;
x0   = [C1*h1_0^3; C2*h2_0^2];

%% Symulacja
dt = 0.01;
t  = 0:dt:1000;
N  = numel(t);
x  = zeros(2,N);                 % x = [V1; V2]
x(:,1) = x0;

%% Sterowanie F1 jako wektor (z opóźnieniem tau)
F1_in = F1_0*ones(1,N);
F1    = F1_0*ones(1,N);
for i = 1:length(t_skok)
    F1_in(t >= t_skok(i)) = F1_wart(i);
    F1(t >= t_skok(i) + tau) = F1_wart(i);
end

%% Zakłócenie FD jako wektor (bez opóźnienia)
FD = FD0*ones(1,N);
for i = 1:length(tFD_skok)
    FD(t >= tFD_skok(i)) = FD_wart(i);
end

% RK4
for n = 1:N-1
    k1 = model_procesu(x(:,n),          F1(n), FD(n), C1, C2, alpha1, alpha2);
    k2 = model_procesu(x(:,n)+dt*k1/2,  F1(n), FD(n), C1, C2, alpha1, alpha2);
    k3 = model_procesu(x(:,n)+dt*k2/2,  F1(n), FD(n), C1, C2, alpha1, alpha2);
    k4 = model_procesu(x(:,n)+dt*k3,    F1(n), FD(n), C1, C2, alpha1, alpha2);
    x(:,n+1) = x(:,n) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

%% Wyniki
h1 = (x(1,:)/C1).^(1/3);
h2 = (x(2,:)/C2).^(1/2); % zmienna regulowana

figure;
subplot(3,1,1);
plot(t, F1_in, 'Color', [0 0.447 0.741], 'LineWidth', 2); hold on;
plot(t, F1, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 0.8);
grid on;
ylabel('F [cm^3/s]'); title('Sterowanie F_{1,in} i przepływ F_1 (z opóźnieniem)');
legend('F_{1,in}', 'F_1', 'Location', 'best');

subplot(3,1,2);
plot(t, FD, 'Color', [0.929 0.694 0.125], 'LineWidth', 2); grid on;
ylabel('F_D [cm^3/s]'); title('Zakłócenie F_D');

subplot(3,1,3);
plot(t, h2, 'Color', [0.466 0.674 0.188], 'LineWidth', 2); grid on;
xlabel('Czas [s]'); ylabel('h_2 [cm]'); title('Zmienna regulowana h_2');

%% Model
function dx = model_procesu(x, F1, FD, C1, C2, a1, a2)
% Model nieliniowy, x = [V1; V2]
    h1 = (x(1)/C1)^(1/3);
    h2 = (x(2)/C2)^(1/2);
    F2 = a1*sqrt(h1);
    F3 = a2*sqrt(h2);
    dx = [F1 + FD - F2;
          F2 - F3];
end