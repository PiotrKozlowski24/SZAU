clear; clc; close all;

%% Stałe
C1 = 0.6;
C2 = 0.9;
alpha1 = 28;
alpha2 = 27;
FD  = 10;    % cm^3/s, stały dodatkowy dopływ
tau = 50;    % s, opóźnienie wejścia F1_in

%% Plan skoków sterowania F1_in
t_skok  = [0   150  300  450];   % s, momenty zadania skoków
F1_wart = [60  80   40   60];    % cm^3/s, wartości po skoku

%% Symulacja
dt = 0.1;
t  = 0:dt:600;
N  = numel(t);
x  = zeros(2,N);                 % x = [V1; V2]

%% Sterowanie F1 jako wektor (z opóźnieniem tau)
F1 = zeros(1,N);
for i = 1:length(t_skok)
    F1(t >= t_skok(i) + tau) = F1_wart(i);
end

% RK4
for n = 1:N-1
    k1 = model_procesu(x(:,n),          F1(n), FD, C1, C2, alpha1, alpha2);
    k2 = model_procesu(x(:,n)+dt*k1/2,  F1(n), FD, C1, C2, alpha1, alpha2);
    k3 = model_procesu(x(:,n)+dt*k2/2,  F1(n), FD, C1, C2, alpha1, alpha2);
    k4 = model_procesu(x(:,n)+dt*k3,    F1(n), FD, C1, C2, alpha1, alpha2);
    x(:,n+1) = x(:,n) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

%% Wyniki
h1 = (x(1,:)/C1).^(1/3);
h2 = (x(2,:)/C2).^(1/2); % zmienna regulowana

figure;
subplot(2,1,1);
plot(t, F1, 'LineWidth', 1.5); grid on;
ylabel('F_1 [cm^3/s]'); title('Sterowanie (z opóźnieniem)');
subplot(2,1,2);
plot(t, h2, 'LineWidth', 1.5); grid on;
xlabel('Czas [s]'); ylabel('h_2'); title('Zmienna regulowana');