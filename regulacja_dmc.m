clear;
clc;
close all;

%% Stale i parametry obiektu
C1 = 0.6;
C2 = 0.9;
alpha1 = 28;
alpha2 = 27;

FD0 = 10;          % [cm^3/s] nominalny dodatkowy doplyw (punkt pracy)
tau = 50;          % [s] opoznienie wejscia F1
dt = 1;            % [s] krok probkowania

%% Punkt pracy
F1_0 = 60;         % [cm^3/s] nominalne sterowanie
F2_0 = F1_0 + FD0;

h1_0 = (F2_0 / alpha1)^2;
h2_0 = (F2_0 / alpha2)^2;

% Zmienne stanu: objetosci cieczy w zbiornikach
x0 = [C1*h1_0^3;
      C2*h2_0^2];

%% Czas symulacji
t_max = 2000;      % [s]
K = round(t_max / dt);
t = (0:K)' * dt;

%% Odpowiedz skokowa obiektu liniowego (dla FD0)
D = 200;           % horyzont dynamiki
s = odp_skok(C1, C2, alpha1, alpha2, ...
             FD0, F1_0, tau, dt, D);

%% Wykres odpowiedzi skokowej s
ts = (1:D)' * dt;

figure('Name', 'Odpowiedź skokowa obiektu');
stairs(ts, s, '--', 'Color', [0 0.447 0.741], 'LineWidth', 0.8); hold on;
scatter(ts, s, 12, [0.85 0.325 0.098], 'filled');
grid on;
xlabel('Czas t [s]');
ylabel('s');
title('Odpowiedź skokowa obiektu');
legend('odpowiedź skokowa', 'punkty próbkowania', 'Location', 'southeast');
    
%% Parametry regulatora DMC
N = 100;           % horyzont predykcji
Nu = 10;           % horyzont sterowania
lambda = 10;      % wspolczynnik kary za zmiany sterowania

% Macierze regulatora liczone raz
[Ke, Mp] = dmc_macierze(D, N, Nu, lambda, s);

%% Ograniczenia sterowania
F1_min = 0;        % [cm^3/s]
F1_max = Inf;

%% Liczba probek odpowiadajaca opoznieniu
delay_steps = round(tau / dt);

if abs(delay_steps*dt - tau) > 1e-9
    error('Opoznienie tau musi byc wielokrotnoscia kroku dt.');
end

%% Scenariusze
% Scenariusz 1: skoki wartosci zadanej, FD stale
t_zad = [0 100 300 500 1000 1600];        % [s] chwile skokow
w_zad = [h2_0 8 5 10 3 h2_0];             % [cm] kolejne wartosci zadane h2

% Scenariusz 2: dwie wartosci zadane (punkt pracy i inna), po 2 skoki FD dla kazdej
t_zad2 = [0 1100];                        % [s] chwile zmiany wartosci zadanej
w_zad2 = [h2_0 13];                        % [cm] punkt pracy, potem inna wartosc

t_FD = [0 200 500 800 2300 3100 ];             % [s] chwile skokow FD
w_FD = [10 15 5 10 15 5];                   % [cm^3/s] kolejne wartosci FD
% h2_zad = h2_0: skoki FD w 250 s i 600 s
% h2_zad = 9:    skoki FD w 1300 s i 1650 s

% Przebiegi w czasie dla obu scenariuszy (kolumna = scenariusz)
h2_zad_all = h2_0 * ones(K+1, 2);
FD_all = FD0 * ones(K+1, 2);

for i = 1:length(t_zad)
    h2_zad_all(t >= t_zad(i), 1) = w_zad(i);
end

for i = 1:length(t_zad2)
    h2_zad_all(t >= t_zad2(i), 2) = w_zad2(i);
end

for i = 1:length(t_FD)
    FD_all(t >= t_FD(i), 2) = w_FD(i);
end

% Wyniki
h2_nl_all = zeros(K+1, 2);
F1_cmd_all = zeros(K+1, 2);
F1_applied_all = zeros(K+1, 2);

%% Petla po scenariuszach
for sc = 1:2

    h2_zad_vec = h2_zad_all(:, sc);
    FD_vec = FD_all(:, sc);

    % Inicjalizacja symulacji
    x_nl = zeros(2, K+1);
    x_nl(:,1) = x0;

    % Historie w wartosciach bezwzglednych (start z punktu pracy)
    u_hist = F1_0 * ones(D,1);
    y_hist = h2_0 * ones(D,1);

    % Sterowanie zadane przez regulator i doprowadzone do obiektu
    F1_cmd = F1_0 * ones(K+1,1);
    F1_applied = F1_0 * ones(K+1,1);

    %% Petla symulacyjna
    for n = 1:K

        %% 1. Pomiar wyjscia obiektu nieliniowego
        h2_meas = sqrt(max(x_nl(2,n), 0) / C2);
        y_hist = [y_hist(2:end); h2_meas];

        %% 2. Sterowanie z regulatora DMC
        y_zad = h2_zad_vec(n);

        u = regulator_dmc(Ke, Mp, u_hist, y_hist, y_zad);

        % Ograniczenie rzeczywistego przeplywu F1
        u = min(max(u, F1_min), F1_max);

        F1_cmd(n) = u;
        u_hist = [u_hist(2:end); u];
        
        %% 3. Opoznienie wejscia
        idx_opoz = n - delay_steps;

        if idx_opoz >= 1
            F1_in = F1_cmd(idx_opoz);
        else
            F1_in = F1_0;
        end

        F1_applied(n) = F1_in;

        %% 4. Obiekt nieliniowy, RK4 (FD zmienne w czasie)
        FD = FD_vec(n);

        k1 = model_procesu(x_nl(:,n), ...
                           F1_in, FD, C1, C2, alpha1, alpha2);
        k2 = model_procesu(x_nl(:,n) + dt*k1/2, ...
                           F1_in, FD, C1, C2, alpha1, alpha2);
        k3 = model_procesu(x_nl(:,n) + dt*k2/2, ...
                           F1_in, FD, C1, C2, alpha1, alpha2);
        k4 = model_procesu(x_nl(:,n) + dt*k3, ...
                           F1_in, FD, C1, C2, alpha1, alpha2);

        x_nl(:,n+1) = x_nl(:,n) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);

        % Zabezpieczenie przed ujemnymi objetosciami
        x_nl(:,n+1) = max(x_nl(:,n+1), 0);

    end

    % Uzupelnienie ostatnich probek do wykresow
    F1_cmd(end) = F1_cmd(end-1);
    F1_applied(end) = F1_applied(end-1);

    % Zapamietanie wynikow scenariusza
    h2_nl_all(:,sc) = sqrt(max(x_nl(2,:), 0) ./ C2).';
    F1_cmd_all(:,sc) = F1_cmd;
    F1_applied_all(:,sc) = F1_applied;

end

%% Wykres 1: skoki wartosci zadanej
figure('Name', 'DMC - skoki wartości zadanej');
 
subplot(2,1,1);
plot(t, h2_nl_all(:,1), 'b', 'LineWidth', 1.5);
hold on;
stairs(t, h2_zad_all(:,1), 'r--', 'LineWidth', 1.2);
yline(h2_0, 'k:', 'LineWidth', 1.2);
grid on;
xlabel('Czas t [s]');
ylabel('Poziom h_2 [cm]');
title('Odpowiedź układu regulacji DMC - skoki wartości zadanej');
legend('h_2 - wartość regulowana', 'h_{2,zad} - wartość zadana', ...
       'h_{2,0} - punkt pracy');
 
subplot(2,1,2);
stairs(t, F1_cmd_all(:,1), 'b', 'LineWidth', 1.3);
hold on;
yline(F1_0, 'k:', 'LineWidth', 1.2);
grid on;
xlabel('Czas t [s]');
ylabel('Przepływ F_{1,in} [cm^3/s]');
title('Sterowanie regulatora');
legend('F_{1,in}', 'F_{1,0} - punkt pracy');
 
%% Wykres 2: skoki zaklocenia FD
figure('Name', 'DMC - skoki zaklocenia FD');
 
subplot(3,1,1);
plot(t, h2_nl_all(:,2), 'b', 'LineWidth', 1.5);
hold on;
stairs(t, h2_zad_all(:,2), 'r--', 'LineWidth', 1.2);
yline(h2_0, 'k:', 'LineWidth', 1.2);
grid on;
xlabel('Czas t [s]');
ylabel('Poziom h_2 [cm]');
title('Odpowiedź układu regulacji DMC - skoki zakłócenia F_D');
legend('h_2 - wartość regulowana', 'h_{2,zad} - wartość zadana', ...
       'h_{2,0} - punkt pracy');
 
subplot(3,1,2);
stairs(t, F1_cmd_all(:,2), 'b', 'LineWidth', 1.3);
hold on;
yline(F1_0, 'k:', 'LineWidth', 1.2);
grid on;
xlabel('Czas t [s]');
ylabel('Przepływ F_{1,in} [cm^3/s]');
title('Sterowanie regulatora');
legend('F_{1,in}', 'F_{1,0} - punkt pracy');
 
subplot(3,1,3);
stairs(t, FD_all(:,2), 'k', 'LineWidth', 1.3);
hold on;
yline(FD0, 'k:', 'LineWidth', 1.2);
grid on;
xlabel('Czas t [s]');
ylabel('F_D [cm^3/s]');
title('Zakłócenie F_D');
legend('F_D', 'F_{D,0} - punkt pracy');

%% Informacje koncowe
fprintf('Poczatkowy poziom h2: %.4f cm\n', h2_0);
fprintf('Scenariusz 1 - koncowe h2: %.4f cm (zadane %.4f cm)\n', ...
        h2_nl_all(end,1), h2_zad_all(end,1));
fprintf('Scenariusz 2 - koncowe h2: %.4f cm (zadane %.4f cm)\n', ...
        h2_nl_all(end,2), h2_zad_all(end,2));