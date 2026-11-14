clear;
clc;
close all;

%% Stałe
C1 = 0.6;
C2 = 0.9;
alpha1 = 28;
alpha2 = 27;
FD0 = 10;    % cm^3/s, początkowa wartość zakłócenia
tau = 50;    % s, opóźnienie wejścia F1_in

%% Punkt pracy (warunek początkowy)
F1_0 = 60;
F2_0 = F1_0 + FD0;
h1_0 = (F2_0/alpha1)^2;
h2_0 = (F2_0/alpha2)^2;


%% Granice zbiorów rozmytych
F1_fuzzy_min = F1_0 * 0.3;
F2_fuzzy_min = F1_fuzzy_min + FD0;
h2_fuzzy_min = (F2_fuzzy_min/alpha2)^2;

F1_fuzzy_max = F1_0 * 1.7;
F2_fuzzy_max = F1_fuzzy_max + FD0;
h2_fuzzy_max = (F2_fuzzy_max/alpha2)^2;


%% Podział na modele lokalne
n_lok = 2; % ilość modeli lokalnych, od 2 do 5
h2_srodki = linspace(h2_fuzzy_min, h2_fuzzy_max, n_lok);

F_lok  = alpha2 * sqrt(h2_srodki);  
F1_lok = F_lok - FD0;
h1_lok = (F_lok / alpha1).^2;

%% Test modelu rozmytego
% 1. W rownowadze kazdego modelu lokalnego pochodne powinny byc ~0
fprintf('\nTest 1: rownowagi modeli lokalnych (n = %d)\n', n_lok);
for i = 1:n_lok
    x_rown = [C1*h1_lok(i)^3; C2*h2_srodki(i)^2];   % stan rownowagi modelu i
    dx_test = model_rozmyty(x_rown, F1_lok(i), FD0, FD0, ...
                            C1, C2, alpha1, alpha2, F1_lok, h2_srodki);
    fprintf('  model %d: h2 = %6.2f cm, F1 = %6.2f -> dx = [%9.2e; %9.2e]\n', ...
            i, h2_srodki(i), F1_lok(i), dx_test(1), dx_test(2));
end





% %% Rysowanie funkcji przynaleznosci
% h2_siatka = 0:0.01:20;     % [cm] punkty do wykresu
% 
% % wiersz = punkt siatki, kolumna = zbior rozmyty
% mu_siatka = zeros(length(h2_siatka), n_lok);
% for k = 1:length(h2_siatka)
%     mu_siatka(k, :) = przynaleznosc(h2_siatka(k), h2_srodki, "trojkat");
% end
% 
% fig = figure('Name', sprintf('Funkcje przynależności, n = %d', n_lok));
% plot(h2_siatka, mu_siatka, 'LineWidth', 1.5);
% hold on;
% xline(h2_srodki, 'k:');
% grid on;
% ylim([0 1.1]);
% xlabel('h_2 [cm]');
% ylabel('\mu');
% title(sprintf('Funkcje przynależności, n = %d', n_lok));
% legend(arrayfun(@(i) sprintf('zbiór %d', i), 1:n_lok, ...
%        'UniformOutput', false), 'Location', 'best');
% 
% %% Zapis do PNG
% folder = 'wykresy_png';
% if ~exist(folder, 'dir')
%     mkdir(folder);
% end
% exportgraphics(fig, fullfile(folder, sprintf('funkcje_przynaleznosci_n%d.png', n_lok)), ...
%                'Resolution', 300);