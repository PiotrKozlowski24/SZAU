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
n_lok = 5; % ilość modeli lokalnych, od 2 do 5
h2_srodki = linspace(h2_fuzzy_min, h2_fuzzy_max, n_lok);

F_lok  = alpha2 * sqrt(h2_srodki);  
F1_lok = F_lok - FD0;
h1_lok = (F_lok / alpha1).^2;   

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