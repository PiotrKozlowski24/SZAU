function [mu] = przynaleznosc(h2, h2_srodki, rodzaj_funkcji)
% Zwraca wektor przynaleznosci mu (1 x n) wartosci h2 do zbiorow rozmytych
% o srodkach h2_srodki
    if strcmp(rodzaj_funkcji, "trojkat")
        % przyciecie h2 do zakresu srodkow 
        h2 = min(max(h2, h2_srodki(1)), h2_srodki(end));

        koszyk = discretize(h2, h2_srodki);
        mu = zeros(1, length(h2_srodki));

        % funkcja trojkatna
        dystans = h2_srodki(koszyk + 1) - h2_srodki(koszyk);
        mu(koszyk + 1) = (h2 - h2_srodki(koszyk)) / dystans;
        mu(koszyk)     = (h2_srodki(koszyk + 1) - h2) / dystans;
    else
        error('Nieobslugiwany rodzaj funkcji przynaleznosci: "%s"', rodzaj_funkcji);
    end
end