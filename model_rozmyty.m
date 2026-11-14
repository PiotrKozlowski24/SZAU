function dx = model_rozmyty(x, F1, FD, FD0, C1, C2, a1, a2, F1_lok, h2_srodki)
    h2 = (max(x(2),0)/C2)^(1/2);
    mu = przynaleznosc(h2, h2_srodki, 'trojkat');
    n = length(h2_srodki);
    dx = zeros(2,1);
    for i = 1:n
        dx_lok = model_liniowy(x,F1,FD,F1_lok(i),FD0,C1,C2,a1,a2);
        dx = dx_lok * mu(i) + dx;
    end

end
