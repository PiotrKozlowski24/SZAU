function dx = model_procesu(x, F1, FD, C1, C2, a1, a2)
% Model dwoch zbiornikow, x = [V1; V2]
% Zwraca pochodne dx = [dV1/dt; dV2/dt]

    h1 = (max(x(1),0)/C1)^(1/3);
    h2 = (max(x(2),0)/C2)^(1/2);

    F2 = a1*sqrt(h1);
    F3 = a2*sqrt(h2);

    dx = [F1 + FD - F2;
          F2 - F3];
end