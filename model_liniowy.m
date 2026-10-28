function dx = model_liniowy(x, F1, F1_0, FD, C1, C2, a1, a2)
% Model dwoch zbiornikow zlinearyzowany w punkcie pracy F1_0 i FD, x = [V1; V2]
% Zwraca pochodne dx = [dV1/dt; dV2/dt]
    % punkt pracy (stan ustalony)
    F2_0 = F1_0 + FD;
    h1_0 = (F2_0/a1)^2;
    h2_0 = (F2_0/a2)^2;
    V1_0 = C1*h1_0^3;
    V2_0 = C2*h2_0^2;
    % wspolczynniki Taylora: k2 = dF2/dV1, k3 = dF3/dV2
    k2 = a1/6 * C1^(-1/6) * V1_0^(-5/6);
    k3 = a2/4 * C2^(-1/4) * V2_0^(-3/4);
    % pochodne w zmiennych odchylkowych
    dV1 = x(1) - V1_0;
    dV2 = x(2) - V2_0;
    dx = [(F1 - F1_0) - k2*dV1;
          k2*dV1 - k3*dV2];
end