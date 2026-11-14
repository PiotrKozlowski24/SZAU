function dx = model_liniowy(x, F1, FD, F1_0, FD0, C1, C2, a1, a2)
% Model zlinearyzowany w punkcie pracy F1_0 i FD0, x = [V1; V2]
    F2_0 = F1_0 + FD0;
    h1_0 = (F2_0/a1)^2;
    h2_0 = (F2_0/a2)^2;
    V1_0 = C1*h1_0^3;
    V2_0 = C2*h2_0^2;
    k2 = a1/6 * C1^(-1/6) * V1_0^(-5/6);
    k3 = a2/4 * C2^(-1/4) * V2_0^(-3/4);
    dV1 = x(1) - V1_0;
    dV2 = x(2) - V2_0;
    dx = [(F1 - F1_0) + (FD - FD0) - k2*dV1;
          k2*dV1 - k3*dV2];
end