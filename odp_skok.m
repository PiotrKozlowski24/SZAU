
function s = odp_skok(C1, C2, alpha1, alpha2, FD, F1_0, tau, dt, D)
    %% Punkt pracy
    F2_0 = F1_0 + FD;

    h1_0 = (F2_0 / alpha1)^2;
    h2_0 = (F2_0 / alpha2)^2;

    V1_0 = C1 * h1_0^3;
    V2_0 = C2 * h2_0^2;

    %% Wspolczynniki modelu liniowego
    k2 = alpha1/6 * C1^(-1/6) * V1_0^(-5/6);
    k3 = alpha2/4 * C2^(-1/4) * V2_0^(-3/4);

    %% Czas
    t = (0:D-1) * dt;

    %% Jednostkowy skok wejscia z opoznieniem
    dF1_tau = zeros(1, D);
    dF1_tau(t >= tau) = 1;

    %% Symulacja modelu liniowego metoda RK4
    dx = zeros(2, D);

    for n = 1:D-1
        u = dF1_tau(n);

        k1v = model_lin(dx(:,n), u, k2, k3);
        k2v = model_lin(dx(:,n) + dt*k1v/2, u, k2, k3);
        k3v = model_lin(dx(:,n) + dt*k2v/2, u, k2, k3);
        k4v = model_lin(dx(:,n) + dt*k3v, u, k2, k3);

        dx(:,n+1) = dx(:,n) + ...
            dt/6 * (k1v + 2*k2v + 2*k3v + k4v);
    end

    %% Odpowiedz skokowa poziomu h2
    s = (dx(2,:) / (2*C2*h2_0)).';

end

%% Model liniowy w zmiennych przyrostowych
function dx = model_lin(x, du, k2, k3)

    dx = [du - k2*x(1);
          k2*x(1) - k3*x(2)];

end
