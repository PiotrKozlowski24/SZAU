function u = regulator_dmc(Ke, Mp, u_hist, y_hist, y_zad)

    N = size(Mp, 1);
    D = size(Mp, 2) + 1;

    % przeszle przyrosty: dup(j) = u(k-j) - u(k-j-1), j = 1..D-1
    uh  = u_hist(:);
    dup = flipud(diff(uh(end-D+1:end)));

    % odpowiedz swobodna i uchyb
    y0 = y_hist(end)*ones(N,1) + Mp*dup;
    e  = y_zad(:).*ones(N,1) - y0;

    % przyrost sterowania i sterowanie
    du = Ke*e;
    u  = uh(end) + du;
end