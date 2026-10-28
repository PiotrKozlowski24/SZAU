function [Ke, Mp] = dmc_macierze(D, N, Nu, lambda, s)
    s = s(:);
    s_ext = [s; s(D)*ones(N,1)];
    
    % macierz dynamiczna M (N x Nu)
    M = zeros(N, Nu);
    for j = 1:Nu
        M(j:N, j) = s_ext(1:N-j+1);
    end

    % macierz Mp (N x D-1)
    [I, J] = ndgrid(1:N, 1:D-1);
    Mp = s_ext(I+J) - s(J);

    % pierwszy wiersz macierzy wzmocnien
    K = (M'*M + lambda*eye(Nu)) \ M';
    Ke = K(1,:);
end