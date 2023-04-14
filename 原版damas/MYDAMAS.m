function DAMAS_result = MYDAMAS(DAS_result, a, maxIter)

    temp = abs(DAS_result);

    [N_X, N_Y, N_mic] = size(a);
    en = reshape(a, N_X * N_Y, N_mic);
    A = (abs(en * en') .^ 2) ./ N_mic ^ 2;

    Q = zeros(size(temp)); Q0 = temp;

    for i = 1:maxIter

        for n = 1:N_X * N_Y
            Q(n) = max(0, temp(n) - A(n, 1:n - 1) * Q(1:n - 1)' ...
                - A(n, n + 1:end) * Q0(n + 1:end)');
        end
        Q0 = Q;
    end

    DAMAS_result = Q;

end
