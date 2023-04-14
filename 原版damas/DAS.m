function [DAS_result, e, CSM] = DAS(N, z0, f, rn, source, SNR)
    N_mic = size(rn, 1);
    c = 343;
    omega = 2 * pi * f;
    DAS_result = zeros(N, N);
    scan_range_X = linspace(-2, 2, N);
    scan_range_Y = linspace(2, -2, N);
    [X, Y] = meshgrid(scan_range_X, scan_range_Y);

    d = zeros(N, N, N_mic);
    e = zeros(N, N, N_mic);
    S = zeros(N, N, N_mic);

    d0 = sqrt(X .^ 2 + Y .^ 2 + z0 ^ 2);

    for n = 1:N_mic
        d(:, :, n) = sqrt((X - rn(n, 1)) .^ 2 + (Y - rn(n, 2)) .^ 2 + z0 ^ 2);
        e(:, :, n) = (d(:, :, n) ./ d0) .* exp(-1j * omega .* d(:, :, n) ./ c);
        S(:, :, n) = (d0 ./ d(:, :, n)) .* exp(-1j * omega .* d(:, :, n) ./ c) + 10 ^ (-SNR / 10) * (rand(N, N) + 1j * rand(N, N));
    end

    CSM = zeros(N_mic, N_mic);

    for k = 1:size(source, 1)
        disp("size of S");
        disp(size(squeeze(S(source(k, 2), source(k, 1), :))))
        CSM = CSM + squeeze(S(source(k, 2), source(k, 1), :)) * squeeze(S(source(k, 2), source(k, 1), :))';
    end

    for index1 = 1:length(X)

        for index2 = 1:length(Y)
            a = squeeze(e(index1, index2, :));
            DAS_result(index1, index2) = dot(a, CSM * a) / N_mic ^ 2;
        end

    end

end
