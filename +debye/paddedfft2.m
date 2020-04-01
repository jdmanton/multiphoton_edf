function [matrix_fft2, matrix_fft] = paddedfft2(M, N, matrix)

% Pad and transform in x, then pad and transform in y for efficiency over padding both ways at the start

fft_x = [zeros((N - 2*M - 1) / 2, 2*M + 1); matrix; zeros((N - 2*M - 1) / 2, 2*M + 1)];
matrix_fft = fftshift(fft(fftshift(fft_x, 1)), 1);
matrix_fft = matrix_fft((N - 2*M -1) / 2 + 1:(N - 2*M -1) / 2 + 1 + 2*M, :);

matrix_fft2 = [zeros(2*M + 1, (N - 2*M -1) / 2), matrix_fft, zeros(2*M + 1,(N - 2*M -1) / 2)];
matrix_fft2 = fftshift(fft(fftshift(matrix_fft2.', 1)), 1).';

matrix_fft2 = matrix_fft2(:, (N - 2*M - 1) / 2 + 1:(N - 2*M - 1) / 2 + 1 + 2*M);

end
