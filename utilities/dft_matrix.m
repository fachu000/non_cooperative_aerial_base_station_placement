function U = dft_matrix(N)
% unitary N x N DFT matrix

U = 1/sqrt(N)*fft(eye(N),N,1);

end