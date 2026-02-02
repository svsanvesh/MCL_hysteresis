function real_psi = lowpass_fil(p,cut_off)
psi_k = fft2(p);
cutoff_freq = cut_off;
% Create a circular low-pass filter mask
[h, w] = size(p);
[X1, Y1] = meshgrid(1:w, 1:h);
center_x = floor(w / 2) + 1;
center_y = floor(h / 2) + 1;
distance = sqrt((X1 - center_x).^2 + (Y1 - center_y).^2);
mask = double(distance >= cutoff_freq);
% figure;
% imagesc(mask)
% Apply the mask in the frequency domain
psi_k_masked = psi_k.* mask;

% Perform the inverse FFT
real_psi = real(ifft2(psi_k_masked));
end