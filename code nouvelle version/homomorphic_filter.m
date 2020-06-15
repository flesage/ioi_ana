function im2 = homomorphic_filter(im,sig_ph,sig_pb,alpha,beta)
I = double(im);
I = log(1 + I);

M = 2*size(I,1) + 1;
N = 2*size(I,2) + 1;

sigma = sig_ph;

[X, Y] = meshgrid(1:N,1:M);
centerX = ceil(N/2);
centerY = ceil(M/2);
gaussianNumerator = (X - centerX).^2 + (Y - centerY).^2;
H = exp(-gaussianNumerator./(2*sigma.^2));
H = 1 - H;

Hemphasis = alpha + beta*H;

If = fft2(I, M, N);
Iout = real(ifft2(Hemphasis.*If));
Iout = Iout(1:size(I,1),1:size(I,2));

im2 = uint16(imgaussfilt(((exp(Iout) - 1)),sig_pb));
end

