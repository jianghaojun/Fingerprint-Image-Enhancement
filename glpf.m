function H = glpf(D0,M)
% Create a Gaussian low pass filter

[DX, DY] = meshgrid(1:M);
D2 = (DX-M/2-1).^2+(DY-M/2-1).^2;
H = exp(-D2/(2*D0*D0));