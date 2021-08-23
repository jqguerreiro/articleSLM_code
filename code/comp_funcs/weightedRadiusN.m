function outputArg = weightedRadiusN(I,weights)
%Calculates the weighted radius
%   I=input Irradiance matrix from the CCD Camera (must be double image)
%   outputArg= ponderated radius value
%   Calculates de centroid position, distances of each pixel to it,
%   normalizes the input matrix,
%   multiplies the I_norm*distances and calculates the mean of it.

%% NORMALIZE I
new_I = im2double(I).*255;
%Maximum value of irradiance
M = max(new_I(:));
%Normalized Irradiance matrix
I_norm = new_I./M;
%% D_ponderada and outputArg
D_ponderada = I_norm.*weights;
outputArg = sqrt(mean2(D_ponderada(new_I>0)));
end