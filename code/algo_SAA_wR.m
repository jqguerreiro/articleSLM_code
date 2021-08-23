%% SIMULATED ANNEALING ALGORITHM WRN OPTIMIZATION
clear all; clc;
addpath(strcat(pwd,'\otslm-master'));
addpath(strcat(pwd,'\comp_funcs'));

%% Parameters for the SLM
side = 16;                                      %choose 8 (8x8) or 16 (16x16)
beam_size_in_SLM = 2.048;
off_set = [-165,60];

%Parameters for the SLM
    pixel_size = 8E-3;                          %mm
    slm_side_size = beam_size_in_SLM;           %mm
    slm_side_real = slm_side_size/pixel_size;   %number of pixels of slm side according to beam_size_in_SLM
    factor = ceil(slm_side_real/side);          %factor to use repelem in phase matrix
    slm_side = factor*side;                     %actual number of pixels of slm side according to side value chosen
    slm_side1 = slm_side*pixel_size;            %mm

slm_center = [1080,1920]./2 + off_set;          %Adjust center of slm matrix if needed
offset = round(slm_center - slm_side/2);

%Setting up the CMOS Camera
v = videoinput('avtmatlabadaptor64_r2009b', 1, 'F7M0_Mono8_752x480');
triggerconfig(v, 'manual');
w=warning("query","last"); warning("off",w.identifier); clc;

%Setting up the SLM
slm = otslm.utils.ScreenDevice(2, 'pattern_type', 'phase', ...
    'prescaledPatterns', false, ...
    'size', [slm_side,slm_side], ...
    'offset', offset, ...
    'value_range', {0:255});

%% Capture initial image
% Show the ScreenDevice figure
slm.showRaw();
w=warning("query","last"); warning("off",w.identifier); clc; pause(1);
start(v)
slm.showRaw(); pause(1);
I0 = getsnapshot(v);
stop(v)

%% Choose the weights (optional) or Load weights, below
figure; imagesc(I0); colorbar;
centroid = ginput(1);
radius = ginput(1);
close();

%Weighted Radius Weights
[X,Y]= meshgrid(1:size(I0,2),1:size(I0,1)); %X=cols; Y=rows
distances = sqrt((X-centroid(1)).^2 + (Y-centroid(2)).^2);
weights = (distances./max(distances(:))).*100;
rad = pdist2(centroid, radius);
center = round(centroid);
[xgrid, ygrid] = meshgrid(1:size(I0,2), 1:size(I0,1));
mask = ((xgrid-center(1)).^2 + (ygrid-center(2)).^2) <= rad.^2;
weights(mask) = 1;

%Plot the selected weights
Ishow = I0;
new_mask = ( (((xgrid-center(1)).^2 + (ygrid-center(2)).^2) <= (rad.^2 + 200)) & (((xgrid-center(1)).^2 + (ygrid-center(2)).^2) >= (rad.^2 - 200)) );
Ishow(new_mask) = uint8(255);
figure; imagesc(Ishow); colorbar; caxis([0 max(I0(:))]); hold on;
plot(centroid(1),centroid(2),'r+','MarkerSize',5,'LineWidth',1);

%% Save the weights
save("weights_wR.mat", "center", "centroid", "rad", "radius", "weights");

%% Load the weights 
load("weights_wR.mat");

%% SAA
%Show the ScreenDevice figure
slm.showRaw();
%Calculate cost function value
phase0 = zeros(side);
value0 = weightedRadiusN(I0, weights);

%Initial parameters
iters = side^2 * 10;
best_phase=phase0; best_value=value0; best_I=I0;
top_phase=phase0; top_value=value0; top_I=I0;
%Save values for plotting
it = 0; iter = it; values = value0;

%Required inputs
segments = side^2;
percentage = 0.5;
n_x = round(percentage*segments);
p = 120E-3;

%Specific inputs
alfa = 0.97; T0 = 0.1; k = 0;
if side == 8
    n_markov = 8;
elseif side == 16
    n_markov = 32;
end

tic

while it<iters
    
    k=k+1;
    Tk=(alfa^k)*T0;
    for i=1:n_markov
        it=it+1;                                          %Iteration number
        new_phase=best_phase;
        %Introduce random disturbance to phase
        idx=randperm(segments,n_x);
        new_phase(idx)=rand(1,n_x);
        %Send phase to SLM and acquire image
        slm.showRaw(repelem(uint8(new_phase.*255),factor,factor)); pause(p);
        new_I = getsnapshot(v);
        %Calculates the cost function value
        new_value = weightedRadiusN(new_I,weights);
        %Probability of acceptance
        if ~isnan(new_value)
            if new_value < best_value
                %update values
                best_value = new_value;
                best_phase = new_phase;
                if new_value < top_value
                    top_value = new_value;
                    top_phase = new_phase;
                    top_I = new_I;
                    disp("Improved!");
                end
                %Save values for plotting
                iter = [iter, it];
                values = [values, best_value];
            else
                nR=rand; P=P_Metropolis(best_value,new_value,Tk);
                if nR < P
                    %update values
                    best_value = new_value;
                    best_phase = new_phase;
                    %Save values for plotting
                    iter = [iter, it];
                    values = [values, best_value];
                end
            end
        end
        disp(it);
    end
    
end

t=toc;
stop(v);
slm.close();

%Save values for plotting
iter = [iter, it];
values = [values, best_value];

if top_value < value0
    disp("Optimization made")
else
    disp("No better results")
end

%Overwrite values
best_value = top_value;
best_I = top_I;
best_phase = top_phase;

%% Results
figure(1);
subplot(1,2,1); imagesc(I0); colorbar; caxis([0 max(best_I(:))]);
subplot(1,2,2); imagesc(best_I); colorbar; caxis([0 max(best_I(:))]);
figure(2); imshow(best_phase); colorbar;

%% Save the results
name = "test_SAA_wR_"+string(side)+"x"+string(side)+".mat";
save(name, "side", "beam_size_in_SLM", "off_set", "I0", "best_I", "phase0", ...
    "value0", "best_phase", "best_value", "iter", "values", ...
    "alfa", "T0", "n_markov", "t");