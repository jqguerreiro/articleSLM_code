%% FOUR-ELEMENT DIVISION ALGORITHM WA OPTIMIZATION
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

%Weighted Average Weights
[X,Y]= meshgrid(1:size(I0,2),1:size(I0,1)); %X=cols; Y=rows
distances = sqrt((X-centroid(1)).^2 + (Y-centroid(2)).^2);
weights = abs(distances - max(distances(:)));
weights = (weights./max(weights(:))).*100;
rad = pdist2(centroid, radius);
center = round(centroid);
[xgrid, ygrid] = meshgrid(1:size(I0,2), 1:size(I0,1));
mask = ((xgrid-center(1)).^2 + (ygrid-center(2)).^2) <= rad.^2;
weights(mask) = max(weights(:));

%Plot the selected weights
Ishow = I0;
new_mask = ( (((xgrid-center(1)).^2 + (ygrid-center(2)).^2) <= (rad.^2 + 200)) & (((xgrid-center(1)).^2 + (ygrid-center(2)).^2) >= (rad.^2 - 200)) );
Ishow(new_mask) = uint8(255);
figure; imagesc(Ishow); colorbar; caxis([0 max(I0(:))]); hold on;
plot(centroid(1),centroid(2),'r+','MarkerSize',5,'LineWidth',1);

%% Save the weights
save("weights_wA.mat", "center", "centroid", "rad", "radius", "weights");

%% Load the weights 
load("weights_wA.mat");

%% FEDA
%Show the ScreenDevice figure
slm.showRaw();
%Calculate cost function value
phase0 = zeros(side);
value0 = weightedAverage(I0, weights);

%Initial parameters
best_phase=phase0; best_value=value0; best_I=I0;
%Save values for plotting
it = 0; iter = it; values = value0;

%Required inputs
segments = side^2;
percentage = 0.5;
n_x = round(percentage*segments);
p = 120E-3;

tic

for n=1:log2(side)
    
    if n==1
        [r,c]=fedaCore(1:side,1:side);
        for i=1:size(r,1)
            for k=0:(1/10):0.9
               it=it+1;                                   %Iteration number
               new_phase=best_phase;
               new_phase(r(i,:),c(i,:))=k;
               %Send phase to SLM and acquire image
               slm.showRaw(repelem(uint8(new_phase.*255),factor,factor)); pause(p);
               new_I = getsnapshot(v);
               %Calculates the cost function value
               new_value = weightedAverage(new_I,weights);
               %Criteria of acceptance
               if ~isnan(new_value) && new_value>best_value
                   %Update values
                   best_phase = new_phase;
                   best_value = new_value;
                   best_I = new_I;
                   %Save values for plotting
                   iter = [iter, it]; 
                   values = [values, best_value];
                   disp("Improved!");
               end
            end
            if i==1
            w=warning("query","last"); warning("off",w.identifier);
            end
        end
    else
        for i=1:size(r,1)
            [r_new,c_new] = fedaCore(r(i,:),c(i,:));
            for m=1:size(r_new,1)
               for k=0:(1/10):0.9
                   it=it+1;                               %Iteration number
                   new_phase=best_phase;
                   new_phase(r_new(m,:),c_new(m,:))=k;
                   %Send phase to SLM and acquire image
                   slm.showRaw(repelem(uint8(new_phase.*255),factor,factor)); pause(p);
                   new_I = getsnapshot(v);
                   %Calculates the cost function value
                   new_value = weightedAverage(new_I,weights);
                   %Criteria of acceptance
                   if ~isnan(new_value) && new_value>best_value
                       %Update values
                       best_phase = new_phase;
                       best_value = new_value;
                       best_I = new_I;
                       %Save values for plotting
                       iter = [iter, it]; 
                       values = [values, best_value];
                       disp("Improved!");
                   end
               end
            end
            if i==1
                r_save = r_new;
                c_save = c_new;
            elseif i==size(r,1)
                r = [r_save;r_new];
                c = [c_save;c_new];
            else
                r_save = [r_save;r_new];
                c_save = [c_save;c_new];
            end
        end
    end
    
end

t=toc;
stop(v);
slm.close();

%Save values for plotting
iter = [iter, it];
values = [values, best_value];

if best_value > value0
    disp("Optimization made")
else
    disp("No better results")
end

%% Results
figure(1);
subplot(1,2,1); imagesc(I0); colorbar; caxis([0 max(best_I(:))]);
subplot(1,2,2); imagesc(best_I); colorbar; caxis([0 max(best_I(:))]);
figure(2); imshow(best_phase); colorbar;

%% Save the results
name = "test_FEDA_wA_"+string(side)+"x"+string(side)+".mat";
save(name, "side", "beam_size_in_SLM", "off_set", "I0", "best_I", "phase0", ...
    "value0", "best_phase", "best_value", "iter", "values", "t");