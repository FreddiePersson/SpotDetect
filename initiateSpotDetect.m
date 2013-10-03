function [coords] = SMT_initiateSpotDetect(IP, noiseTh)

% Function to run SMT_spotDetect on a stack IP(:, :, x) with x planes. All
% used options are set here in the code except for the Noise threshold which
% is an input parameter. Comments below is valid for cell population images 
%
% When using it with actualPlaneNoise it is sensitive to the amount of
% black (no cells) in each frame. Less total cellarea means more detected dots
% (lower effective threshold), and the opposite.
% Cropped images just containing cells use noiseTH = 3. 
%
% The best approach is to detect the areas without cells first and dont use 
% that in the noiseTH calculation (carried out in SMT_getWaveletNoiseLevels.m).
%
% Using it with firstPlaneNoise use noiseTH = 3. It is now sensitive to
% changing SNR instead in the same manner. Just as stable wave detection (SWD).
%
% Possible values to send to the spotDetector:
% 'Levels'        : followed by the number of wavelet planes to calculate.
%                   Default value is 3.
% 'Plane'         : followed by an integer saying what wavelet plane to
%                   detect in. If 0 it uses a sum of wavelet planes 2:levels.
%                   Default is 2.
% 'noiseTH'       : followed by the threshold to use for the wavelet coeff 
%                   significance, sign = noiseTH*noiselevels. Default is 3.
% 'simNoise'      : if given, the significance level is determined from
%                   simulated gaussian noise (Starck et al 1995).
% 'actualPlaneNoise'   : if given the significance level will be estimated from the 
%                   actual plane, using a numeric prefactor (Olivo-Marin 2002). 
%                   This is the default method.
% 'firstPlaneNoise'    : if given the significance level will be estimated from
%                   the first wavelet plane. 
% 'noiseROI'      : only use a subimage to estimate the noise. 
% 'autoDetectColony'    : same as noiseROI but detects a dense bacterial 
%                   cell colony to use automatically using the 7th wavelet plane. 
% 'noiseROI'      : only use a subimage to estimate the noise. 
% 'Jeffrey'       : if given the Jeffrey non-informative prior will be used
%                   for selecting significant values instead of the common
%                   hard threshold.
% 'Display'       : if given, displays a varity of images through different
%                   stages. Can lead to many open figures.
% 'maxIter'       : followed by the max number of iterations for the iterative 
%                   significant wavelet coefficient filtering. Default is
%                   10.
% 'convLimit'     : The convergence criteria for the iterative significant 
%                   wavelet coefficient filtering. Default is 0.001.
% 'Iter'          : if given use iterative filtering on the image.  
% 'multiScale'    : if given the multiscale product in wavelet space will
%                   be used.
% 'neighbourhood' : COMMENTED OUT!!! 
%                   followed by the size (one odd integer) of the square 
%                   kernel around the pixel to be considered for thresholds. 
%                   Set to 0 it the whole image statistics is to be used.
%                   Default is 9.


display = false;
do_noiseROI = false;
do_autoDetectColony = true;
plane = 2;

%%
coords = [];

% Mark out the noise ROI in the first image
if do_noiseROI
    figure; imagesc(rawData(:, :, 1));
    rect = getrect;
    noiseROI = [ceil(rect(2)), floor(rect(2)+rect(4)), ceil(rect(1)), floor(rect(1)+rect(3))];
%     figure; imagesc(rawData(noiseROI(1):noiseROI(2), noiseROI(3):noiseROI(4), 1));
end


for ind=1:size(IP, 3); 
    if do_noiseROI & do_autoDetectColony
        error('initiateSpotDetect: You can only use one way of selecting ROI.')
    elseif do_noiseROI
        [mask, n]=SMT_spotDetect(IP(:, :, ind), 'Plane', plane, 'actualPlaneNoise', 'noiseTH', noiseTh, 'noiseROI', noiseROI);
    elseif do_autoDetectColony
        [mask, n]=SMT_spotDetect(IP(:, :, ind), 'Plane', plane, 'actualPlaneNoise', 'noiseTH', noiseTh, 'autoDetectColony');
    else
        [mask, n]=SMT_spotDetect(IP(:, :, ind), 'Plane', plane, 'actualPlaneNoise', 'noiseTH', noiseTh);
    end
    
    stats = regionprops(mask, IP(:, :, ind), 'WeightedCentroid');
    for ind2=1:length(stats)
        coord=stats(ind2).WeightedCentroid;
        coords = [coords; coord, ind];
    end
    
    if display %&& ind == 24;
        figure;
        hold on
        imagesc(IP(:, :, ind));
        colormap(gray)
        for ind2=1:length(stats)
            hold on
            coord=stats(ind2).WeightedCentroid;
            circle(coord(1), coord(2), 5);
        end
        hold off
    end; 
    
end;
