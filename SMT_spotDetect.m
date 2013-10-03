function [mask, n]=SMT_spotDetect(I, varargin)
%% [mask, n]=SMT_spotDetect(I, varargin)
%
% Detects local intensity regions, dots, by using selected wavelet planes or 
% wavelet multiscale products and iterative filtering with significant coefficients. 
% In addition, various morphological operations can be performed to clean up the resulting data. 
% For more reading consult:
% J.-L. Starck, F. Murtagh, A. Bijaoui, "Multiresolution Support Applied to 
% Image Filtering and Restoration", Graph. Models Image Process, 57:5 (1995)
% J.-C. Olivo-Marin, "Extraction of spots in biological images using 
% multiscale products", Pattern Recognition, 35 (2002).
%
% Input: 
% I               : The data that should be operated on.
%
% Output: 
% A struct containing information on the detected spots (spots) and the filtered 
% data used to find the spots.
% 
% options:
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
%                   This is the default method defined in SMT_getWaveletNoiseLevels.
% 'firstPlaneNoise'    : if given the significance level will be estimated from
%                   the first wavelet plane. 
% 'noiseROI'      : only use a subimage to estimate the noise. 
% 'autoDetectColony'    : as noiseROI but automatically identifies a dense
%                   bacterial colony.
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
% 
%

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMT_spotDetect.m, Detects dots by wavelet filtering.
% =========================================================================
% 
% Copyright (C) 2012 Fredrik Persson
% 
% E-mail: freddie.persson@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.


%% set and read options
levels = 3;
% neighbourhood = 9;
plane = 2;
do_display = false;
do_multiscale = false;
do_iter = false;
extraArgs_wF = {''};

if(nargin>1)        % parse options
    kmax = nargin;   % stopping criterion
    if iscell(varargin{1})
        varargin = varargin{1};     % fulhack
        kmax = length(varargin)+1;
    end
    % argument counter
    k = 1; 
    while(k<kmax)
        option=varargin{k};
        if(strcmpi(option,'Display'))
            do_display = true;
            extraArgs_wF{end+1} = option;
            k=k+1;
        elseif(strcmpi(option,'multiScale'))
            do_multiscale = true;
            k=k+1;
        elseif(strcmpi(option,'Iter'))
            do_iter = true;
            k=k+1;
        elseif(strcmpi(option,''))   % fulhack

            k=k+1;
        elseif(strcmpi(option,'Levels'))
            if(~isempty(varargin{k+1}))
                extraArgs_wF{end+1} = option;
                levels=varargin{k+1};
                extraArgs_wF{end+1} = levels;
                if(~isnumeric(levels) || levels<=2 || levels~=round(levels))
                    error('SMT_spotDetect: Levels option must be followed by a positive integer larger than 2.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'noiseROI'))
            if(~isempty(varargin{k+1}))
                extraArgs_wF{end+1} = option;
                roi = varargin{k+1};
                extraArgs_wF{end+1} = roi;
                if(~isnumeric(roi) | length(roi)~=4 | roi~=round(roi))
                    error('SMT_spotDetect: noiseROI option must be followed by a positive integer array of length 4.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'Plane'))
            if(~isempty(varargin{k+1}))
                plane = varargin{k+1};
                if(~isnumeric(plane) || plane<0 || plane~=round(plane) || plane>levels)
                    error('SMT_spotDetect: Plane option must be followed by a positive integer smaller than levels.')
                end
            end
            k=k+2;  
%         elseif(strcmpi(option,'neighbourhood'))
%             if(~isempty(varargin{k+1}))
%                 neighbourhood=varargin{k+1};
%                 if(~isnumeric(neighbourhood) || neighbourhood<0 || neighbourhood~=round(neighbourhood) || (mod(neighbourhood+1, 2)&neighbourhood~= 0))
%                     error('SMT_spotDetect: neighbourhood option must be followed by a positive odd integer.')
%                 end
%             end
%             k=k+2;  
        elseif(strcmpi(option,'maxIter'))
            if(~isempty(varargin{k+1}))
                extraArgs_wF{end+1} = option;
                extraArgs_wF{end+1} = varargin{k+1};
                if(~isnumeric(varargin{k+1}) || varargin{k+1}<=0 || varargin{k+1}~=round(varargin{k+1}))
                    error('SMT_spotDetect: maxIter option must be followed by a positive integer.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'convLimit'))
            if(~isempty(varargin{k+1}))
                extraArgs_wF{end+1} = option;
                extraArgs_wF{end+1} = varargin{k+1};
                if(~isnumeric(varargin{k+1}) || varargin{k+1}<=0 || varargin{k+1}>=1)
                    error('SMT_spotDetect: convLimit option must be followed by a positive number smaller than 1.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'noiseTH'))
            if(~isempty(varargin{k+1}))
                extraArgs_wF{end+1} = option;
                extraArgs_wF{end+1} = varargin{k+1};
                if(~isnumeric(varargin{k+1}) || varargin{k+1}<=0)
                    error('SMT_spotDetect: noiseTH option must be followed by a positive number.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'simNoise')) 
            extraArgs_wF{end+1} = option;
            k=k+1;
        elseif(strcmpi(option,'actualPlaneNoise'))
            extraArgs_wF{end+1} = option;
            k=k+1;    
        elseif(strcmpi(option,'firstPlaneNoise')) 
            extraArgs_wF{end+1} = option;
            k=k+1;  
        elseif(strcmpi(option,'Jeffrey'))
            extraArgs_wF{end+1} = option;
            k=k+1;
        elseif(strcmpi(option,'autoDetectColony'))
            extraArgs_wF{end+1} = option;
            k=k+1;
        else
            error(['SMT_spotDetect: option ' option ' not recognized.'])
        end
    end
end
   
%% initiate input and output data
I_orig = double(I);

% Run a variace stabilizing transform to account for most noise being
% Poisson distributed. Can be run also for addidative Poisson and Gaussian noise, 
% but then the gauss distributiona and detector gain has to be known.
originalImg = SMT_genAnscombe(I_orig);


%%
% Send on to wavelet filtering
[filtWP, multiSuppMask] = SMT_waveletFilt(originalImg, extraArgs_wF);

%%

if plane ==0
    filtImg = sum(filtWP(:, :, 2:levels), 3);
else
    filtImg = filtWP(:, :, plane);
end
if do_display
    if plane ==0
        figure('Name', ['Sign coeff of the summed wavelet planes (2:' num2str(levels) ')' ]);
        filtImg = sum(filtWP(:, :, 2:levels), 3);
    else
        figure('Name', ['Sign coeff of the selected wavelet plane (' num2str(plane) ')']);
        filtImg = filtWP(:, :, plane);
    end
    imagesc(filtImg);
    if do_multiscale
        figure('Name', ['Multiscale support mask']);
        imagesc(multiSuppMask);
    end
end


%% create a binary mask

mask = zeros(size(originalImg));
mask(filtImg>0) = 1;

% if neighbourhood>0
%     % use neighbourhood filtering (local image statistics)
%     [imgAvg, imgStd] = SMT_imgStats2D(originalImg, neighbourhood);
% else
%     % use global statistics
%     imgAvg = mean(originalImg(:))*ones(size(originalImg));
%     imgStd = std2(originalImg)*ones(size(originalImg));
% end
%
% imgTH = imgAvg + noiseTH*imgStd;
% % make mask based on noise threshold alone
% mask(filtImg>imgTH) = 1;

if do_multiscale
    mask = mask.*multiSuppMask;
end

if do_display
    figure('Name', ['Masked raw image, prior to morphological processing']);
    imagesc(I_orig.*mask);
end

%% clean up the mask

% % mask = bwmorph(mask, 'fill');
mask = bwmorph(mask, 'clean');
mask = bwmorph(mask, 'hbreak');
mask = bwmorph(mask, 'spur');
mask = bwmorph(mask, 'clean');
mask = bwmorph(mask, 'thicken');

% extra steps including an erode step
if 0
    mask = bwmorph(mask, 'erode');
    mask = bwmorph(mask, 'hbreak');
    mask = bwmorph(mask, 'spur');
    mask = bwmorph(mask, 'clean');
    mask = bwmorph(mask, 'thicken');
end

if do_display
    figure('Name', ['Masked image, after morphological processing']);
    imagesc(filtImg.*mask);
end

[~, n] = bwlabel(mask);


end
