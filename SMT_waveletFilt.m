function [signWP, multiSuppMask]=SMT_waveletFilt(I, varargin)
%% [filteredImg, multiSuppMask]=SMT_waveletFilt(I, varargin)
%
% Filters input data by applying a possibly iterative wavelet denoising approach using
% the A Trous wavelet transform. The iterative filtering is based on identifying
% significant wavelet coefficients and rebuilding the image from them. Then 
% take the residual (Image-reconstructed image) and do the same for that and
% add it to the old reconstructed image until the last reconstructed image
% only contains non significant data...
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
% An [m,n,levels+1], where [m, n]=size(I). The levels+1 image plane
% is the last smoothed image.
% 
% options:
% 'Levels'        : followed by the number of wavelet planes to calculate.
%                   Default value is 3.
% 'Iter'          : if given use iterative filtering on the image. 
% 'Display'       : if given, the filtered image and the multiscale support are displayed.
% 'maxIter'       : followed by the max number of iterations. Default is
%                   10.
% 'convLimit'     : The convergence criteria for the iterations. Default is
%                   0.001
% 'simNoise'      : if given, the significance level is determined from
%                   simulated gaussian noise (Starck et al 1995).
% 'actualPlaneNoise'   : if given the significance level will be estimated from the 
%                   actual plane, using a numeric prefactor (Olivo-Marin 2002). 
%                   This is the default method defined in SMT_getWaveletNoiseLevels.
% 'firstPlaneNoise'    : if given the significance level will be estimated from
%                   the first wavelet plane. 
% 'noiseTH'       : followed by the threshold to use for the wavelet coeff 
%                   significance, sign = TH*noiselevels. Default is 3.
% 'noiseROI'      : followed by a [y1, y2, x1, x2] array for the ROI. 
% 'autoDetectColony'    : automatically identifies a dense bacterial colony
%                   to use for noise estimation.
% 'Jeffrey'       : if given the Jeffrey non-informative prior will be used
%                   for selecting significant values instead of the common
%                   hard threshold.
%

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMT_waveletFilt.m, Performs iterative wavelet filtering
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

%% start of actual code

%% read options
levels = 3;
maxIter = 10;
do_display = false;
do_iter = false;
singlePlane = 0;
convLimit = 0.001;
extraArgs_gSWL = {''};
extraArgs_gSWC = {''};
extraArgs_AT = {''};


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
            extraArgs_gSWL{end+1} = option;
            extraArgs_AT{end+1} = option;
            k=k+1;
        elseif(strcmpi(option,''))   % fulhack
            k=k+1;
        elseif(strcmpi(option,'Levels'))
            if(~isempty(varargin{k+1}))
                extraArgs_gSWL{end+1} = option;
                extraArgs_gSWC{end+1} = option;
                extraArgs_AT{end+1} = option;
                levels = varargin{k+1};
                extraArgs_gSWL{end+1} = levels;
                extraArgs_gSWC{end+1} = levels;
                extraArgs_AT{end+1} = levels;
                if(~isnumeric(levels) || levels<=2 || levels~=round(levels))
                    error('SMT_waveletFilt: Levels option must be followed by a positive integer larger than 2.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'noiseROI'))
            if(~isempty(varargin{k+1}))
                extraArgs_gSWL{end+1} = option;
                roi = varargin{k+1};
                extraArgs_gSWL{end+1} = roi;
                if(~isnumeric(roi) | length(roi)~=4 | roi~=round(roi))
                    error('SMT_waveletFilt: noiseROI option must be followed by a positive integer array of length 4.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'Iter'))
                do_iter = true;
            k=k+1;  
        elseif(strcmpi(option,'maxIter'))
            if(~isempty(varargin{k+1}))
                maxIter=varargin{k+1};
                if(~isnumeric(maxIter) || maxIter<=0 || maxIter~=round(maxIter))
                    error('SMT_waveletFilt: maxIter option must be followed by a positive integer.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'convLimit'))
            if(~isempty(varargin{k+1}))
                convLimit=varargin{k+1};
                if(~isnumeric(convLimit) || convLimit<=0 || convLimit>=1)
                    error('SMT_waveletFilt: convLimit option must be followed by a positive number smaller than 1.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'noiseTH'))
            if(~isempty(varargin{k+1}))
                extraArgs_gSWC{end+1} = option;
                extraArgs_gSWC{end+1} = varargin{k+1};
                if(~isnumeric(varargin{k+1}) || varargin{k+1}<=0)
                    error('SMT_waveletFilt: noiseTH option must be followed by a positive number.')
                end
            end
            k=k+2;
        elseif(strcmpi(option,'Jeffrey'))
            extraArgs_gSWC{end+1} = option;
            k=k+1;
        elseif(strcmpi(option,'simNoise')) 
            extraArgs_gSWL{end+1} = option;
            k=k+1;
        elseif(strcmpi(option,'actualPlaneNoise'))
            extraArgs_gSWC{end+1} = option;
            k=k+1;    
        elseif(strcmpi(option,'firstPlaneNoise')) 
            extraArgs_gSWC{end+1} = option;
            k=k+1;  
        elseif(strcmpi(option,'autoDetectColony')) 
            extraArgs_gSWL{end+1} = option;
            k=k+1;  
        else
            error(['SMT_waveletFilt: option ' option ' not recognized.'])
        end
    end
end
   
%% modify the input data
I_orig = double(I);
originalImg = I_orig;

%% Filtering using significant wavelet coefficients

noiseLevels = SMT_getWaveletNoiseLevels(originalImg, extraArgs_gSWL);  
signWP = SMT_getSignWaveletCoeffs(originalImg, noiseLevels, extraArgs_gSWC);

if do_iter
    % iterative filtering
    %[Starck et al, "Multiresolution support applied to image filtering and restoration", (1995)]
    signImg = sum(signWP(:, :, 1:levels+1), 3);
    residualImg = originalImg - signImg;
    sigma_res0 = std2(residualImg);
    
    epsilon = 1;
    count = 1;
    while (epsilon > convLimit && count < maxIter)
        signWP = SMT_getSignWaveletCoeffs(residualImg, noiseLevels, extraArgs_gSWC);
        signRes = sum(signWP(:, :, 1:levels+1), 3);
        signImg = signImg + signRes; % add significant residuals
        residualImg = originalImg - signImg;
        sigma_res1 = std2(residualImg);
        epsilon = abs(sigma_res0/sigma_res1 - 1);
        sigma_res0 = sigma_res1;
        count = count+1;
    end
    signWP = SMT_ATrous(signImg, extraArgs_AT);
end
    
%% Multiscale product of wavelets
% The multiresolution support is given by the product in the wavelet domain.
multiSuppImg = abs(prod(signWP(:, :, 1:levels), 3));
multiSuppMask = zeros(size(multiSuppImg));
multiSuppMask(multiSuppImg>std2(multiSuppImg))=1;
end

%%%%%%%%%%%%%%%%
% Subfunctions %
%%%%%%%%%%%%%%%%

