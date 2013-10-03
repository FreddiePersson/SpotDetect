function WP=SMT_ATrous(I, varargin)
%% WP=SMT_ATrous(I, varargin)
%
% Applies the A Trous wavelet transform and gives the 'levels' first wavelets 
% from the image I. This wavelet acts by applying a smoothing filter (often gaussian)
% on the image and the 1st wavelet plane is the difference the original image
% and the first smoothed image. Then the smoothed image is filtered by an even 
% coarser smoothing filter etc. 
% For more reading consult:
% J.-L. Starck, F. Murtagh, A. Bijaoui, "Multiresolution Support Applied to 
% Image Filtering and Restoration", Graph. Models Image Process, 57:5 (1995)
% J.-C. Olivo-M-Marin, "Extraction of spots in biological images using 
% multiscale products", Pattern Recognition, 35 (2002).
%
% Input: 
% I               : The image that should be operated on.
%
% Output: 
% An [m,n,levels+1], where [m, n]=size(I). The levels+1 image plane
% is the last smoothed image.
% 
% options:
% 'Levels'        : followed by the number of wavelet planes to calculate.
%                   Default value is 3.
% 'Display'       : if given, the wavelets are displayed in individual plots.
%

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMT_ATrous.m, Performs the A Trous wavelet transform
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



%% read options
levels = 3;
do_display = false;

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
            k=k+1;
        elseif(strcmpi(option,''))     % fulhack
            
            k=k+1;
        elseif(strcmpi(option,'Levels'))
            if(~isempty(varargin{k+1}))
                levels=varargin{k+1};
                if(~isnumeric(levels) || levels<=2 || levels~=round(levels))
                    error('SMT_ATrous: Levels option must be followed by a positive integer larger than 2.')
                end
            end
            k=k+2;
        else
            error(['SMT_ATrous: option ' option ' not recognized.'])
        end
    end
end
    
%% define variables

[m, n] = size(I);


originalImg = I;
lastImg = I;
img = I;
WP = zeros(m, n, levels + 1);
baseKernel = (1/16)*[1 4 6 4 1];
kernel = [];


%% Make the wavelet planes

if do_display
    figure('Name', ['Wavelet planes']);
end

for ind=1:levels
    if ind>1
        kernel = zeros(1,2*length(kernel)-1);
        kernel(1:2^(ind-1):end) = baseKernel;
    else
        kernel = baseKernel;
    end
    
    % Pad the matrix before to avoid zero-padding
    pad = (length(kernel)-1)/2;
    lastImgPad = padarray(lastImg, [pad, pad], 'replicate');
    
    img = conv2(kernel, kernel, lastImgPad, 'valid');
    
    waveletPlane = double(lastImg)-double(img);
    
    if do_display
        subplot(1, levels, ind);
        imagesc(waveletPlane);
    end
    
    WP(:, :, ind) = waveletPlane;
    lastImg = img;
end
WP(:, :, end) = img;
end




