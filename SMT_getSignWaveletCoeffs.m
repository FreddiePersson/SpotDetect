function WP = SMT_getSignWaveletCoeffs(I, noiseLevels, varargin)
%% resultImg = SMT_getSignWaveletCoeffs(I, noiseLevels, varargin)
%
% Gets the significant wavelet coefficients provided the given significance levels. 
% It can either use a hard threshold or the softer Jeffreys non-informative
% prior based threshold.
%
% Input: 
% I               : The data that should be operated on.
%
% Output: 
% An [m,n,levels], where [m, n]=size(I).
% 
% options:
% 'Levels'        : followed by the number of wavelet planes to use.
%                   Default value is 3.
% 'Jeffrey'       : if given, the softer threshold based on Jeffreys 
%                   non-informative prior will be used.
% 'noiseTH'     : followed by the threshold to use for the wavelet coeff 
%                   significance, sign = TH*noiselevels. Default is 3.
% 'actualPlaneNoise'   : if given the significance level will be estimated from the 
%                   actual plane, using a numeric prefactor (Olivo-Marin 2002). 
%                   This is the default method.
% 'firstPlaneNoise'    : if given the significance level will be estimated from
%                   the first wavelet plane. 
%
%

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMT_getSignWaveletCoeffs.m, Performs a general version of the Anscombe transform 
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

%% parse input
do_jeffrey = false;
extraArgs_AT = {''};
levels = 3;
noiseTH = 3;
do_actualPlaneNoise = true;
do_firstPlaneNoise = false;
 
if(nargin>2)        % parse options
    kmax = nargin-1;   % stopping criterion
    if iscell(varargin{1})
        varargin = varargin{1};     % fulhack
        kmax = length(varargin)+1;
    end
    % argument counter
    k = 1;
    while(k<kmax)
        option=varargin{k};
        if(strcmpi(option,'Jeffrey'))
            do_jeffrey = true;
            k=k+1;
        elseif(strcmpi(option,''))   % fulhack
            
            k=k+1;
        elseif(strcmpi(option,'Levels'))
            if(~isempty(varargin{k+1}))
                extraArgs_AT{end+1} = option;
                levels = varargin{k+1};
                extraArgs_AT{end+1} = levels;
                if(~isnumeric(levels) || levels<=2 || levels~=round(levels))
                    error('SMT_getSignWaveletCoeffs: Levels option must be followed by a positive integer larger than 2.')
                end
            end
            k=k+2; 
        elseif(strcmpi(option,'actualPlaneNoise'))
            do_actualPlaneNoise = true;
            do_firstPlaneNoise = false;
            k=k+1;    
        elseif(strcmpi(option,'firstPlaneNoise')) 
            do_firstPlaneNoise = true;
            do_actualPlaneNoise = false;
            k=k+1; 
        elseif(strcmpi(option,'noiseTH'))
            if(~isempty(varargin{k+1}))
                noiseTH = varargin{k+1};
                if(~isnumeric(noiseTH) || noiseTH<=0)
                    error('SMT_getSignWaveletCoeffs: noiseTH option must be followed by a positive number.')
                end
            end
            k=k+2; 
        else
            error(['SMT_getSignWaveletCoeffs: option ' option ' not recognized.'])
        end
    end
end
 
%% start of actual code

tempImg = zeros(size(I));
WP = SMT_ATrous(I, extraArgs_AT);

if do_actualPlaneNoise
    
    for ind = 1:levels
        waveletPlane = WP(:,:,ind);
        if do_jeffrey
            A = (WaveletPlane.^2-noiseTH*noiseLevels(ind)^2);    % Jeffreys non informative prior threshold
            A(A<0) = 0;
            tempImg = (1./WaveletPlane).*A;
        else
            tempImg(abs(waveletPlane) >= noiseTH*noiseLevels(ind)) = 1; % hard threshold
        end
        
        WP(:, :, ind) = waveletPlane.*tempImg;
    end
elseif do_firstPlaneNoise
    for ind = 1:levels
        waveletPlane = WP(:,:,ind);
        if do_jeffrey
            A = (WaveletPlane.^2-noiseTH*noiseLevels(1)^2);    % Jeffreys non informative prior threshold
            A(A<0) = 0;
            tempImg = (1./WaveletPlane).*A;
        else
            tempImg(abs(waveletPlane) >= noiseTH*noiseLevels(1)) = 1; % hard threshold
        end
        
        WP(:, :, ind) = waveletPlane.*tempImg;
    end
end
end

