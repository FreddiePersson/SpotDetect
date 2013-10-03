function T = genAnscombe(I, varargin)
%% T=SMT_genAnscombe(I, meanGauss, stdGauss, gain)
%
% Implementation of both the traditional Anscombe transform and a variant of 
% the Anscombe variance stabilizing transform that takes into account Poisson 
% noise from the detector AND additative Gaussian read-out
% noise. The variance stabilizing formulas aim at transforming poisson
% distributions, where the variance is linked to the mean, to gauss
% distributions with a constant variance for easier denoising etc. 
% The generalized transform is taken from:
% F. Murtagh, J.-L. Starck, A. Bijaoui, "Image restoration with noise supression 
% using a multiresolution support", Astrology and Astrophysics, 112 (1995)
%
%
% Input: 
% I               : The data that should be operated on.
%
% Possible extra input:
% meanGauss       : The mean of the gaussian distributed noise/data
% stdGauss        : The std dev of the gaussian distributed noise/data
% gain            : the detector gain
%
% Output: 
% An array of the same size as the input data 'I'.
% 
% options:
%

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMT_genAnscombe.m, Performs a general version of the Anscombe transform 
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
I = double(I);

if nargin == 1 
    % Run the transform for 'pure' Poisson noise.
    T = 2*sqrt(I+3/8);
elseif nargin == 4
    % Run the transform for additative Poisson and Gaussian noise.
    meanGauss = varargin{1};
    stdGauss = varargin{2};
    gain = varargin{3};
    T = (2/gain)*sqrt(gain*double(I) + (3/8) + stdGauss^2-gain*meanGauss);

end




