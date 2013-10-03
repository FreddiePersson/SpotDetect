function mask = getColony(IM)
%
% function IMbinary=getColony(IMfluor)
%
% Detects the cell colony in the image IM
% The output IMb is a logic image with
% IMb(i,j)=true for coordinates (i,j) inside the 
% colony and IMb(i,j)=false outside the colony.
%
% The debug argument should be removed
% before the function is implemented in the code. 
% It is temporary for you to see how well the code works.
% 

%% SETTINGS
% Using a debug mode that displays the colony segmentation
debug = true;

% The waveletplane that is used for detecting
% the cell colony
wplane=7;

%% ACTUAL CODE

IM=double(IM);

if debug
    figure
    imagesc(IM)
end

% Calculate the wavelet planes
WP = SMT_ATrous(IM, 'Levels', 7);
actWP = WP(:, :, wplane);

% Segment the colony with the Otsu threshold method
actWP_gray = mat2gray(actWP);
th = graythresh(actWP_gray);
actWP_th = actWP_gray>th;

% Fill the binary image. This is done because the
% binary image of the colony may contain holes.
regProp=regionprops(actWP_th,'FilledImage','BoundingBox');

filledMask = regProp.FilledImage;
boundBox = regProp.BoundingBox;

jj = ceil(boundBox(1));
ii = ceil(boundBox(2));

dj = boundBox(3)-1;
di = boundBox(4)-1;

mask = zeros(size(IM));

% Assign the filled binary bounding box image to a binary image
% with the same size as IM.
mask(ii:ii+di,jj:jj+dj) = filledMask;

if debug
    disp('When running many files and batches turn off the debug flag in getColony.m')
    figure;imagesc(IM.*~mask)
    figure;imagesc(IM.*mask)
    figure;imagesc(WP(:, :, 6).*mask)
    size(find(mask==0))
    size(find(WP(:, :, 2).*mask==0))
end

end

