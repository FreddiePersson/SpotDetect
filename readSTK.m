
function [rawData, timeBetwFrames, z] = readSTK(varargin)

if nargin == 0
    % Read in the filename and path
    [filename, pathname] = uigetfile({'*.stk'}, 'Select image stack file (.stk)');
    if ( filename == 0 )
        disp('Error! No (or wrong) file selected!')
        return
    end
    full_filename = [ pathname, filename ];
else
    full_filename = varargin{1};
end

% Open the file in read mode
fid = fopen(full_filename, 'r', 'l');

% Jump in 4 bytes in the file
fseek(fid, 4, 'bof');

% Get offset to the first image directory
firstImDirOffset = fread(fid, 1, 'uint32');

% Go to the first image directory
fseek(fid, firstImDirOffset, 'bof');

% Read in the number of directory fields is in each directory.
numDirFields = fread(fid, 1, 'uint16');

% Loop through the directory fields to find the height, width, # images and
% strip offset.
for ind = 0:numDirFields
   
    % Read the field tag
    fieldTag = fread(fid, 1, 'uint16');
    % Read the field type
    fieldType = fread(fid, 1, 'uint16');
    % Read the number of values in the field
    numVal = fread(fid, 1, 'uint32');
    % Read the values in the field
    val = fread(fid, 1, 'uint32');
    
    % Assign the values according to the field tags.
    switch fieldTag
        
        case 256
            imWidth = val;
        case 257
            imHeight = val;
        case 33629
            numImages = numVal;
            fPos = ftell(fid);
            % Go to the list and the 3rd field 
            fseek(fid, val+3*4, 'bof');
            % Read in the times in microseconds for all frames
            times = fread(fid, numVal, 'long', 5*4); % Skips 5 fields
            timeBetwFrames = (times(end)-times(1))/numVal;
            fseek(fid, fPos, 'bof');
    end
end

% Manually go and read in the z positions of the frames
fPos = ftell(fid);
fseek(fid, 0, 'eof');
fSize = ftell(fid);
fseek(fid, fPos, 'bof');
fzPos = fSize-numImages*12-286;
fseek(fid, fzPos, 'bof');
zNum = fread(fid, numImages, 'long', 4);
fseek(fid, fzPos+4, 'bof');
zDen = fread(fid, numImages, 'long', 4);
z = double(zNum./zDen)*1000;   % in nm

% Go to the position of the first image frame in the file.
fseek(fid, 8, 'bof');

% Read in the images since in .stk the images are just stored directly after each other its just to read in all the data.
rawData = fread(fid, numImages*imWidth*imHeight, 'uint16');

fclose(fid);

% Turn and change the matrix so its correctly displayed as image in Matlab
rawData = permute(reshape(rawData, imWidth, imHeight, numImages), [2, 1, 3]);
rawData = flipdim(rawData, 1);
end