% Redas in a tif or STK file to variable rawData in the workspace

% Get filename and path with "uigetfile"
[filename, pathname] = uigetfile({'*.stk; *.tif', 'Images (*.tif,*.stk)'; '*.tif',  'TIF files (*.tif)'; '*.stk', 'STK files (*.stk)'; '*.*', 'All Files (*.*)'}, 'Select image file (STK or TIF)');

if ( filename == 0 )
    disp('Error! No (or wrong) file selected!')
    filename = 0;
    pathname = 0;
    return
end
full_filename = [ pathname, filename ];

if ~isempty(strfind(filename, '.tif'))
    
    % Calc. new variables
    imInfo = imfinfo(full_filename);
    imWidth = imInfo(1).Width;
    imHeight = imInfo(1).Height;
    stackSize = numel(imInfo);
    
    rawData = zeros(imHeight, imWidth, stackSize);
    
    t = Tiff(full_filename, 'r');
    for inFrame = 1:stackSize
        inFrame;
        rawData(:, :, inFrame) = t.read();
        if inFrame<stackSize
            t.nextDirectory();
        end
    end
    t.close();
    clear 'im*' 't' 'inFrame' 'filename' 'pathname' 'stackSize';
    
else
    [rawData, timeBetweenFrames] = readSTK(full_filename);
    stackSize = size(rawData, 3);
timeBetweenFrames;
end