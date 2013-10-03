
function [average sigma] = SMT_imgStats2D(I, kernelSize)

if mod(kernelSize+1, 2)
    error('SMT_imgStats2D: the kernel size should be an odd integer.');
end;

pad = (kernelSize-1)/2;
Img = padarray(I, [pad pad], 'replicate');
Img2 = Img.^2;

kernel = ones(kernelSize, kernelSize);


% Making use of the std dev for a sample: sigma = sqrt[(sum(x^2)-n*[sum(x)/n]^2)/(n-1)]
M = filter2(kernel, Img, 'valid');
M2 = filter2(kernel, Img2, 'valid');

sigma = sqrt((M2-(M.^2)*(1/kernelSize^2))/(kernelSize^2-1));

average = M/(kernelSize^2);       

end
