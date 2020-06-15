function [fcMapBlend hFig] = pat_overlay_blend(anatomical, fcMap, varargin)
% Blends a top layer (fcMap) on the base layer (anatomical).
% The parts of the top layer where base layer is light become lighter, the parts
% where the base layer is dark become darker. Depending on the value of the base
% layer, one gets a linear interpolation between black, the top layer and white.
% SYNTAX
% [fcMapBlend h] = pat_overlay_blend(anatomical, fcMap, brainMask, fcMapRange,...
%                                   alphaRange, fcColorMap, figIntensity)
% INPUTS
% anatomical    base layer, usually an anatomical image
% fcMap         top layer, usually a functional image
% [brainMask]   binary image, where pixels to display are marked with 1's.
%               Default: Show all pixels
% [fcMapRange]  2-element vector with lower and upper limits of range to display
%               Default: min & max of top layer
% [alphaRange]  transparency range, it can be a 2-element vector with lower and
%               upper limits of range or a 4-element vector for separate ranges.
%               Default: min & max of top layer
% [fcColorMap]  Colormap to use for the top layer, default: jet(256)
% [figIntensity] Weight of the base layer, between 0 and 1.
% OUTPUTS
% fcMapBlend    RGB overlay blend
% hFig          Handle to the figure
%_______________________________________________________________________________
% Copyright (C) 2013 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

%% Handling optional inputs
% Colormap levels
nColorLevels = 256;
% only want 5 optional inputs at most
numVarArgs = length(varargin);
if numVarArgs > 5
    error('pat12:pat_overlay_blend:TooManyInputs', ...
        'requires at most 5 optional inputs: brainMask, fcMapRange, alphaRange, fcColorMap, figIntensity');
end
% set defaults for optional inputs
optArgs = {ones(size(anatomical)) [min(fcMap(:)) max(fcMap(:))] ...
    [min(fcMap(:)) max(fcMap(:))] jet(nColorLevels) 0.7};
% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% now put these defaults into the optArgs cell array, and overwrite the ones
% specified in varargin.
optArgs(newVals) = varargin(newVals);
% Place optional args in memorable variable names
[brainMask fcMapRange alphaRange fcColorMap figIntensity] = optArgs{:};

% Correct alpha range
if numel(alphaRange) ~= 2 && numel(alphaRange) ~= 4
    alphaRange = [min(fcMap(:)) max(fcMap(:))];
    fprintf('alphaRange = [%0.4f %0.4f]\n', alphaRange(1), alphaRange(2));
end

%% Prepare anatomical and functional images
% Convert anatomical image to grayscale (weighted by figIntensity)
anatomicalGray      = figIntensity .* mat2gray(anatomical);
anatomicalGray      = repmat(anatomicalGray,[1 1 3]);
% Convert functional image to RGB
fcMapGray           = mat2gray(fcMap, fcMapRange); % Fix range for correlation maps
fcMapIdx            = gray2ind(fcMapGray, nColorLevels);
fcMapRGB            = ind2rgb(fcMapIdx, fcColorMap);
% Set transparency according to mask and pixels range
pixelMask = false(size(brainMask));
% Single range (2-element vector)
if numel(alphaRange) == 2,
    pixelMask(fcMap > alphaRange(1) & fcMap < alphaRange(2)) = true;
elseif numel(alphaRange) == 4,
    % Double range (4-element vector)
    pixelMask(fcMap > alphaRange(1) & fcMap < alphaRange(2)) = true;
    pixelMask(fcMap > alphaRange(3) & fcMap < alphaRange(4)) = true;
end
fcMapRGB(repmat(~brainMask | ~pixelMask,[1 1 3])) = 0.5;

%% Apply overlay blend algorithm
fcMapBlend = 1 - 2.*(1 - anatomicalGray).*(1 - fcMapRGB);
fcMapBlend(anatomicalGray<0.5) = 2.*anatomicalGray(anatomicalGray<0.5).*fcMapRGB(anatomicalGray<0.5);

%% Display new figure
h = imshow(fcMapBlend, 'InitialMagnification', 'fit', 'border', 'tight');
hFig = gcf;
set(hFig, 'color', 'k')
% Allow printing of black background
set(hFig, 'InvertHardcopy', 'off');

% EOF