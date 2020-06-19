%Generate grid of frequency values for an image of resolution resX,resY.
% Outputs:
% -wx : Array of horizontal frequency values.
% -wy : Array of vertical frequency values.
% -xC : Horizontal index of the zero frequency ( wx(:,xC) =0 / xAxis oriented to the right).
% -yC : Vertical index of the zero frequency ( wy(yC,:) = 0 / yAxis oriented upwards).

function [wx,wy,xC,yC] = GenerateFrequencyGrid(resX,resY)
    [wx,wy] = meshgrid(1:resX,1:resY);
    xC = ceil((resX+1)/2);
    yC = ceil((resY+1)/2);
    wx=double((wx(:,1:end)-xC)/(resX-1));
    wy=double((yC-wy(:,1:end))/(resY-1));
end