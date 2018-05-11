function score=obj_func(x, mask, penalty, columnsInImage, rowsInImage)

centerX=x(1);
centerY=x(2);
radius=x(3);

% centerX=round(x(1));
% centerY=round(x(2));
% radius=round(x(3));

[imageSizeX, imageSizeY] = size(mask);

% [columnsInImage, rowsInImage ] = meshgrid(1:imageSizeY, 1:imageSizeX);
circlePixels = (rowsInImage - x(1)).^2 + (columnsInImage - x(2)).^2 <= x(3).^2;

%create circle pixels
% circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;

masked_pixels= mask(circlePixels);

score= sum(masked_pixels)- penalty*(length(masked_pixels)-sum(masked_pixels));
% score=-score;  %Minimize it
end
