function I=create_image(param)

[rowsInImage columnsInImage ] = meshgrid(1:param.Image_X, 1:param.Image_Z);

%input in 9nm/pixel unit
centerX=param.centerX;
centerZ=param.centerZ;
radius=param.radius;%9nm/pixel
circlePixels=(rowsInImage-centerX).^2+(columnsInImage-centerZ).^2<=radius.^2;
circlePixels_sub=(rowsInImage-centerX).^2+(columnsInImage-centerZ).^2<=(radius-1).^2;
I=circlePixels-circlePixels_sub;
% [nz nx]=size(I(:,:));
% I_new=imresize(I, [round(nz/2.78) round(nx)],'bilinear'); % z res 25nm/pixel x,y res 9nm/pixel

end