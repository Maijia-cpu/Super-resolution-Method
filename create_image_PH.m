%% create pinhole image
function I=create_image_PH(param)

rowsInImage=linspace(1,param.Image_X,param.Image_X);
columnsInImage=linspace(1,param.Image_Y,param.Image_Y);
pagesInImage=linspace(1,param.Image_Z,param.Image_Z);

[xm,ym] = ndgrid(rowsInImage,columnsInImage) ;

I=zeros(param.Image_X,param.Image_Y,param.Image_Z);
%I=uint16(I);
Radius = param.radius;
centerX=param.centerX;
centerY=param.centerY;
circle= (xm-centerX).^2 + (ym-centerY).^2;
I(:,:,round((param.Image_Z+1)/2))= double(circle <= Radius.^2);

end
