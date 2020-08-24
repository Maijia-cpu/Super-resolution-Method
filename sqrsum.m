function [para_coeff,sqminz]=sqrsum(z,radius,index,D_Airy,lambda,angle0,sd_ag)

load sim.mat % angle included variable angle
exp_line=exp_line';
[max_ze,layere]=max(expz);%layere max of exp data
angle0=angle0;
sd_ag=sd_ag;

param.Image_X=327;%Height of the image,9nm/pixel consider PSF xsize 180 and zero padding
param.radius=radius;%input radius 9nm/pixel
param.centerX=round((param.Image_X+1)/2);
param.Image_Z=2*floor((12*250*2/9+2*radius+2000/9)/2)+1;%9 nm/pixel PSF axial,odd number
param.centerZ=round((param.Image_Z+1)/2);
image=create_image(param);% change to 9nm/pixel in lateral and 9nm/pixel in z
image=image';%first dimenstion lateral, second axial

psf=PSFGL(z,index,D_Airy,lambda);% 9nm/pixel in lateral and 9nm/pixel in z
result=DL2.CONV(image,psf,'c');

result_avgz=mean(result);%get z profile
[max_result,resz]=max(result_avgz);%resz max z index
[para_coeff,sqminz]=scan(result,layere,resz,expz,exp_line,angle0,sd_ag);

end