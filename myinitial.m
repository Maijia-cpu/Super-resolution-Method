clear 
clc

javaaddpath /gpfs/loomis/scratch60/fas/howard/ml2542/file/DeconvolutionLab_2.jar
javaaddpath /gpfs/loomis/scratch60/fas/howard/ml2542/file/JTransforms-3.1-with-dependencies.jar
% javaaddpath /Users/maijialiao/Documents/biophysics/data/result/test_on_exp_image/Formal_data/matlab_configuration/DL/DeconvolutionLab_2.jar
% javaaddpath /Users/maijialiao/Documents/biophysics/data/result/test_on_exp_image/Formal_data/matlab_configuration/DL/JTransforms-3.1-with-dependencies.jar 

load branch_info.mat % contain zvalue FWHM angle sd_angle
angle0=angle0;
sd_ag=sd_ag;
 %% main part
  rng('shuffle');

   a1=zvalue;%9nm/pixel
   b1=zvalue+600;% 9nm/pixel
   z=(b1-a1)*rand()+a1;  
   
   a5=1.33;
   b5=1.45;
   index=(b5-a5)*rand()+a5;
   
   a6=1.2;% I used 1.6AU, but there might be uncertainty
   b6=2.0;
   D_Airy=(b6-a6)*rand()+a6;
   
  a2=((FWHM-200)/2-100)/9;
  b2=((FWHM-200)/2+150)/9;
  if a2>4
  radius=(b2-a2)*rand()+a2;%9nm/pixel
  else
  a2=4;
  radius=(b2-a2)*rand()+a2; 
  end
  
  a7=500;% emission filter 500-550nm
  b7=550;
  lambda=(b7-a7)*rand()+a7;
   
  jsq=0;

  num_exp=8;
%% change to linear varying T  
for j=1:1:600
if j<200
T(j)=(0.2-j*0.1/200)*0.25;  %%%tunable
elseif (j>199)
T(j)=(0.1-0.06/400*(j-200))*0.25;
end
end

[para_coeff,h]=sqrsum(z,radius,index,D_Airy,lambda,angle0,sd_ag);

            param(jsq+1,1)=radius;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda;
            param(jsq+1,6)=h;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda;
            record(jsq+1,6)=h;           
            across(jsq+1,:)=para_coeff;%coeff

  jsqn=0;
  
while jsqn<1500
number=floor(5*rand);
jsqn=jsqn+1;
randnum(jsqn)=number;
end

   V=rand(1500,1);
   U=rand(1500,1);
   A=rand(1500,1);
   B=rand(1500,1);
   C=rand(1500,1);
   D=rand(1500,1);
   M=rand(1500,1);
   K=rand(1500,1);
   LA=rand(1500,1);
   LB=rand(1500,1);

save('initial.mat')%change from -v7.3 to -v6, assume -v6 is faster since -v7.3 takes time in compression


