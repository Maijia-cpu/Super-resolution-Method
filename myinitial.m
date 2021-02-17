clear 
clc

%%% convolution plugin downloaded from http://bigwww.epfl.ch/deconvolution/deconvolutionlab2/
javaaddpath /gpfs/loomis/project/howard/ml2542/configuration/DeconvolutionLab_2.jar
javaaddpath /gpfs/loomis/project/howard/ml2542/configuration/JTransforms-3.1-with-dependencies.jar

load branch_info.mat 
% branch_info.mat contains: 
% 1.zvalue (upper bound of distance relative to coverslip,unit: 9nm/pixel, usually zvalues varies from 2000-4000); 
% 2. FWHM (measured FWHM from the dendrite image, unit:nm) 
% 3. angle (some dendrites are tilted relative to horizontal, angle describes the relative angle between dendrites and horizontal) 
% 4. sd_angle (standard deviation of meausred angle)

angle0=angle0;
sd_ag=sd_ag;

%%% range of tuning parameters
  rng('shuffle');

   % range of axial distance relative to coverslip
   a1=1000;%9nm/pixel
   b1=zvalue;% 9nm/pixel
   z=(b1-a1)*rand()+a1; 
   
   % range of refractive index
   a5=1.33;
   b5=1.55;
   index=(b5-a5)*rand()+a5;
   
   % range of pinhole size
   a6=1.6;% SDC pinhole size 1.6AU, but there might be uncertainty
   b6=3;
   D_Airy=(b6-a6)*rand()+a6;
   
     % range of estimated radius
   a2=((FWHM-200)/2-100)/9;
   b2=((FWHM-200)/2+150)/9;
     if a2>2
      radius=(b2-a2)*rand()+a2;%9nm/pixel
     else
      a2=2;
      radius=(b2-a2)*rand()+a2; 
     end
    
    % range of emission wavelength
    a7=500;% emission filter 500-550nm
    b7=550;
    lambda=(b7-a7)*rand()+a7;

    num_exp=8;% number of image stacks usually in the range of 8-12

%%% T: "temperature" in MC simulations. Please determine according to simulation output value h  
    for j=1:1:1000
      if j<200
      T(j)=(0.2-j*0.1/200)*0.25;  
      elseif (j>199)
      T(j)=(0.1-0.06/400*(j-200))*0.25;
      end
    end

 %%% initial values preparations   
            jsq=0;
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

 %%% random numbers used in MC simulations           
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

save('initial.mat')


