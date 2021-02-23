function [para_coeff_m,sqminz_m]=scan(result,layere,resz,expz,exp_line,angle0,sd_ag)
    %% give delta x and delta z reasonable range
    delta=160/9;%9nm/pixel; 160 here is used step size in experiments
    num_exp=length(expz);
    delta_s=delta*(zeros(40,1)+1);
    rng('shuffle');
   
    x_frac=rand(800,1500);
    z_frac=rand(500,1500);
    b_frac=rand(500,2600);
  
    j=0;
    jsq=0;
dif_im=zeros(30,1);
para_coeff=zeros(4,30);
%% loop
%  tic
while j<30
    j=j+1;
startz=resz-ceil((layere+0.7)*delta+1)+round(3.5*z_frac(300:339,j+1000)*delta);% tuning range for the starting layer position, 0.7 and 3.5 are tunable based on experimental condition
zlayer=round(linspaceNDim(startz,startz+(num_exp-1)*delta_s,num_exp));

for i=1:40    
    imagez=result(:,zlayer(i,:));% 9nm/pixel resolution  
    im_c=zeros(327,1);
    md_im=(1+length(im_c))/2;
    im_c(md_im)=1;
    IAx=(27*x_frac(500+i,j+1000)+1):1:(27*x_frac(500+i,j+1000)+300);%start of x position 
    image_z=imagez(int16(IAx),:);%layer taken by experiment
    imc=im_c(int16(IAx));
    %% consider rotation
    N=12*20;
    imagezp=repmat(image_z,1,1,N);
    imagez_p=permute(imagezp,[1 3 2]);
    na0=sd_ag/5;
    na=round(na0);  
   for ag=0:1:na
       jsq=jsq+1;
    angle_rd=angle0+(ag-na/2)*5;  
    im_r= imrotate(imagez_p,angle_rd);% counterclockwise
    mask=ones(12,12,1);%9nm/pixel to 108nm/pixel
    imager_cg=convn(im_r,mask,'valid');%coarse-grain
    imr_cg=imager_cg(1:12:end,1:12:end,:);

%% find center point of rotated image
c_mat=repmat(imc,1,N);
c_mat_r=imrotate(c_mat,angle_rd);
c_mat_r_cg=conv2(c_mat_r,mask,'valid');
cr_cg=c_mat_r_cg(1:12:end,1:12:end);
[ycm,xcm]=find(cr_cg);
xm=mean(xcm);
ym=mean(ycm);
   %% rotate perpendicular to line
   th=-(angle_rd/180*pi);
   R = [cos(th) -sin(th) ;sin(th) cos(th)];
   v1=[0;1];
   v2 = R*v1;%clockwise

   for il=1:num_exp
     for  t=-3:1:3
         y_m=ym+t*sin(th);
         x_m=xm+t*cos(th);
         xk1=[x_m-v2(1)*12 x_m+v2(1)*12];
         yk1=[y_m-v2(2)*12 y_m+v2(2)*12];       
%          hold on
%          plot(x_m,y_m,'ro')
%          line(xk1,yk1)
         profile_all((t+4),:)=improfile(imr_cg(:,:,il),xk1,yk1,25); 
     end 
         imagecg(:,il)=sum(profile_all,1);
   end 
    imagecg_norm=double(imagecg/max(imagecg(:)));
    %% add background
    bg_im=0.15*b_frac(100+i,1001+j:2000+j)';
%     tic
    [diff_min]=sqr_min(imagecg_norm,bg_im,expz,exp_line);
%     toc
    [sq_min I_min]=min(diff_min);
    record_data(jsq,1)=sq_min;
    record_data(jsq,2)=I_min;
    record_data(jsq,3)=i;
    record_data(jsq,4)=angle_rd;
    end
end 

[minsq jsq_min]=min(record_data(:,1));% jsq_min-row
jsq_min=jsq_min(1);
dif_im(j)=minsq;
I_num=record_data(jsq_min,2);
imin=record_data(jsq_min,3);
para_coeff(1,j)=x_frac(500+imin,j+1000);
para_coeff(2,j)=z_frac(300+imin-1,j+1000);
para_coeff(3,j)=b_frac(100+imin,j+1001+I_num-1);
para_coeff(4,j)=record_data(jsq_min,4);
clear record_data jsq jsq_min
jsq=0;
end
%  toc
[sqminz_m,I_m]=min(dif_im);
para_coeff_m=para_coeff(:,I_m);

end
