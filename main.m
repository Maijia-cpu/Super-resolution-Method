clear 
clc

javaaddpath /gpfs/loomis/scratch60/fas/howard/ml2542/file/DeconvolutionLab_2.jar
javaaddpath /gpfs/loomis/scratch60/fas/howard/ml2542/file/JTransforms-3.1-with-dependencies.jar 
% javaaddpath /home/ml2542/project/MATLAB/DeconvolutionLab_2.jar
% javaaddpath /home/ml2542/project/MATLAB/JTransforms-3.1-with-dependencies.jar % %  lambda=530;

load initial.mat
num_exp=8;


while jsq<600

num=randnum(jsq+1);% random generate 0 or 1 or 2; 0:zvalue 1:radius 2:index,3 Airy
disp(num);
if num==0
    z1=(b1-a1)*A(jsq+500)+a1;%9nm/pixel
    [para_coeff1,h1]=sqrsum(z1,radius,index,D_Airy,lambda,angle0,sd_ag);
    jsq=jsq+1;
    Z=exp(-h1/T(jsq))/exp(-h/T(jsq));
    disp(jsq);
    if Z<1
        num1=B(500+jsq);
        if Z>num1
            z=z1;
            h=h1;
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z1;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda;
            param(jsq+1,6)=h1;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z1;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda;
            record(jsq+1,6)=h1;
            
            across(jsq+1,:)=para_coeff1;%coeff
        else
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda;
            param(jsq+1,6)=h;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z1;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda;
            record(jsq+1,6)=h1;
 
            across(jsq+1,:)=across(jsq,:);%layer int difference
        end
    else
            z=z1;
            h=h1;
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z1;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda;
            param(jsq+1,6)=h1;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z1;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda;
            record(jsq+1,6)=h1;
            
            across(jsq+1,:)=para_coeff1;%coeff
    end
        
elseif num==1 %radius 9nm/pixel
    radius1=(b2-a2)*C(jsq+500)+a2;
    [para_coeff1,h1]=sqrsum(z,radius1,index,D_Airy,lambda,angle0,sd_ag);
    jsq=jsq+1;
    Z=exp(-h1/T(jsq))/exp(-h/T(jsq));
    disp(jsq);
    if Z<1
        num1=D(jsq+500);
        if Z>num1
            radius=radius1;
            h=h1;
            param(jsq+1,1)=radius1;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda;
            param(jsq+1,6)=h1;
            
            record(jsq+1,1)=radius1;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda;
            record(jsq+1,6)=h1;
            
            across(jsq+1,:)=para_coeff1;%coeff
        else
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda;
            param(jsq+1,6)=h;
            
            record(jsq+1,1)=radius1;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda;
            record(jsq+1,6)=h1;

            across(jsq+1,:)=across(jsq,:);%layer int difference
        end
    else
            radius=radius1;
            h=h1;
            param(jsq+1,1)=radius1;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda;          
            param(jsq+1,6)=h;
            
            record(jsq+1,1)=radius1;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda;
            record(jsq+1,6)=h1;
            
            across(jsq+1,:)=para_coeff1;%coeff
    end
elseif num==2
     index1=(b5-a5)*U(jsq+500)+a5;
     [para_coeff1,h1]=sqrsum(z,radius,index1,D_Airy,lambda,angle0,sd_ag);
     jsq=jsq+1;
     Z=exp(-h1/T(jsq))/exp(-h/T(jsq));
     disp(jsq);
    if Z<1
        num1=V(jsq+500);
        if Z>num1
            index=index1;
            h=h1;
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index1;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda;  
            param(jsq+1,6)=h1;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index1;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda;  
            record(jsq+1,6)=h1;
            
            across(jsq+1,:)=para_coeff1;%coeff
        else
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda;  
            param(jsq+1,6)=h;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index1;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda;  
            record(jsq+1,6)=h1;

            across(jsq+1,:)=across(jsq,:);%layer int difference

        end
    else
            index=index1;
            h=h1;
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index1;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda;  
            param(jsq+1,6)=h1;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index1;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda;
            record(jsq+1,6)=h1;
            
            across(jsq+1,:)=para_coeff1;%coeff

    end
elseif num==3
    D_Airy1=(b6-a6)*M(jsq+500)+a6;
    [para_coeff1,h1]=sqrsum(z,radius,index,D_Airy1,lambda,angle0,sd_ag);
    jsq=jsq+1;
    Z=exp(-h1/T(jsq))/exp(-h/T(jsq));
    disp(jsq);
    if Z<1
        num1=K(jsq+500);
        if Z>num1
            D_Airy=D_Airy1;
            h=h1;
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy1;
            param(jsq+1,5)=lambda; 
            param(jsq+1,6)=h1;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy1;
            record(jsq+1,5)=lambda; 
            record(jsq+1,6)=h1;
            
            across(jsq+1,:)=para_coeff1;%coeff
            
        else
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda; 
            param(jsq+1,6)=h;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy1;
            record(jsq+1,5)=lambda; 
            record(jsq+1,6)=h1;

            across(jsq+1,:)=across(jsq,:);%layer int difference

        end
    else
            D_Airy=D_Airy1;
            h=h1;
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy1;
            param(jsq+1,5)=lambda; 
            param(jsq+1,6)=h1;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy1;
            record(jsq+1,5)=lambda; 
            record(jsq+1,6)=h1;
            
            across(jsq+1,:)=para_coeff1;%coeff

    end
else 
    lambda1=(b7-a7)*LA(jsq+500)+a7;
    [para_coeff1,h1]=sqrsum(z,radius,index,D_Airy,lambda1,angle0,sd_ag);
    jsq=jsq+1;
    Z=exp(-h1/T(jsq))/exp(-h/T(jsq));
    disp(jsq);
    if Z<1
        num1=LB(jsq+500);
        if Z>num1
            lambda=lambda1;
            h=h1;
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda1;
            param(jsq+1,6)=h1;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda1;
            record(jsq+1,5)=h1;
            
            across(jsq+1,:)=para_coeff1;%coeff
            
        else
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
            record(jsq+1,5)=lambda1;
            record(jsq+1,6)=h1;

            across(jsq+1,:)=across(jsq,:);%layer int difference

        end
    else
            lambda=lambda1;
            h=h1;
            param(jsq+1,1)=radius;
            param(jsq+1,2)=z;
            param(jsq+1,3)=index;
            param(jsq+1,4)=D_Airy;
            param(jsq+1,5)=lambda1;
            param(jsq+1,6)=h1;
            
            record(jsq+1,1)=radius;
            record(jsq+1,2)=z;
            record(jsq+1,3)=index;
            record(jsq+1,4)=D_Airy;
            record(jsq+1,5)=lambda1;
            record(jsq+1,6)=h1;
            
            across(jsq+1,:)=para_coeff1;%coeff

    end   
end
save('initial.mat')

end