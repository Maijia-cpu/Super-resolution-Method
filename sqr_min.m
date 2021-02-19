%% obtained the least squared difference between experimental image and simulated image
function [sq_min]=sqr_min(imagecg_norm,bg_im,expz,exp_line)
  
    dnum=length(bg_im);
    image_norm=repmat(imagecg_norm,[1,1,dnum]);
    bg_im_mat1=repmat(bg_im,[1,25]);
    bg_im_mat2=repmat(bg_im_mat1,[1,1,8]);
    bg_im_mat=permute(bg_im_mat2,[2 3 1]);  
    imagenorm=image_norm+bg_im_mat;

[maxline,line_layer]=max(imagenorm);
[aj]=find(~maxline);
[ajn]=find(maxline);%nonzero
num=length(expz);
if length(ajn)==num*dnum
 imagenorm_nm=imagenorm./maxline;
else
 imagenorm_nm(:,aj)=imagenorm(:,aj);
 imagenorm_nm(:,ajn)=imagenorm(:,ajn)./(maxline(:,ajn));
end

expline_norm=repmat(exp_line,[1,1,dnum]);

%% calculate difference
diff_image=(imagenorm_nm-expline_norm).^2;
diff_image_s1=sum(diff_image);

weight_expz=expz./sum(expz);
weight_norm=repmat(weight_expz,[1,1,dnum]);

h_df=diff_image_s1.*weight_norm;

sq_min=5*squeeze(sum(h_df));
end
