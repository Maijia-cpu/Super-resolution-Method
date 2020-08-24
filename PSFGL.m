%% Demo for the MicroPSF, a fast and accurate approximation of the Gibson-Lanni model
% If the parameters are not assigned, default setting will be loaded, see
% MicroscPSF for details.
%
%   Reference:
%       [1] Gibson, S.F. and Lanni, F., 1992.
%           Experimental test of an analytical model of aberration in an
%           oil-immersion objective lens used in three-dimensional light
%           microscopy. JOSA A, 9(1), pp.154-166.
%       [2] Li, J., Xue, F. and Blu, T. Fast and accurate 3D PSF
%           computation for fluorescence microscopy. JOSA A. Accepted.
%
%   Copyright � Jizhou Li, Feng Xue and Thierry Blu, 2017
%   Update date: 4 May, 2017

function psftd=PSFGL(z,index,D_Airy,lambda)
addpath('Utilities/');

params.NA=1.2;
params.lambda = 488e-9;
params.M = 60;
params.resLateral = 9e-9;
params.resAxial = 9e-9;
params.pZ = z*9*10^(-9);%9nm res
size_z=(round(params.pZ/params.resAxial))*2+400*2+1;%  more space in addition to particle
params.size = [251 251 size_z];
params.ns = index;
params.ni0=1.33;
params.ni = 1.33;
params.tg0 = 170e-6;
params.tg = 170e-6;

tic;
PSFi = MicroscPSF(params); %max psf independent of depth z when index=1.33 decrease when 
t = toc;
[psfmax lmax]=max(squeeze(PSFi(126,126,:)));

z_p=min([400,lmax])-1;
PSF_c=PSFi(:,:,(lmax-z_p):(lmax+z_p));

paramsem.NA=1.2;
paramsem.lambda = 1e-9*lambda;
paramsem.M = 60;
paramsem.resLateral = 9e-9;
paramsem.resAxial = 9e-9;
paramsem.pZ = z*9*10^(-9);%10nm res
size_zem=(round(paramsem.pZ/paramsem.resAxial))*2+400*2+1;
paramsem.size = [251 251 size_zem];
paramsem.ns = index;
paramsem.ni0=1.33;
paramsem.ni = 1.33;
paramsem.tg0 = 170e-6;
paramsem.tg = 170e-6;

tic;
PSFem = MicroscPSF(paramsem);
t = toc;

PSF_em=PSFem(:,:,(lmax-z_p):(lmax+z_p));

D_PH=1.22*params.lambda/params.NA*D_Airy/params.resLateral;%diamter
sz=size(PSF_em);
parami.Image_X=sz(1);
parami.Image_Y=sz(2);
parami.Image_Z=sz(3);
parami.radius=D_PH/2;
parami.centerX=(parami.Image_X+1)/2;
parami.centerY=(parami.Image_Y+1)/2;
PH=create_image_PH(parami);
PSF_phg=DL2.CONV(PSF_em,PH,'c');% new version of deconvolution software, work!

psf_3d=PSF_phg.*PSF_c;
psf_2d=squeeze(mean(psf_3d,2));
psf2d=flip(psf_2d,2);%reverse direction, left close to coverslip,%determined from tetraspeck beads
psf_avgz=mean(psf_2d);%get z profile
[max_result,pz]=max(psf_avgz);%resz max z index
psf_pz=min([pz,(length(psf_avgz)-pz)])-1;%center at maximum
psf_td=psf_2d(:,(pz-psf_pz):(pz+psf_pz));
psftd=psf_td/max(psf_td(:));
end

