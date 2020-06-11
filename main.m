%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These codes can only be used for academic research freely. 
% For other purposes, please contact Yanmei Liang (ymliang@nankai.edu.cn).
%
% This file simulates symmetrical illumination based color Fourier
% ptychographic microscopy.
% 
% ref
% Zhang M, Yang D, Liang Y, Color Fourier ptychographic microscopy based on symmetrical illumination and wavelength multiplexing, J Opt, 22, 065604 (2020).
% 
% 
% Several files borrowed Lei Tian's open source.
% ref
% Tian L, Li X, Ramchandran K, Waller L, Multiplexed coded illumination for Fourier Ptychography with an LED array microscope, Biomed Opt Express, 5, 2376-2389 (2014).
%
%
% last modified on 6/05/2020
% by Muyang Zhang, m765421@126.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

%% Light Source Param
waveLength = [0.632e-6; 0.532e-6; 0.47e-6;];
spsize     = 4.8e-6/4;
psize      = spsize/4;
NA         = 0.13;
k0         = 2*pi./waveLength;
%% LED Array Param Square
LED_num_per_cir = [1 8 16 24 32 40 48 56]; % LED for each circle
snum            = 7;
arraysize       = 2*snum+1;
LED_num_sum     = sum(LED_num_per_cir);
LEDgap          = 10;
H               = 190;
zr              = -20e-6;%set defocus distance
zg              = -20e-6;
zb              = -20e-6;
syn_ill         = 'sin';%'sym' or  'sin'
zr_c            = 0;%defocus distance used in this file
zg_c            = 0;
zb_c            = 0;
Img_num         = 225;
LED_num_per_lit = zeros([1,Img_num]);

%% LED Wave Vector Square
xlocation=zeros(3,(2*snum+1)^2);
ylocation=zeros(3,(2*snum+1)^2);
for ch = 1:3
    for i=1:arraysize
        xlocation(ch,1+arraysize*(i-1):15+arraysize*(i-1))...
            = (-(arraysize-1)/2:1:(arraysize-1)/2)*LEDgap;
        ylocation(ch,1+arraysize*(i-1):15+arraysize*(i-1))...
            = ((arraysize-1)/2-(i-1))*LEDgap;
    end
end
kx_relative = -sin(atan(xlocation/H));
ky_relative = -sin(atan(ylocation/H));
ch_r_seq = gseq_dir_assign(arraysize, 'S');%
ch_g_seq = gseq_dir_assign(arraysize, 'S');
ch_b_seq = gseq_dir_assign(arraysize, 'S');

%% Simulation source
objectAmplitude = double(imread('fl3.tiff'));
% objectAmplitude(:,:,end) = [];
objectAmplitude = imresize(objectAmplitude,[600 600]);
objectAmplitude = mat2gray(objectAmplitude)*255;
phase           = double(imread('flower.tiff'));
phase           = imresize(phase,[600 600]);
phase           = mat2gray(phase);
phase           = (0.2*pi*phase - 0.1*pi)*0;%  ignore phase or you can keep it
object          = objectAmplitude.*exp(1i.*phase);
figure;imshow(uint8(abs(object)));title('Input Intensity');
imwrite(mat2gray(abs(object)),'Intensity.tiff');
figure;imshow(angle(object),[]);title('Input Phase');
imwrite(mat2gray(angle(object)),'Phase.tiff');

%% Genarate LR images with multi LED
[m, n, dim]     = size(object);
ml              = m/(spsize/psize);
nl              = n/(spsize/psize);
imSeqLowRes     = zeros(ml, nl, LED_num_sum, 3);
kx              = kx_relative;
ky              = ky_relative;
kx(1,:)         = k0(1)*kx_relative(1,:);
kx(2,:)         = k0(2)*kx_relative(2,:);
kx(3,:)         = k0(3)*kx_relative(3,:);
ky(1,:)         = k0(1)*ky_relative(1,:);
ky(2,:)         = k0(2)*ky_relative(2,:);
ky(3,:)         = k0(3)*ky_relative(3,:);
dkx             = 2*pi/(psize*n);
dky             = 2*pi/(psize*m);
cutoffFrequency = NA*k0;
kmax            = pi/spsize;
[kxm, kym]      = meshgrid(-kmax:kmax/((nl-1)/2):kmax,-kmax:kmax/((ml-1)/2):kmax);
CTF_r           = ((kxm.^2+kym.^2)<cutoffFrequency(1)^2);
CTF_g           = ((kxm.^2+kym.^2)<cutoffFrequency(2)^2);
CTF_b           = ((kxm.^2+kym.^2)<cutoffFrequency(3)^2);
CTF             = cat(3, CTF_r, CTF_g, CTF_b);
kzm_r           = sqrt(k0(1)^2-kxm.^2-kym.^2);
kzm_g           = sqrt(k0(2)^2-kxm.^2-kym.^2);
kzm_b           = sqrt(k0(3)^2-kxm.^2-kym.^2);
pupil_r         = exp(1i.*zr.*real(kzm_r)).*exp(-abs(zr).*abs(imag(kzm_r)));
pupil_g         = exp(1i.*zg.*real(kzm_g)).*exp(-abs(zg).*abs(imag(kzm_g)));
pupil_b         = exp(1i.*zb.*real(kzm_b)).*exp(-abs(zb).*abs(imag(kzm_b)));
pupil           = cat(3, pupil_r, pupil_g, pupil_b);
objectFT_r      = fftshift(fft2(object(:,:,1)));%
objectFT_g      = fftshift(fft2(object(:,:,2)));%
objectFT_b      = fftshift(fft2(object(:,:,3)));%
objectFT        = cat(3, objectFT_r, objectFT_g, objectFT_b);
HR_CTF_base     = zeros(m, n, LED_num_sum, 3);
for ch = 1:3
    for tt = 1:LED_num_sum
        kxc = round((n+1)/2+kx(ch,tt)/dkx);
        kyc = round((m+1)/2+ky(ch,tt)/dky);
        kyl = round(kyc-(ml-1)/2);kyh = round(kyc+(ml-1)/2);
        kxl = round(kxc-(nl-1)/2);kxh = round(kxc+(nl-1)/2);
        imSeqLowFT                         = (ml/m)^2*objectFT(kyl:kyh,kxl:kxh,ch).*CTF(:,:,ch).*pupil(:,:,ch);
        imSeqLowRes(:,:,tt,ch)             = abs(ifft2(ifftshift(imSeqLowFT)));
        HR_CTF_base(kyl:kyh,kxl:kxh,tt,ch) = CTF(:,:,ch);
    end
end
LR_images = imSeqLowRes.^2;
pic_l = zeros([ml,nl,3]);
pic_l(:,:,1) = imSeqLowRes(:,:,113,1);
pic_l(:,:,2) = imSeqLowRes(:,:,113,2);
pic_l(:,:,3) = imSeqLowRes(:,:,113,3);
imwrite(uint8(pic_l),'low_center.tiff','Resolution',600);

%% LED light pattern and syn images
base = cell(1,length(LED_num_per_cir));
base{1} = 1;
for i = 2:length(LED_num_per_cir)
    base{i} = sum(LED_num_per_cir(1:(i-1))) + (1:LED_num_per_cir(i));
end
syn_img = zeros([ml, nl, Img_num, 3]);
Ns      = cell(1,Img_num);
scale   = cell(1,Img_num);
light_pattern = cell(1,Img_num);
flag    = zeros([1,Img_num]);
num = 1;

kx(1,:)         = k0(1)*kx_relative(1,:);
kx(2,:)         = k0(2)*kx_relative(2,:);
kx(3,:)         = k0(3)*kx_relative(3,:);
ky(1,:)         = k0(1)*ky_relative(1,:);
ky(2,:)         = k0(2)*ky_relative(2,:);
ky(3,:)         = k0(3)*ky_relative(3,:);
for i = 1:length(LED_num_per_cir)
    for j = 1:LED_num_per_cir(i)
        light_pattern{num}   = [ch_r_seq(num) ch_g_seq(num) ch_b_seq(num)];
        
        if i<=1 %central, do not syn 
            LED_num_per_lit(num) = 1;
            syn_img(:,:,num,:)   = cat(3,LR_images(:,:,ch_r_seq(num),1),LR_images(:,:,ch_g_seq(num),2),LR_images(:,:,ch_b_seq(num),3));
            flag(num)            = 1;
        else
            LED_num_per_lit(num) = 3;
            syn_img(:,:,num,1)   = LR_images(:,:,ch_r_seq(num),1)+LR_images(:,:,ch_g_seq(num),2)+LR_images(:,:,ch_b_seq(num),3);
            flag(num)            = 0;
        end
        
        Ns{num}          = [round((m+1)/2+ky(1,light_pattern{num}(1))/dky)' round((n+1)/2+kx(1,light_pattern{num}(1))/dkx)';
                            round((m+1)/2+ky(2,light_pattern{num}(2))/dky)' round((n+1)/2+kx(2,light_pattern{num}(2))/dkx)';
                            round((m+1)/2+ky(3,light_pattern{num}(3))/dky)' round((n+1)/2+kx(3,light_pattern{num}(3))/dkx)';];
        scale{num}       = ones([1,LED_num_per_lit(num)])/LED_num_per_lit(num);
        num              = num + 1;       
    end
end
if syn_ill == 'sym'
    for i = 2:length(LED_num_per_cir)
        for j = 1:LED_num_per_cir(i)
            if j>(LED_num_per_cir(i)/2)
                continue;
            end
            temp = syn_img(:,:,base{i}(j),1) + syn_img(:,:,base{i}(j+LED_num_per_cir(i)/2),1);
            syn_img(:,:,base{i}(j),1) = temp/2;
            syn_img(:,:,base{i}(j+LED_num_per_cir(i)/2),1) = temp/2;
        end
    end
end
%% Reconstruct
pupil_r         = exp(1i.*zr_c.*real(kzm_r)).*exp(-abs(zr_c).*abs(imag(kzm_r)));
pupil_g         = exp(1i.*zg_c.*real(kzm_g)).*exp(-abs(zg_c).*abs(imag(kzm_g)));
pupil_b         = exp(1i.*zb_c.*real(kzm_b)).*exp(-abs(zb_c).*abs(imag(kzm_b)));
pupil           = cat(3, pupil_r, pupil_g, pupil_b);
F = @(x) fftshift(fft2(x));
Ft = @(x) ifft2(ifftshift(x));
row = @(x) x(:).';
cen0 = [round((m+1)/2) round((n+1)/2)];

downsamp = @(x,cen) x(round(cen(1)-(ml-1)/2):round(cen(1)+(ml-1)/2),...
    round(cen(2)-(nl-1)/2):round(cen(2)+(nl-1)/2));

tr   = double(CTF(:,:,1)).*F(sqrt(syn_img(:,:,1,1)));
Or   = padarray(tr,[(m-ml)/2 (n-nl)/2]);
tg   = double(CTF(:,:,2)).*F(sqrt(syn_img(:,:,1,2)));
Og   = padarray(tg,[(m-ml)/2 (n-nl)/2]);
tb   = double(CTF(:,:,3)).*F(sqrt(syn_img(:,:,1,3)));
Ob   = padarray(tb,[(m-ml)/2 (n-nl)/2]);
O    = cat(3, Or, Og, Ob);
P    = double(CTF);
Ps   = double(CTF);
loop = 3;
H0   = pupil;
OP_alpha = 1;
OP_beta = 1e3;

for iter = 1:loop
    for j = 1:Img_num
        Psi0 = zeros(ml,nl,LED_num_per_lit(j));
        Psi_scale = zeros(ml,nl,LED_num_per_lit(j));
        cen = zeros(2,LED_num_per_lit(j));
        scale0 = zeros(LED_num_per_lit(j),1);
        if flag(j)
            for p = 1:3
                kyc = Ns{j}(p,1);
                kxc = Ns{j}(p,2);
                kyl = round(kyc-(ml-1)/2);kyh = round(kyc+(ml-1)/2);
                kxl = round(kxc-(nl-1)/2);kxh = round(kxc+(nl-1)/2);
                cen(1) = kyc;
                cen(2) = kxc;
                scale0 = scale{j}(1);
                Psi0(:,:,1) = O(kyl:kyh,kxl:kxh,p).*P(:,:,p).*H0(:,:,p);
                Psi_scale(:,:,1) = sqrt(scale0)*Psi0(:,:,1);
                I_mea = syn_img(:,:,j,p);
                psi0 = Ft(Psi_scale);
                I_est = sum(abs(psi0).^2,3);
                Psi = Proj_Fourier_v2(psi0, I_mea, I_est, scale0, F);
                Psi2 = Psi;
                dPsi = Psi2-Psi0;
                Omax = abs(O(cen0(1),cen0(2),1));
                P2 = @(O,P,dpsi,Omax,cen,Ps)...
                GDUpdate_Multiplication_rank1(O,P,dpsi,Omax,cen,Ps,...
                OP_alpha,OP_beta);
                [O(:,:,p),P(:,:,p)] = P2(O(:,:,p),P(:,:,p),dPsi./repmat(H0(:,:,p),[1,1,LED_num_per_lit(j)]),Omax,cen+1,Ps(:,:,p));
            end
        else
            for p = 1:LED_num_per_lit(j)
                kyc = Ns{j}(p,1);
                kxc = Ns{j}(p,2);
                kyl = round(kyc-(ml-1)/2);kyh = round(kyc+(ml-1)/2);
                kxl = round(kxc-(nl-1)/2);kxh = round(kxc+(nl-1)/2);
                cen(1,p) = kyc;
                cen(2,p) = kxc;
                scale0(p) = scale{j}(p);
                Psi0(:,:,p) = O(kyl:kyh,kxl:kxh,p).*P(:,:,p).*H0(:,:,p);
                Psi_scale(:,:,p) = sqrt(scale0(p))*Psi0(:,:,p);
            end
            I_mea = syn_img(:,:,j,1);
            psi0 = Ft(Psi_scale);
            I_est = sum(abs(psi0).^2,3);
            Psi = Proj_Fourier_v2(psi0, I_mea, I_est, scale0, F);
            Psi2 = cat(3, Psi(:,:,3), Psi(:,:,1:2));
            dPsi = Psi2-Psi0;
            Omax = abs(O(cen0(1),cen0(2)));
            P2 = @(O,P,dpsi,Omax,cen,Ps)...
                GDUpdate_Multiplication_rank_rm(O,P,dpsi,Omax,cen,Ps,...
                OP_alpha,OP_beta);
            [O,P] = P2(O,P,dPsi./H0,Omax,cen+1,Ps);
        end
    end  
    iter
    iter = iter+1;
end
Ir = abs(ifft2(ifftshift(O(:,:,1))));
Ig = abs(ifft2(ifftshift(O(:,:,2))));
Ib = abs(ifft2(ifftshift(O(:,:,3))));
figure,imshow(cat(3,Ir/mean2(Ir)*mean2(pic_l(:,:,1)),Ig/mean2(Ig)*mean2(pic_l(:,:,2)),Ib/mean2(Ib)*mean2(pic_l(:,:,3)))/255),title('Reconstruction')
imwrite(uint8(cat(3,Ir/mean2(Ir)*mean2(pic_l(:,:,1)),Ig/mean2(Ig)*mean2(pic_l(:,:,2)),Ib/mean2(Ib)*mean2(pic_l(:,:,3)))),'Reconstruction.tiff','Resolution',600);
