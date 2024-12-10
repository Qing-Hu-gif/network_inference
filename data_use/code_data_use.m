
%% 1. calculate WT R0 r0 J0 based on ode
% using ode method to calculate R0 r0 J0 at unperturbed E, H, M

clear; clc; close all
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\WT'

s0=load('EMT9_all1_EHM_k03_0012_kd3_03.txt');
para=[0.6 0.09 1.66 0.03 0.09 1.66 0.035 0.5 0.5];   
load J_sym1_0.mat
syms T s S R3 z Z R2 E N
% calculate R0 r0 J0
ode_JEHM=zeros(length(para),length(para),3); % 1 E, 2 H, 3 M
ode_REHM=zeros(length(para),length(para),3); % 1 E, 2 H, 3 M
ode_REHM_ln=zeros(length(para),length(para),3); % 1 E, 2 H, 3 M
ode_REHM_dln=zeros(length(para),length(para),3); % 1 E, 2 H, 3 M
ode_rEHM=zeros(length(para),length(para),3); % 1 E, 2 H, 3 M
ode_rEHM_dln=zeros(length(para),length(para),3); % 1 E, 2 H, 3 M
for i=1:3
    ode_JEHM(:,:,i)=double(subs(J_sym_0,[T,s,S,R3,z,Z,R2,E,N],s0(i,2:10)));
    ode_REHM_dln(:,:,i)=inv(diag(s0(i,2:10)))*inv(ode_JEHM(:,:,i))*diag(s0(i,2:10))*diag(para);
    ode_rEHM_dln(:,:,i)=-inv(diag(s0(i,2:10)))*inv(diag(diag(ode_JEHM(:,:,i))))*ode_JEHM(:,:,i)*diag(s0(i,2:10));
end
ode_J_one=zeros(49,3);
ode_R_one_dln=zeros(49,3);
ode_r_one_dln=zeros(49,3);
for i=1:3
    ode_J_one(:,i)=reshape(ode_JEHM(1:7,1:7,i)',49,1);
    ode_R_one_dln(:,i)=reshape(ode_REHM_dln(1:7,1:7,i)',49,1);
    ode_r_one_dln(:,i)=reshape(ode_rEHM_dln(1:7,1:7,i)',49,1);
end
x0=s0(:,2:8)';
clear  T s S R3 z Z R2 E N J_sym_0 i 
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_use' 
% save ode_J_one.mat ode_J_one;
% save ode_R_one_dln.mat ode_R_one_dln;
save ode_rr_one_dln.mat ode_r_one_dln;


%% 2. calculate based on perturbation data
% using data at different perturbation combiantion to calculate R r J
clear;clc
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\WT'
s0=load('EMT9_all1_EHM_k03_0012_kd3_03.txt');
p0=[0.6 0.09 1.66 0.03 0.09 1.66 0.035 0.5 0.5];

% 20lr_u
x_all=cell(3,9);
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\20lr_u\E'
x_all{1,1}=load('EMT9_kdT_E_10000_20lr_u.txt');x_all{1,2}=load('EMT9_kdss_E_10000_20lr_u.txt');
x_all{1,3}=load('EMT9_kdS_E_10000_20lr_u.txt');x_all{1,4}=load('EMT9_kd3_E_10000_20lr_u.txt');
x_all{1,5}=load('EMT9_kdzz_E_10000_20lr_u.txt');x_all{1,6}=load('EMT9_kdZ_E_10000_20lr_u.txt');
x_all{1,7}=load('EMT9_kd2_E_10000_20lr_u.txt');x_all{1,8}=load('EMT9_kde_E_10000_20lr_u.txt');
x_all{1,9}=load('EMT9_kdn_E_10000_20lr_u.txt');
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\20lr_u\H'
x_all{2,1}=load('EMT9_kdT_H_10000_20lr_u.txt');x_all{2,2}=load('EMT9_kdss_H_10000_20lr_u.txt');
x_all{2,3}=load('EMT9_kdS_H_10000_20lr_u.txt');x_all{2,4}=load('EMT9_kd3_H_10000_20lr_u.txt');
x_all{2,5}=load('EMT9_kdzz_H_10000_20lr_u.txt');x_all{2,6}=load('EMT9_kdZ_H_10000_20lr_u.txt');
x_all{2,7}=load('EMT9_kd2_H_10000_20lr_u.txt');x_all{2,8}=load('EMT9_kde_H_10000_20lr_u.txt');
x_all{2,9}=load('EMT9_kdn_H_10000_20lr_u.txt');
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\20lr_u\M'
x_all{3,1}=load('EMT9_kdT_M_10000_20lr_u.txt');x_all{3,2}=load('EMT9_kdss_M_10000_20lr_u.txt');
x_all{3,3}=load('EMT9_kdS_M_10000_20lr_u.txt');x_all{3,4}=load('EMT9_kd3_M_10000_20lr_u.txt');
x_all{3,5}=load('EMT9_kdzz_M_10000_20lr_u.txt');x_all{3,6}=load('EMT9_kdZ_M_10000_20lr_u.txt');
x_all{3,7}=load('EMT9_kd2_M_10000_20lr_u.txt');x_all{3,8}=load('EMT9_kde_M_10000_20lr_u.txt');
x_all{3,9}=load('EMT9_kdn_M_10000_20lr_u.txt');
x=cell(1,3);R_dln=cell(1,3);r_dln=cell(1,3);J_dln=cell(1,3);
r_one_dln=cell(1,3); R_one_dln=cell(1,3); J_one_dln=cell(1,3);
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_use'
i=1;
for j=1:3
    [x{i,j},R_dln{i,j},r_dln{i,j},J_dln{i,j}]=fun_xRrJ_dln(x_all(j,:),s0(j,2:end-1),p0);
end

n=size(r_dln{1,1},1)-2;
for j=1:3
    idr=zeros(n*n,size(r_dln{i,j},3));idR=idr;idJ=idr;
    for k=1:size(idr,2)
        idr(:,k)=reshape(r_dln{i,j}(1:n,1:n,k)',n*n,1);
        idR(:,k)=reshape(R_dln{i,j}(1:n,1:n,k)',n*n,1);
        idJ(:,k)=reshape(J_dln{i,j}(1:n,1:n,k)',n*n,1);
    end
    r_one_dln{i,j}=idr; R_one_dln{i,j}=idR; J_one_dln{i,j}=idJ;
end

% 50lr_u
x_all=cell(3,9);
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\50lr_u\E'
x_all{1,1}=load('EMT9_kdT_E_10000_50lr_u.txt');x_all{1,2}=load('EMT9_kdss_E_10000_50lr_u.txt');
x_all{1,3}=load('EMT9_kdS_E_10000_50lr_u.txt');x_all{1,4}=load('EMT9_kd3_E_10000_50lr_u.txt');
x_all{1,5}=load('EMT9_kdzz_E_10000_50lr_u.txt');x_all{1,6}=load('EMT9_kdZ_E_10000_50lr_u.txt');
x_all{1,7}=load('EMT9_kd2_E_10000_50lr_u.txt');x_all{1,8}=load('EMT9_kde_E_10000_50lr_u.txt');
x_all{1,9}=load('EMT9_kdn_E_10000_50lr_u.txt');
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\50lr_u\H'
x_all{2,1}=load('EMT9_kdT_H_10000_50lr_u.txt');x_all{2,2}=load('EMT9_kdss_H_10000_50lr_u.txt');
x_all{2,3}=load('EMT9_kdS_H_10000_50lr_u.txt');x_all{2,4}=load('EMT9_kd3_H_10000_50lr_u.txt');
x_all{2,5}=load('EMT9_kdzz_H_10000_50lr_u.txt');x_all{2,6}=load('EMT9_kdZ_H_10000_50lr_u.txt');
x_all{2,7}=load('EMT9_kd2_H_10000_50lr_u.txt');x_all{2,8}=load('EMT9_kde_H_10000_50lr_u.txt');
x_all{2,9}=load('EMT9_kdn_H_10000_50lr_u.txt');
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\50lr_u\M'
x_all{3,1}=load('EMT9_kdT_M_10000_50lr_u.txt');x_all{3,2}=load('EMT9_kdss_M_10000_50lr_u.txt');
x_all{3,3}=load('EMT9_kdS_M_10000_50lr_u.txt');x_all{3,4}=load('EMT9_kd3_M_10000_50lr_u.txt');
x_all{3,5}=load('EMT9_kdzz_M_10000_50lr_u.txt');x_all{3,6}=load('EMT9_kdZ_M_10000_50lr_u.txt');
x_all{3,7}=load('EMT9_kd2_M_10000_50lr_u.txt');x_all{3,8}=load('EMT9_kde_M_10000_50lr_u.txt');
x_all{3,9}=load('EMT9_kdn_M_10000_50lr_u.txt');
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_use'
i=2;
for j=1:3
    [x{i,j},R_dln{i,j},r_dln{i,j},J_dln{i,j}]=fun_xRrJ_dln(x_all(j,:),s0(j,2:end-1),p0);
end
for j=1:3
    idr=zeros(n*n,size(r_dln{i,j},3));idR=idr;idJ=idr;
    for k=1:size(idr,2)
        idr(:,k)=reshape(r_dln{i,j}(1:n,1:n,k)',n*n,1);
        idR(:,k)=reshape(R_dln{i,j}(1:n,1:n,k)',n*n,1);
        idJ(:,k)=reshape(J_dln{i,j}(1:n,1:n,k)',n*n,1);
    end
    r_one_dln{i,j}=idr; R_one_dln{i,j}=idR; J_one_dln{i,j}=idJ;
end

% 10n
x_all=cell(3,9);
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\10n\E'
x_all{1,1}=load('EMT9_kdT_E_10000_10n.txt');x_all{1,2}=load('EMT9_kdss_E_10000_10n.txt');
x_all{1,3}=load('EMT9_kdS_E_10000_10n.txt');x_all{1,4}=load('EMT9_kd3_E_10000_10n.txt');
x_all{1,5}=load('EMT9_kdzz_E_10000_10n.txt');x_all{1,6}=load('EMT9_kdZ_E_10000_10n.txt');
x_all{1,7}=load('EMT9_kd2_E_10000_10n.txt');x_all{1,8}=load('EMT9_kde_E_10000_10n.txt');
x_all{1,9}=load('EMT9_kdn_E_10000_10n.txt');
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\10n\H'
x_all{2,1}=load('EMT9_kdT_H_10000_10n.txt');x_all{2,2}=load('EMT9_kdss_H_10000_10n.txt');
x_all{2,3}=load('EMT9_kdS_H_10000_10n.txt');x_all{2,4}=load('EMT9_kd3_H_10000_10n.txt');
x_all{2,5}=load('EMT9_kdzz_H_10000_10n.txt');x_all{2,6}=load('EMT9_kdZ_H_10000_10n.txt');
x_all{2,7}=load('EMT9_kd2_H_10000_10n.txt');x_all{2,8}=load('EMT9_kde_H_10000_10n.txt');
x_all{2,9}=load('EMT9_kdn_H_10000_10n.txt');
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\10n\M'
x_all{3,1}=load('EMT9_kdT_M_10000_10n.txt');x_all{3,2}=load('EMT9_kdss_M_10000_10n.txt');
x_all{3,3}=load('EMT9_kdS_M_10000_10n.txt');x_all{3,4}=load('EMT9_kd3_M_10000_10n.txt');
x_all{3,5}=load('EMT9_kdzz_M_10000_10n.txt');x_all{3,6}=load('EMT9_kdZ_M_10000_10n.txt');
x_all{3,7}=load('EMT9_kd2_M_10000_10n.txt');x_all{3,8}=load('EMT9_kde_M_10000_10n.txt');
x_all{3,9}=load('EMT9_kdn_M_10000_10n.txt');
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_use'
i=3;
for j=1:3
    [x{i,j},R_dln{i,j},r_dln{i,j},J_dln{i,j}]=fun_xRrJ_dln(x_all(j,:),s0(j,2:end-1),p0);
end
for j=1:3
    idr=zeros(n*n,size(r_dln{i,j},3));idR=idr;idJ=idr;
    for k=1:size(idr,2)
        idr(:,k)=reshape(r_dln{i,j}(1:n,1:n,k)',n*n,1);
        idR(:,k)=reshape(R_dln{i,j}(1:n,1:n,k)',n*n,1);
        idJ(:,k)=reshape(J_dln{i,j}(1:n,1:n,k)',n*n,1);
    end
    r_one_dln{i,j}=idr; R_one_dln{i,j}=idR; J_one_dln{i,j}=idJ;
end

% 20n
x_all=cell(3,9);
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\20n\E'
x_all{1,1}=load('EMT9_kdT_E_10000_20n.txt');x_all{1,2}=load('EMT9_kdss_E_10000_20n.txt');
x_all{1,3}=load('EMT9_kdS_E_10000_20n.txt');x_all{1,4}=load('EMT9_kd3_E_10000_20n.txt');
x_all{1,5}=load('EMT9_kdzz_E_10000_20n.txt');x_all{1,6}=load('EMT9_kdZ_E_10000_20n.txt');
x_all{1,7}=load('EMT9_kd2_E_10000_20n.txt');x_all{1,8}=load('EMT9_kde_E_10000_20n.txt');
x_all{1,9}=load('EMT9_kdn_E_10000_20n.txt');
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\20n\H'
x_all{2,1}=load('EMT9_kdT_H_10000_20n.txt');x_all{2,2}=load('EMT9_kdss_H_10000_20n.txt');
x_all{2,3}=load('EMT9_kdS_H_10000_20n.txt');x_all{2,4}=load('EMT9_kd3_H_10000_20n.txt');
x_all{2,5}=load('EMT9_kdzz_H_10000_20n.txt');x_all{2,6}=load('EMT9_kdZ_H_10000_20n.txt');
x_all{2,7}=load('EMT9_kd2_H_10000_20n.txt');x_all{2,8}=load('EMT9_kde_H_10000_20n.txt');
x_all{2,9}=load('EMT9_kdn_H_10000_20n.txt');
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\20n\M'
x_all{3,1}=load('EMT9_kdT_M_10000_20n.txt');x_all{3,2}=load('EMT9_kdss_M_10000_20n.txt');
x_all{3,3}=load('EMT9_kdS_M_10000_20n.txt');x_all{3,4}=load('EMT9_kd3_M_10000_20n.txt');
x_all{3,5}=load('EMT9_kdzz_M_10000_20n.txt');x_all{3,6}=load('EMT9_kdZ_M_10000_20n.txt');
x_all{3,7}=load('EMT9_kd2_M_10000_20n.txt');x_all{3,8}=load('EMT9_kde_M_10000_20n.txt');
x_all{3,9}=load('EMT9_kdn_M_10000_20n.txt');
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_use'
i=4;
for j=1:3
    [x{i,j},R_dln{i,j},r_dln{i,j},J_dln{i,j}]=fun_xRrJ_dln(x_all(j,:),s0(j,2:end-1),p0);
end
for j=1:3
    idr=zeros(n*n,size(r_dln{i,j},3));idR=idr;idJ=idr;
    for k=1:size(idr,2)
        idr(:,k)=reshape(r_dln{i,j}(1:n,1:n,k)',n*n,1);
        idR(:,k)=reshape(R_dln{i,j}(1:n,1:n,k)',n*n,1);
        idJ(:,k)=reshape(J_dln{i,j}(1:n,1:n,k)',n*n,1);
    end
    r_one_dln{i,j}=idr; R_one_dln{i,j}=idR; J_one_dln{i,j}=idJ;
end
x0=x;r_one_dln0=r_one_dln;R_one_dln0=R_one_dln;J_one_dln0=J_one_dln;
save data_J_one_dln.mat J_one_dln
save data_rr_one_dln.mat r_one_dln 
save data_R_one_dln.mat R_one_dln
save data_x.mat x
clear i j k idJ idR idr x_all




