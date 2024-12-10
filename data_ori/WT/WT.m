%% 1. save the unperturbed three state: E H M
clear;clc
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\WT'
tgf0=1;
k0t=0.06;kt=0.6;jt=0.06;nr2=2;
k0s=0.0006;ks=0.03;js=1.6;nt=2;
kS=17;jS=0.08;nr3=2;

k03=0.0012;k3=0.012;j13=0.15;j23=0.36;nS=2;nZ=2;

k0z=0.003;kz=0.06;jz=3.5;
kZ=17;jZ=0.06;
k02=0.0002;k2=0.012;j12=5;j22=0.2;
ke1=1;j1e=0.2;ke2=0.6;j2e=0.5;
kn1=1;j1n=0.2;kn2=0.6;j2n=0.5;

kdt=0.6;kds=0.09;kdS=1.66;kd3=0.03;kdz=0.09;kdZ=1.66;kd2=0.035;kde=0.5;kdn=0.5;

ntt=8000;dtt=0.1;

for j1=1:100
    T(1)=rand(1);s(1)=rand(1); S(1)=rand(1); R3(1)=rand(1); z(1)=rand(1);
    Z(1)=rand(1); R2(1)=rand(1); E(1)=rand(1); N(1)=rand(1);
    for i=1:ntt
        T(i+1)=T(i)+dtt*(k0t+kt/(1+(R2(i)/jt)^nr2)-kdt*T(i));
        s(i+1)=s(i)+dtt*(k0s+ks*(((T(i)+tgf0)/js)^nt)/(1+((T(i)+tgf0)/js)^nt)-kds*s(i));
        S(i+1)=S(i)+dtt*(kS*s(i)/(1+(R3(i)/jS)^nr3)-kdS*S(i));
        R3(i+1)=R3(i)+dtt*(k03+k3/(1+(S(i)/j13)^nS+(Z(i)/j23)^nZ)-kd3*R3(i));
        z(i+1)=z(i)+dtt*(k0z+kz*((S(i)/jz)^nS)/(1+(S(i)/jz)^nS)-kdz*z(i));
        Z(i+1)=Z(i)+dtt*(kZ*z(i)/(1+(R2(i)/jZ)^nr2)-kdZ*Z(i));
        R2(i+1)=R2(i)+dtt*(k02+k2/(1+(S(i)/j12)^nS+(Z(i)/j22)^nZ)-kd2*R2(i));
        E(i+1)=E(i)+dtt*(ke1/(1+(S(i)/j1e)^nS)+ke2/(1+(Z(i)/j2e)^nS)-kde*E(i));
        N(i+1)=N(i)+dtt*(kn1*((S(i)/j1n)^nS)/(1+(S(i)/j1n)^nS)+kn2*((Z(i)/j2n)^nZ)/(1+(Z(i)/j2n)^nZ)-kdn*N(i));
    end
    sol_all0(:,2)=T(:);sol_all0(:,3)=s(:);sol_all0(:,4)=S(:);sol_all0(:,5)=R3(:);sol_all0(:,6)=z(:);
    sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=0;sol_all0(:,11)=tgf0;
    while sum(abs(round((sol_all0(end, :) - sol_all0((end-500), :)), 4))) ~= 0
        T(1)=sol_all0(end,2);s(1)=sol_all0(end,3); S(1)=sol_all0(end,4); R3(1)=sol_all0(end,5); z(1)=sol_all0(end,6);
        Z(1)=sol_all0(end,7); R2(1)=sol_all0(end,8); E(1)=sol_all0(end,9); N(1)=sol_all0(end,10);
        for i=1:ntt
            T(i+1)=T(i)+dtt*(k0t+kt/(1+(R2(i)/jt)^nr2)-kdt*T(i));
            s(i+1)=s(i)+dtt*(k0s+ks*(((T(i)+tgf0)/js)^nt)/(1+((T(i)+tgf0)/js)^nt)-kds*s(i));
            S(i+1)=S(i)+dtt*(kS*s(i)/(1+(R3(i)/jS)^nr3)-kdS*S(i));
            R3(i+1)=R3(i)+dtt*(k03+k3/(1+(S(i)/j13)^nS+(Z(i)/j23)^nZ)-kd3*R3(i));
            z(i+1)=z(i)+dtt*(k0z+kz*((S(i)/jz)^nS)/(1+(S(i)/jz)^nS)-kdz*z(i));
            Z(i+1)=Z(i)+dtt*(kZ*z(i)/(1+(R2(i)/jZ)^nr2)-kdZ*Z(i));
            R2(i+1)=R2(i)+dtt*(k02+k2/(1+(S(i)/j12)^nS+(Z(i)/j22)^nZ)-kd2*R2(i));
            E(i+1)=E(i)+dtt*(ke1/(1+(S(i)/j1e)^nS)+ke2/(1+(Z(i)/j2e)^nS)-kde*E(i));
            N(i+1)=N(i)+dtt*(kn1*((S(i)/j1n)^nS)/(1+(S(i)/j1n)^nS)+kn2*((Z(i)/j2n)^nZ)/(1+(Z(i)/j2n)^nZ)-kdn*N(i));
        end
        sol_all0(:,2)=T(:);sol_all0(:,3)=s(:);sol_all0(:,4)=S(:);sol_all0(:,5)=R3(:);sol_all0(:,6)=z(:);
        sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=0;sol_all0(:,11)=tgf0;
    end
    sol_all(j1,:)=sol_all0(end,:);
end
sol_all=round(sol_all,4);sol_all=unique(sol_all,'rows');sol_all=sortrows(sol_all,2);
for ii=1:size(sol_all,1)-1
    if sum(abs(sol_all(ii,:)-sol_all(ii+1,:)))<0.01
        sol_all(ii,:)=[0,0,0,0,0,0,0,0,0,0,0];
    end
end
sol_all(all(sol_all==0,2),:)=[];
dlmwrite('EMT9_all1_EHM_k03_0012_kd3_03.txt',sol_all,'-append');


%% 2. save the Jacobian matrix
clear;clc;
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\WT'
tgf0=1;
k0t=0.06;kt=0.6;jt=0.06;nr2=2;
k0s=0.0006;ks=0.03;js=1.6;nt=2;
kS=17;jS=0.08;nr3=2;

k03=0.0012;k3=0.012;j13=0.15;j23=0.36;nS=2;nZ=2;

k0z=0.003;kz=0.06;jz=3.5;
kZ=17;jZ=0.06;
k02=0.0002;k2=0.012;j12=5;j22=0.2;
ke1=1;j1e=0.2;ke2=0.6;j2e=0.5;
kn1=1;j1n=0.2;kn2=0.6;j2n=0.5;

kdt=0.6;kds=0.09;kdS=1.66;kd3=0.03;kdz=0.09;kdZ=1.66;kd2=0.035;kde=0.5;kdn=0.5;
sol_all=ones(1,9);
syms T s S R3 z Z R2 E N
dT=k0t+kt/(1+(R2/jt)^nr2)-kdt*T;
ds=k0s+ks*(((T+tgf0)/js)^nt)/(1+((T+tgf0)/js)^nt)-kds*s;
dS=kS*s/(1+(R3/jS)^nr3)-kdS*S;
dR3=k03+k3/(1+(S/j13)^nS+(Z/j23)^nZ)-kd3*R3;
dz=k0z+kz*((S/jz)^nS)/(1+(S/jz)^nS)-kdz*z;
dZ=kZ*z/(1+(R2/jZ)^nr2)-kdZ*Z;
dR2=k02+k2/(1+(S/j12)^nS+(Z/j22)^nZ)-kd2*R2;
dE=ke1/(1+(S/j1e)^nS)+ke2/(1+(Z/j2e)^nZ)-kde*E;
dN=kn1*((S/j1n)^nS)/(1+(S/j1n)^nS)+kn2*((Z/j2n)^nZ)/(1+(Z/j2n)^nZ)-kdn*N;
J_sym_0=jacobian([dT,ds,dS,dR3,dz,dZ,dR2,dE,dN],[T,s,S,R3,z,Z,R2,E,N]);
save J_sym1_0.mat J_sym_0


