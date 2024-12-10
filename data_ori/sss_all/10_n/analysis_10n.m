%% 10n

%% 1. perturbation paramaters
clear
clc
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all\10n'

p_all=20000;
p=[0.6 0.09 1.66 0.03 0.09 1.66 0.035 0.5 0.5];
para=zeros(p_all,size(p,2));
for i=1:p_all
    for j=1:size(p,2)
        para(i,j)=p(1,j)+normrnd(0,0.1*p(1,j));
    end
end
m=mean(para,1);
s=std(para,1);
sig=1.96;
for j=1:size(para,2)
    id=find( para(:,j)>(m(j)+sig*s(j)) | para(:,j)<(m(j)-sig*s(j)) );
    para(id,:)=[];
end
dlmwrite('EMT9_para_10000_10n.txt',para(1:10000,:),'-append');


%% 2. calculate the perturbed sss
% 1 kdt
clear;clc
para=load('EMT9_para_10000_10n.txt');
ntt=4000;dtt=0.1;
for i1=1:length(para)
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

    kdt=para(i1,1);kds=0.09;kdS=1.66;kd3=0.03;kdz=0.09;kdZ=1.66;kd2=0.035;kde=0.5;kdn=0.5;

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
        sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kdt;
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
            sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kdt;
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
    dlmwrite('EMT9_kdt_all_10000_10n.txt',sol_all,'-append');
end

% 2 kds
clear;clc
para=load('EMT9_para_10000_10n.txt');
ntt=4000;dtt=0.1;
for i1=1:length(para)
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

    kdt=0.6;kds=para(i1,2);kdS=1.66;kd3=0.03;kdz=0.09;kdZ=1.66;kd2=0.035;kde=0.5;kdn=0.5;

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
        sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kds;
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
            sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kds;
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
    dlmwrite('EMT9_kdss_all_10000_10n.txt',sol_all,'-append');
end

% 3 kdS
clear;clc
para=load('EMT9_para_10000_10n.txt');
ntt=4000;dtt=0.1;
for i1=1:length(para)
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

    kdt=0.6;kds=0.09;kdS=para(i1,3);kd3=0.03;kdz=0.09;kdZ=1.66;kd2=0.035;kde=0.5;kdn=0.5;

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
        sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kdS;
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
            sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kdS;
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
    dlmwrite('EMT9_kdS_all_10000_10n.txt',sol_all,'-append');
end

% 4 kd3
clear;clc
para=load('EMT9_para_10000_10n.txt');
ntt=4000;dtt=0.1;
for i1=1:length(para)
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

    kdt=0.6;kds=0.09;kdS=1.66;kd3=para(i1,4);kdz=0.09;kdZ=1.66;kd2=0.035;kde=0.5;kdn=0.5;

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
        sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kd3;
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
            sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kd3;
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
    dlmwrite('EMT9_kd3_all_10000_10n.txt',sol_all,'-append');
end

% 5 kdz
clear;clc
para=load('EMT9_para_10000_10n.txt');
ntt=4000;dtt=0.1;
for i1=1:length(para)
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

    kdt=0.6;kds=0.09;kdS=1.66;kd3=0.03;kdz=para(i1,5);kdZ=1.66;kd2=0.035;kde=0.5;kdn=0.5;

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
        sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kdz;
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
            sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kdz;
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
    dlmwrite('EMT9_kdzz_all_10000_10n.txt',sol_all,'-append');
end

% 6 kdZ
clear;clc
para=load('EMT9_para_10000_10n.txt');
ntt=4000;dtt=0.1;
for i1=1:length(para)
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

    kdt=0.6;kds=0.09;kdS=1.66;kd3=0.03;kdz=0.09;kdZ=para(i1,6);kd2=0.035;kde=0.5;kdn=0.5;

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
        sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kdZ;
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
            sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kdZ;
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
    dlmwrite('EMT9_kdZ_all_10000_10n.txt',sol_all,'-append');
end

% 7 kd2
clear;clc
para=load('EMT9_para_10000_10n.txt');
ntt=4000;dtt=0.1;
for i1=1:length(para)
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

    kdt=0.6;kds=0.09;kdS=1.66;kd3=0.03;kdz=0.09;kdZ=1.66;kd2=para(i1,7);kde=0.5;kdn=0.5;

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
        sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kd2;
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
            sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kd2;
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
    dlmwrite('EMT9_kd2_all_10000_10n.txt',sol_all,'-append');
end

% 8 kde
clear;clc
para=load('EMT9_para_10000_10n.txt');
ntt=4000;dtt=0.1;
for i1=1:length(para)
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

    kdt=0.6;kds=0.09;kdS=1.66;kd3=0.03;kdz=0.09;kdZ=1.66;kd2=0.035;kde=para(i1,8);kdn=0.5;

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
        sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kde;
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
            sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kde;
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
    dlmwrite('EMT9_kde_all_10000_10n.txt',sol_all,'-append');
end
% 9 kdn
clear;clc
para=load('EMT9_para_10000_10n.txt');
ntt=4000;dtt=0.1;
for i1=1:length(para)
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

    kdt=0.6;kds=0.09;kdS=1.66;kd3=0.03;kdz=0.09;kdZ=1.66;kd2=0.035;kde=0.5;kdn=para(i1,9);

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
        sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kdn;
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
            sol_all0(:,7)=Z(:);sol_all0(:,8)=R2(:);sol_all0(:,9)=E(:);sol_all0(:,10)=N(:);sol_all0(:,1)=i1;sol_all0(:,11)=kdn;
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
    dlmwrite('EMT9_kdn_all_10000_10n.txt',sol_all,'-append');
end


%% 3. classify the three sss for these parameters to E H M file

clear;clc;
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all_10n'
sss_all=load('EMT9_kdt_all_10000_10n.txt'); % change, for the nine sss file
sss_all0=sss_all;
sss_array=tabulate(sss_all(:,1));sss_array=sortrows(sss_array,-2);
sss_mean = mean(sss_all(:,2:10));sss_std = std(sss_all(:,2:10));maxs=[max(sss_all(:,2:10))];% 3sigma处理数据
min3si=sss_mean-3*sss_std;max3si=sss_mean+3*sss_std;
for j=2:10
     num=find(sss_all(:,j)>max3si(j-1) | sss_all(:,j)<min3si(j-1));
     sss_all(num,:)=[];
end
sss_nor=zscore(sss_all(:,2:10));
y=pdist(sss_nor,'euclidean');z=linkage(y,'average');
TT=cluster(z,'maxclust',3); % modified, 4 maybe 3 2 
sss_all(:,12)=TT;
phen1=find(TT==1);sss_phen1=sss_all(phen1,:);
phen2=find(TT==2);sss_phen2=sss_all(phen2,:);
phen3=find(TT==3);sss_phen3=sss_all(phen3,:);
% phen4=find(TT==4);sss_phen4=sss_all(phen4,:);
% phen5=find(TT==5);sss_phen5=sss_all(phen5,:);
[coeff,score,latent,tsquared,explained,mu]=pca(sss_nor);
sss_all(:,12)=sss_all(:,1);
snum=tabulate(sss_all(:,1)); % sss in every paramaters
snum1=find(snum(:,2)==1); % the paramaters when sss number=1
for i=1:length(snum1)
    s1=sum(snum(1:snum1(i)-1,2))+1;sss_num1(i,:)=sss_all(s1,:); % sss when sss number=1
end
snum2=find(snum(:,2)==2); % the paramaters when sss number=2
for i=1:length(snum2)
    s2=sum(snum(1:snum2(i)-1,2))+1;sss_num2(2*i-1,:)=sss_all(s2,:);sss_num2(2*i,:)=sss_all(s2+1,:); % sss when sss number=2
end
snum3=find(snum(:,2)==3); % the paramaters when sss number=3
for i=1:length(snum3)
    s3=sum(snum(1:snum3(i)-1,2))+1;sss_num3(3*i-2,:)=sss_all(s3,:);
    sss_num3(3*i-1,:)=sss_all(s3+1,:);sss_num3(3*i,:)=sss_all(s3+2,:); % sss when sss number=3
end
phen=[length(phen1),length(phen2),length(phen3)];
phenfre=phen./sum(phen);
[sss_phen1(1,2:10);sss_phen2(1,2:10);sss_phen3(1,2:10)] 

% % based on some information (figure) to determine the E H M file; for example
% cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all_10n\E'
% dlmwrite('EMT9_kdT_E_10000_10n.txt',sss_phen1,'-append');
% cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all_10n\H'
% dlmwrite('EMT9_kdT_H_10000_10n.txt',sss_phen2,'-append');
% cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_ori\sss_all_10n\M'
% dlmwrite('EMT9_kdT_M_10000_10n.txt',sss_phen3,'-append');
