function [x,R,r,J]=fun_xRrJ_dln(x_all,s0,p0)

n=size(x_all,2);
x_num=zeros(1,n);
for i=1:n
    x_num(1,i)=size(x_all{1,i},1);
end
num_min=min(x_num);
x_rand=zeros(num_min,n);
for i=1:n
    x_rand(:,i)=randi(x_num(1,i),num_min,1);
end
% get x
x=zeros(n,n+1,num_min);
for i=1:n
    x(i,:,:)=x_all{1,i}(x_rand(:,i),[2:n+1,end-1])';
end
id=all((abs(x(:,end,:)-p0(1:n)')>0.0001*ones(n,1))==1,1);
x=x(:,:,id);
% calculate R
num_k=size(x,3);
R=zeros(n,n,num_k);
for k=1:num_k
    for i=1:n
        R(:,i,k)=((x(i,1:n,k)-s0(1:n))./(x(i,1:n,k)+s0(1:n)).*(x(i,end,k)+p0(i))./(x(i,end,k)-p0(i)))';
    end
end
% remove det(R)=0
id=zeros(num_k,1); % the det of R
for i=1:num_k
    id(i,1)=det(R(:,:,i));
end
dis=find(abs(id)<10^(-14));
x(:,:,dis)=[]; R(:,:,dis)=[];
% remove det(diag(inv(R)))=0
num_k=size(x,3);
R_inv=zeros(n,n,num_k);
R_inv_diag=zeros(n,n,num_k);
id=zeros(num_k,1);
for k=1:num_k
    R_inv(:,:,k)=inv(R(:,:,k));
    R_inv_diag(:,:,k)=diag(diag(R_inv(:,:,k)));
    id(k,1)=det(R_inv_diag(:,:,k));
end
dis=find(abs(id)<10^(-14));
x(:,:,dis)=[]; R(:,:,dis)=[]; R_inv(:,:,dis)=[]; R_inv_diag(:,:,dis)=[];
% calculat r J
num_k=size(x,3);
r=zeros(n,n,num_k);
J=zeros(n,n,num_k);
for i=1:num_k
    r(:,:,i)=-inv(R_inv_diag(:,:,i))*R_inv(:,:,i);
    J(:,:,i)=diag(s0(1:n))*diag(p0(1:n))*R_inv(:,:,i)*inv(diag(s0(1:n)));
end
r_mean=mean(r,3); r_std=std(r,[],3);
sig=3;
for j=1:n
    for i=1:n
        id=find(r(i,j,:)>(r_mean(i,j)+sig*r_std(i,j)) | ...
            r(i,j,:)<(r_mean(i,j)-sig*r_std(i,j)));
        x(:,:,id)=[]; R(:,:,id)=[]; R_inv(:,:,id)=[]; R_inv_diag(:,:,id)=[]; r(:,:,id)=[]; J(:,:,id)=[];
    end
end

end