function [iF,iG,overlap]=infid_avg(H0,Hu,Hc,unc_tot,psi0,psi_tg,P,time_grid,c)
unc_num_tot=size(unc_tot,1);
cert_num=length(H0);
unc_num=length(Hu);
ctrl_num=length(Hc);
unc=unc_tot(:,1:unc_num);
unc_ctrl=unc_tot(:,unc_num+1:unc_num+ctrl_num);
bin_num=length(time_grid)-1;
d=size(psi_tg,2);
if d>1
    psi_tg=P*psi_tg*P';
else
    psi_tg=P*psi_tg;
end
c=reshape(c,[bin_num,ctrl_num]);
iF=zeros(1,unc_num_tot);
if nargout<2
parfor i=1:unc_num_tot
    unc_temp=unc(i,:);
    unc_ctrl_temp=unc_ctrl(i,:);
    psi=psi0;
for j=1:bin_num
    dt=time_grid(j+1)-time_grid(j);
    H_tot=sparse(0);
    for k=1:cert_num
        H_tot=H_tot+H0(k).op*H0(k).ft(time_grid(j)+dt/2);
    end
    for k=1:unc_num
        H_tot=H_tot+Hu(k).op*Hu(k).ft(time_grid(j)+dt/2)*unc_temp(k);
    end
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*Hc(k).ft(time_grid(j)+dt/2)*(1+unc_ctrl_temp(k))*c(j,k);
    end
    psi=expm(-1i*H_tot*dt)*psi;
end
    Ovl=trace(psi_tg'*psi)/d;
    iF(i)=1-abs(Ovl)^2;
end
else
iG=zeros(ctrl_num*bin_num,unc_num_tot);
overlap=zeros(bin_num,unc_num_tot);
parfor i=1:unc_num_tot
    unc_temp=unc(i,:);
    unc_ctrl_temp=unc_ctrl(i,:);
    iG_temp=zeros(bin_num,ctrl_num);
    overlap_temp=zeros(bin_num,1);
    D=cell(bin_num,ctrl_num);
    psi_fw=cell(1,bin_num+1);
    psi_bw=cell(1,bin_num+1);
    psi_fw{1}=psi0;
    psi_bw{1}=psi_tg;
for j=1:bin_num
    %forward propagation
    dt=time_grid(j+1)-time_grid(j);
    H_tot=sparse(0);
    for k=1:cert_num
        H_tot=H_tot+H0(k).op*H0(k).ft(time_grid(j)+dt/2);
    end
    for k=1:unc_num
        H_tot=H_tot+Hu(k).op*Hu(k).ft(time_grid(j)+dt/2)*unc_temp(k);
    end
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*Hc(k).ft(time_grid(j)+dt/2)*(1+unc_ctrl_temp(k))*c(j,k);
    end
    psi_fw{j+1}=expm(-1i*H_tot*dt)*psi_fw{j};
    for k=1:ctrl_num
        C1=Hc(k).op*Hc(k).ft(time_grid(j)+dt/2)*(1+unc_ctrl_temp(k));
        C2=H_tot*C1-C1*H_tot;
        D{j,k}=-1i*dt*C1+(dt^2/2)*C2;
%         C3=(H0{i}+H_int)*C2-C2*(H0{i}+H_int);
%         D{k,j}=-1i*Delta_t*C1+(Delta_t^2/2)*C2+1i*(Delta_t^3/6)*C3;
    end
    %backward propagation
    dt=time_grid(bin_num+2-j)-time_grid(bin_num+1-j);
    H_tot=sparse(0);
    for k=1:cert_num
        H_tot=H_tot+H0(k).op*H0(k).ft(time_grid(bin_num+1-j)+dt/2);
    end
    for k=1:unc_num
        H_tot=H_tot+Hu(k).op*Hu(k).ft(time_grid(bin_num+1-j)+dt/2)*unc_temp(k);
    end
    for k=1:ctrl_num
        H_tot=H_tot+Hc(k).op*Hc(k).ft(time_grid(bin_num+1-j)+dt/2)*(1+unc_ctrl_temp(k))*c(bin_num+1-j,k);
    end
    psi_bw{j+1}=expm(1i*H_tot*dt)*psi_bw{j};
end
    Ovl=trace(psi_bw{1}'*psi_fw{bin_num+1})/d;
    iF(i)=1-abs(Ovl)^2;
for j=1:bin_num
    overlap_temp(j)=trace((psi_bw{bin_num-j+2})'*psi_fw{j})/d;
    for k=1:ctrl_num
        iG_temp(j,k)=-(2/d)*real(trace(psi_bw{bin_num-j+2}'*D{j,k}*psi_fw{j})*conj(Ovl));
    end
end
iG(:,i)=iG_temp(:);
overlap(:,i)=overlap_temp;
end
iG=mean(iG,2);
end
iF=mean(iF);
end