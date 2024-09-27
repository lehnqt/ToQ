function psi=state_evol_lanczos(H0,Hu,Hc,unc_tot,psi0,time_grid,c,numK)
unc_num_tot=size(unc_tot,1);
cert_num=length(H0);
unc_num=length(Hu);
ctrl_num=length(Hc);
unc=unc_tot(:,1:unc_num);
unc_ctrl=unc_tot(:,unc_num+1:unc_num+ctrl_num);
bin_num=length(time_grid)-1;
c=reshape(c,[bin_num,ctrl_num]);
if nargout<2
for i=1:unc_num_tot
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
    psi=expvcpu(dt,H_tot,psi,numK);
end
end
end
end