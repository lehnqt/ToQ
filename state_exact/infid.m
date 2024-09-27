function [iF,c]=infid(M0,Mf,c0,ctg,time_grid,f)
bin_num=length(time_grid)-1;
const_num=length(M0);
ctrl_num=length(Mf);
f=reshape(f,[bin_num,ctrl_num]);
c=c0;
for j=1:bin_num
    dt=time_grid(j+1)-time_grid(j);
    M_tot=sparse(0);
    for k=1:const_num
        M_tot=M_tot+M0(k).op*M0(k).ft(time_grid(j)+dt/2);
    end
    for k=1:ctrl_num
        M_tot=M_tot+Mf(k).op*Mf(k).ft(time_grid(j)+dt/2)*f(j,k);
    end
    c=expm(-1i*M_tot*dt)*c;
end
    Ovl=c'*ctg;
    iF=1-Ovl;
end