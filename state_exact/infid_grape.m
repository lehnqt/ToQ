function [iF,iG,c]=infid_grape(M0,Mf,c0,ctg,time_grid,f)
bin_num=length(time_grid)-1;
const_num=length(M0);
ctrl_num=length(Mf);
f=reshape(f,[bin_num,ctrl_num]);
if nargout<2
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
    ovl=c'*ctg;
    iF=1-ovl;
else
    iG=zeros(bin_num,ctrl_num);
    D=cell(bin_num,ctrl_num);
    cfw=cell(1,bin_num+1);
    cbw=cell(1,bin_num+1);
    cfw{1}=c0;
    cbw{1}=ctg;
    for j=1:bin_num
            dt=time_grid(j+1)-time_grid(j);
            M_tot=sparse(0);
        for k=1:const_num
            M_tot=M_tot+M0(k).op*M0(k).ft(time_grid(j)+dt/2);
        end
        for k=1:ctrl_num
            M_tot=M_tot+Mf(k).op*Mf(k).ft(time_grid(j)+dt/2)*f(j,k);
        end
        %forward propagation
        cfw{j+1}=expm(-1i*M_tot*dt)*cfw{j};
        %derivative matrix
        for k=1:ctrl_num
            C1=Mf(k).op*Mf(k).ft(time_grid(j)+dt/2);
            C2=M_tot*C1-C1*M_tot;
            D{j,k}=-1i*dt*C1+(dt^2/2)*C2;
        end
        %backward propagation
        dt=time_grid(bin_num+2-j)-time_grid(bin_num+1-j);
        M_tot=sparse(0);
        for k=1:const_num
            M_tot=M_tot+M0(k).op*M0(k).ft(time_grid(bin_num+1-j)+dt/2);
        end
        for k=1:ctrl_num
            M_tot=M_tot+Mf(k).op*Mf(k).ft(time_grid(bin_num+1-j)+dt/2)*f(bin_num+1-j,k);
        end
        cbw{j+1}=expm(1i*M_tot*dt)*cbw{j};
    end
    ovl=cbw{1}'*cfw{bin_num+1};
    iF=1-ovl;
    for j=1:bin_num
        for k=1:ctrl_num
            iG(j,k)=-cbw{bin_num-j+2}'*D{j,k}*cfw{j};
        end
    end
    iG=iG(:);
    c=cfw{bin_num+1};
end


