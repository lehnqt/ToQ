function [H0,Hu,Hc,O,psi0]=gen_H_kron(H_cert,H_unc,H_ctrl,obs,init_state,sys_num,sys_dim)
cert_num=length(H_cert);
ctrl_num=length(H_ctrl);
unc_num=length(H_unc);
H0=struct('op',cell(1,cert_num),'ft',cell(1,cert_num));
Hu=struct('op',cell(1,unc_num),'ft',cell(1,unc_num));
Hc=struct('op',cell(1,ctrl_num),'ft',cell(1,ctrl_num));
M=sys_dim;
for j=1:cert_num
H0(j).ft=H_cert(j).ft;
H0(j).op=sparse(prod(M),prod(M));
    for k=1:size(H_cert(j).op,1)
        H_temp=1;
        count=0;
        for l=1:sys_num
            if ismember(l,H_cert(j).sys)
                count=count+1;
            H_temp=kron(H_temp,H_cert(j).op{k,count});
            else
                H_temp=kron(H_temp,speye(M(l)));
            end
        end
        H0(j).op=H0(j).op+H_temp;
    end
        if ~issparse(H0(j).op)
        fprintf('warning: H0 %d is not sparse\n',j);
        end
end
for j=1:unc_num
Hu(j).ft=H_unc(j).ft;
Hu(j).op=sparse(prod(M),prod(M));
    for k=1:size(H_unc(j).op,1)
        H_temp=1;
        count=0;
        for l=1:sys_num
            if ismember(l,H_unc(j).sys)
                count=count+1;
            H_temp=kron(H_temp,H_unc(j).op{k,count});
            else
                H_temp=kron(H_temp,speye(M(l)));
            end
        end
        Hu(j).op=Hu(j).op+H_temp;
    end
        if ~issparse(Hu(j).op)
            fprintf('warning: Hu %d is not sparse\n',j);
        end
end
for j=1:ctrl_num
Hc(j).ft=H_ctrl(j).ft;
Hc(j).op=sparse(prod(M),prod(M));
    for k=1:size(H_ctrl(j).op,1)
        H_temp=1;
        count=0;
        for l=1:sys_num
            if ismember(l,H_ctrl(j).sys)
                count=count+1;
            H_temp=kron(H_temp,H_ctrl(j).op{k,count});
            else
                H_temp=kron(H_temp,speye(M(l)));
            end
        end
        Hc(j).op=Hc(j).op+H_temp;
    end
    if ~issparse(Hc(j).op)
       fprintf('warning: Hc %d is not sparse\n',j);
    end
end
obs_num=length(obs);
O=struct('op',cell(1,obs_num));
for j=1:obs_num
O(j).op=sparse(prod(M),prod(M));
    for k=1:size(obs(j).op,1)
        O_temp=1;
        count=0;
        for l=1:sys_num
            if ismember(l,obs(j).sys)
                count=count+1;
            O_temp=kron(O_temp,obs(j).op{k,count});
            else
                O_temp=kron(O_temp,speye(M(l)));
            end
        end
        O(j).op=O(j).op+O_temp;
    end
        if ~issparse(O(j).op)
        fprintf('warning: O %d is not sparse\n',j);
        end
end
if isempty(init_state)
    psi0=eye(prod(M));
else
psi0=0;
for j=1:size(init_state,1)
    psi_temp=1;
    for k=1:sys_num
        Id=eye(M(k));
        psi_temp=kron(psi_temp,Id(:,init_state(j,k)+1));
    end
    psi0=psi0+psi_temp;
end
if issparse(psi0)
       fprintf('warning: psi0 is not full\n');
end
end