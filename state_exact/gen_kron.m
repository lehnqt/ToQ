function Hc=gen_kron(H_ctrl,sys_num,sys_dim)
ctrl_num=length(H_ctrl);
Hc=struct('op',cell(1,ctrl_num),'f',cell(1,ctrl_num));
M=sys_dim;
for j=1:ctrl_num
Hc(j).f=H_ctrl(j).f;
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
end