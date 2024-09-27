function [iF,maxD]=mps_infid_T_nograd(H0,Hc,x,mps0,mpstg,sv_min,D,midstep,nt)
c=x(1:end-1);
T=x(end);
n=length(mps0);
nc=length(Hc);
nbin=length(c)/(nc);
c=reshape(c,[nbin,nc]);
d=size(mps0{1},2);
dt=T/(nbin*nt);
mps=mps0;
maxD=1;
g2=cell(1,n-1);
for j=1:n-1
        h=H0{j};
        gate=expm(-1i*dt*h);
        gate=reshape(gate,[d,d,d,d]);
        g2{j}=gate;
end
for k=1:nbin
     %gate construction
   g1=cell(1,n);
    for j=1:n
       h=zeros(d);
       for jc=1:nc
            for js=1:length(Hc(jc).sys)
                if Hc(jc).sys(js)==j
        h=h+c(k,jc)*Hc(jc).op{js};
                end
            end
       end
       gate=expm(-1i*dt*h);
       g1{j}=gate;
    end 
%forward propagation
for jt=1:nt
%apply 1q terms
for j=1:n
    [mps{j}]=mps_gate_1q(mps{j},g1{j});
end
 %apply odd terms
for j=1:2:n-1
    [mps{j},mps{j+1}]=mps_gate_2q(mps{j},mps{j+1},g2{j},sv_min,D);
end
 %apply even terms
for j=2:2:n-1
    [mps{j},mps{j+1}]=mps_gate_2q(mps{j},mps{j+1},g2{j},sv_min,D);
end
end
if mod(k-1,midstep)==0||k==nbin
mps=mps_normalize(mps);
end
for j=1:n
maxD=max(maxD,max(size(mps{j})));
end
end
ovl=mps_overlap(mpstg,mps);
iF=1-abs(ovl)^2;   
end
