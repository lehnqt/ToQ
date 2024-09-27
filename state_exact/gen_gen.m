function [sigma,kappa,zeta,tau,chi,gen_tot]=gen_gen(N)
sx=sparse([0 1;1 0]);
sy=sparse([0 -1i;1i,0]);
sz=sparse([1,0;0,-1]);
id=speye(2);
sigma=struct('sys',cell(N-1,1),'op',cell(N-1,1),'f',cell(N-1,1));
kappa=struct('sys',cell(N-1,1),'op',cell(N-1,1),'f',cell(N-1,1));
zeta=struct('sys',cell(N-1,1),'f',cell(N-1,1));
tau=struct('sys',cell(1,1),'op',cell(1,1),'f',cell(N-1,1));
chi=struct('sys',cell(1,1),'op',cell(1,1),'f',cell(N-1,1));
gen_tot=struct('sys',cell(3*N-1,1),'op',cell(3*N-1,1),'f',cell(3*N-1,1));
hid=cell(N,N);
for j=1:N
    for k=1:N
        hid{j,k}=id;
    end
end
%chi
h=hid;
for j=1:N  
  h{j,j}=sz;
end
count=1;
chi(count).sys=[1:N];
chi(count).op=h;
chi(count).f=1/(2*sqrt(2));
%tau
h=hid;
for j=1:N  
    for k=j:j+N-2
        h{j,mod(k-1,N)+1}=sz;
    end
end
count=1;
tau(count).sys=[1:N];
tau(count).op=h;
tau(count).f=1/(2*sqrt(2));
%sigma
for n=0:N-2
    h=hid;
    for j=1:N
         h{j,j}=sx;
         for k=j+1:j+n
              h{j,mod(k-1,N)+1}=sz;
         end
              k=j+n+1;
              h{j,mod(k-1,N)+1}=sx;
    end
    sigma(n+1).sys=[1:N];
    sigma(n+1).op=h;
    sigma(n+1).f=1/(2*sqrt(2));
end
%kappa
for n=0:N-2
    h=hid;
    for j=1:N
         h{j,j}=sy;
         for k=j+1:j+n
              h{j,mod(k-1,N)+1}=sz;
         end
              k=j+n+1;
              h{j,mod(k-1,N)+1}=sy;
    end
    kappa(n+1).sys=[1:N];
    kappa(n+1).op=h;
    kappa(n+1).f=1/(2*sqrt(2));
end
%zeta
for n=0:N-2
    h=[hid;hid];
    for j=1:N
         h{j,j}=sx;
         for k=j+1:j+n
              h{j,mod(k-1,N)+1}=sz;
         end
              k=j+n+1;
              h{j,mod(k-1,N)+1}=sy;
    end
    for j=1:N
         h{j+N,j}=sy;
         for k=j+1:j+n
              h{j+N,mod(k-1,N)+1}=sz;
         end
              k=j+n+1;
              h{j+N,mod(k-1,N)+1}=sx;
    end
    zeta(n+1).sys=[1:N];
    zeta(n+1).op=h;
    zeta(n+1).f=1/4;
end
for j=1:N-1
gen_tot(j).sys=sigma(j).sys;
gen_tot(j).op=sigma(j).op;
gen_tot(j).f=sigma(j).f;
end
for j=1:N-1
    gen_tot(j+N-1).sys=kappa(j).sys;
    gen_tot(j+N-1).op=kappa(j).op;
    gen_tot(j+N-1).f=kappa(j).f;
end
for j=1:N-1
    gen_tot(j+2*N-2).sys=zeta(j).sys;
    gen_tot(j+2*N-2).op=zeta(j).op;
    gen_tot(j+2*N-2).f=zeta(j).f;
end
    gen_tot(3*N-2).sys=tau(1).sys;
    gen_tot(3*N-2).op=tau(1).op;
    gen_tot(3*N-2).f=tau(1).f;
    gen_tot(3*N-1).sys=chi(1).sys;
    gen_tot(3*N-1).op=chi(1).op;
    gen_tot(3*N-1).f=chi(1).f;
end