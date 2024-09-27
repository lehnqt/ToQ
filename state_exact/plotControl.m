function plotControl(varargin)
col=linspecer(4);
% inputs:

% ct discrete control, a Nt by 2 real array OR a Nt  by 1 complex array

% T maximum gate time

% titleS the title of the graphs

% example:

% t=linspace(0,12)'; plotControl(exp(1i*t),25)

% t=linspace(0,5)'; plotControl([t/5,exp(-t.^2)],5,'Gauss')

 

ct=varargin{1};

T=1;

titleS='';

for ii=2:nargin

if isnumeric(varargin{ii})

    T=varargin{ii};

elseif ischar(varargin{ii})

    titleS=varargin{ii};

end

end

 

Nt=size(ct,1);

if size(ct,2)==2

    %Nt=Nt*2;

    ct=[ct(:,1)+ 1i*ct(:,2)];

end

 

figure;

t=linspace(0,T,Nt);

plot(t,[real(ct), imag(ct)],'lineWidth',2);

title([titleS,'Control continuous'])

ylim([-1,1]*max(abs(ct))*1.05)

xlabel('Time')

ylabel('amplitude')

legend('real', 'imaginary')

 

figure

subplot(2,1,1);

plot(t,abs(ct),'lineWidth',2);

title([titleS,' Control Spherical Coordinates'])

ylim([-abs(min(abs(ct)))*1.05,max(abs(ct))*1.05])

ylabel('Amplitude')

subplot(2,1,2)

plot(t,angle(ct),'lineWidth',2)

ylabel('Angle (rad)')

 

figure;

dt=T/Nt;

t=linspace(0,T,Nt)+dt/2;

b=bar(t,real(ct),1,'EdgeColor','none','FaceColor',col(4,:));

b.FaceAlpha=0.5;

hold on

b2=bar(t,imag(ct),1,'EdgeColor','none','FaceColor',col(3,:));

b2.FaceAlpha=0.5;

% title([titleS, ' Control discrete'])

legend('$h_1$','$h_2$','interpreter','latex')

% xlabel('Time')

% ylabel('amplitude')
set(gca,'fontsize',14);
box on;
set(gcf,'color','w');
hold off;
end