%Extract epsilon and mu of a thin polymer box filled with Fe3O4 powder
%The cap is at port one

clear all;

[fre,s11_m,s11_p,s21_m,s21_p,s12_m,s12_p,s22_m,s22_p]=textread('1.s2p','%f%f%f%f%f%f%f%f%f','headerlines',8);

% for i=1:800
% if s11_p(i)<0
% s11_p(i)=s11_p(i)+360;
% end
% if s12_p(i)<0
% s12_p(i)=s12_p(i)+360;
% end
% if s21_p(i)<0
% s21_p(i)=s21_p(i)+360;
% end
% if s22_p(i)<0
% s22_p(i)=s22_p(i)+360;
% end
% end


TEaLs11p=10.^(s11_m/20);
TEaLs21p=10.^(s21_m/20);


TEas11a=s11_p;
TEas21a=s21_p;


TEaRs11a=TEas11a*pi/180;
TEaRs21a=TEas21a*pi/180;

% figure;
% plot(fre,TEaRs11a,'r',fre,TEaRs21a,'b');
% legend('Rs11-arg:rad','Rs21-arg:rad');
% grid

TEas11=TEaLs11p.*(cos(TEaRs11a)-1i*sin(TEaRs11a));
TEas21=TEaLs21p.*(cos(TEaRs21a)-1i*sin(TEaRs21a));

%%%%%%%%%%%%%%%%%%%%%%%%%%
Cv=3*10^8;
wa=0.02286;%WR-90
kx=pi/wa;
k0=2*pi*fre/Cv;
ka=sqrt(k0.^2-kx^2);
da=0.004;%thick of the slab under test  被测量样品厚度
d1a=0; %length between slab with detector of S11 靠近port1的距离
t14=0.00968;%1/4 wavelength waveguide thickness   
dta=t14-da; %length between slab with detector of S21   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ra=TEas11.*exp(-2i*d1a*ka);

ta=TEas21.*exp(-1i*(d1a+dta)*ka);

z1a=sqrt(((1+ra).^2-ta.^2)./((1-ra).^2-ta.^2));
z2a=-sqrt(((1+ra).^2-ta.^2)./((1-ra).^2-ta.^2));

delta=0.01;

ephi1a=abs(ta./(1-ra.*(z1a-1)./(z1a+1)))<=1;
ephi2a=abs(ta./(1-ra.*(z2a-1)./(z2a+1)))<=1;

zma=(abs(real(z1a))>=delta).*(z1a.*sign(real(z1a)))+(abs(real(z1a))<delta).*(ephi1a.*z1a+ephi2a.*z2a);

exxa=log(ta./(1-ra.*(zma-1)./(zma+1)));
nma=(imag(exxa)-1i.*real(exxa))./(ka*da);

dna=2*pi./(ka*da);

% figure;
% plot(fre,real(zma),'r-',fre,imag(zma),'b-');
% legend('real(zma)','imag(zma)');
% grid
% 
% figure;
% plot(fre,real(nma),'r-',fre,imag(nma),'b-');
% legend('real(nma)','imag(nma)');
% grid

% axis([1.6 2.4 -5 5]);

na1p=dna+nma;
na1n=-dna+nma;
na2p=dna*2+nma;
na2n=-dna*2+nma;
na3p=dna*3+nma;
na3n=-dna*3+nma;
na4p=dna*4+nma;
na4n=-dna*4+nma;
%nma 是没有额外波长的长度（样品比较小），加减dna是有额外的波长（样品比较大），一个一个试看哪个对
nnnna=nma;%choose phase nma, na1p, na1n
% nnnna(1:98)=na1p(1:98);%for castor oil
% nnnna(99:end)=na1p(99:end);%for mineral oil
% nnnna=zeros(size(nma));
% if n==1
%     nnnna(1:71)=nma(1:71);%choose phase
%     nnnna(72:end)=na1p(72:end);%choose phase
% else
%     nnnna(1:75)=nma(1:75);%choose phase
%     nnnna(76:end)=na1p(76:end);%choose phase
% end

% figure;
% plot(fre,nma,'b',fre,na1p,'b',fre,na2p,'b',fre,na1n,'b',fre,na2n,'b',fre,na3p,'b',fre,na3n,'b',fre, nnnna,'k',fre,-dna/2,'r',fre,dna/2,'r');

% axis([1.6 2.4 -50 50]);
% grid;

uuuua=nnnna.*zma;

%----------------------TEb-------------------------------------------------

TEbLs11p=10.^(s11_m/20);
TEbLs21p=10.^(s21_m/20);

% figure;
% plot(fre,TEbLs11p,'r',fre,TEbLs21p,'b');
% legend('Ls11-linear','Ls21-linear');
% grid

TEbs11a=s11_p;
TEbs21a=s21_p;

TEbRs11a=TEbs11a*pi/180;
TEbRs21a=TEbs21a*pi/180;

% figure;
% plot(fre,TEbRs11a,'r',fre,TEbRs21a,'b');
% legend('Rs11-arg:rad','Rs21-arg:rad');
% grid

TEbs11=TEbLs11p.*(cos(TEbRs11a)-1i*sin(TEbRs11a));
TEbs21=TEbLs21p.*(cos(TEbRs21a)-1i*sin(TEbRs21a));

%%%%%%%%%%%%%%%%%%%%%%%%%%
kb=sqrt(k0.^2-kx^2);
db=da;  %thick of the slab
d1b=d1a; %length between slab with detector of S11
dtb=dta; %length between slab with detector of S21
rb=TEbs11.*exp(-2i*d1b*kb);

tb=TEbs21.*exp(-1i*(d1b+dtb)*kb);

z1b=sqrt(((1+rb).^2-tb.^2)./((1-rb).^2-tb.^2));
z2b=-sqrt(((1+rb).^2-tb.^2)./((1-rb).^2-tb.^2));

deltb=0.01;

ephi1b=abs(tb./(1-rb.*(z1b-1)./(z1b+1)))<=1;
ephi2b=abs(tb./(1-rb.*(z2b-1)./(z2b+1)))<=1;

zmb=(abs(real(z1b))>=deltb).*(z1b.*sign(real(z1b)))+(abs(real(z1b))<deltb).*(ephi1b.*z1b+ephi2b.*z2b);

exxb=log(tb./(1-rb.*(zmb-1)./(zmb+1)));
nmb=(imag(exxb)-1i.*real(exxb))./(kb*db);

dnb=2*pi./(kb*db);

% figure;
% plot(fre,real(zmb),'r-',fre,imag(zmb),'b-');
% legend('real(zmb)','imag(zmb)');
% grid

% figure;
% plot(fre,real(nmb),'r-',fre,imag(nmb),'b-');
% legend('real(nmb)','imag(nmb)');
% grid

% axis([1.6 2.4 -5 5]);

nb1p=dnb+nmb;
nb1n=-dnb+nmb;
nb2p=dnb*2+nmb;
nb2n=-dnb*2+nmb;
nb3p=dnb*3+nmb;
nb3n=-dnb*3+nmb;
nb4p=dnb*4+nmb;
nb4n=-dnb*4+nmb;
%nnnnb考虑了各向异性，这里没用

nnnnb=nmb;%choose phase
% nnnnb(58:end)=nb1p(58:end);%for castor oil
% nnnnb(99:end)=nb1p(99:end);%for mineral oil
% nnnnb=zeros(size(nmb));
% if n==1
%     nnnnb(1:71)=nmb(1:71);%choose phase
%     nnnnb(72:end)=nb1p(72:end);%choose phase
% else
%     nnnnb(1:75)=nmb(1:75);%choose phase
%     nnnnb(76:end)=nb1p(76:end);%choose phase
% end

uuuub=nnnnb.*zmb;

eeeea=(nnnna.^2.*ka.^2+uuuua./uuuub.*kx^2)./k0.^2./uuuua;
eeeeb=(nnnnb.^2.*kb.^2+uuuub./uuuua.*kx^2)./k0.^2./uuuub;

if isempty(find(abs(diff(real(uuuua)))>=0.1, 1))==0%discontinuous
    D=find(abs(diff(real(uuuua)))>=0.1);
    uuuua(1:D)=(nnnna(1:D)-dnb(1:D)).*zma(1:D);
    uuuub(1:D)=(nnnnb(1:D)-dnb(1:D)).*zmb(1:D);
    eeeea(1:D)=((nnnna(1:D)-dnb(1:D)).^2.*ka(1:D).^2+uuuua(1:D)./uuuub(1:D).*kx^2)./k0(1:D).^2./uuuua(1:D);
    eeeeb(1:D)=((nnnnb(1:D)-dnb(1:D)).^2.*kb(1:D).^2+uuuub(1:D)./uuuua(1:D).*kx^2)./k0(1:D).^2./uuuub(1:D);
end

figure
plot(fre/1e9,s11_m,'r',fre/1e9,s21_m,'b',fre/1e9,s12_m,'k',fre/1e9,s22_m,'g','Linewidth',1.7);grid
hold on
legend('mag(S_1_1)','mag(S_2_1)','mag(S_1_2)','mag(S_2_2)');
title('Measured S-Parameter Magnitudes','Fontsize',14,'Fontname','Helvetica');
xlabel('Frequency (GHz)','Fontsize',14,'Fontname','Helvetica');
ylabel('dB','Fontsize',14,'Fontname','Helvetica');
set(gca,'Fontsize',14,'Fontname','Helvetica');
%xlim([8 12]);
saveas(gcf,'CastorOil_Type4_Sp_Mag_New.fig');

figure
plot(fre/1e9,s11_p,'r',fre/1e9,s21_p,'b',fre/1e9,s12_p,'k',fre/1e9,s22_p,'g','Linewidth',1.7);grid
hold on
legend('ang(S_1_1)','ang(S_2_1)','ang(S_1_2)','ang(S_2_2)');
title('Measured S-Parameter Phases','Fontsize',14,'Fontname','Helvetica');
xlabel('Frequency (GHz)','Fontsize',14,'Fontname','Helvetica');
ylabel('degree','Fontsize',14,'Fontname','Helvetica');
set(gca,'Fontsize',14,'Fontname','Helvetica');
%xlim([8 12]);
saveas(gcf,'CastorOil_Type4_Sp_Ph_New.fig')

figure
plot(fre/1e9,real(eeeea),'r',fre/1e9,imag(eeeea),'b','Linewidth',1.7);grid
hold on
legend('Re\{\epsilon\}','Im\{\epsilon\}');
title('Extracted Permittivity','Fontsize',14,'Fontname','Helvetica');
xlabel('Frequency (GHz)','Fontsize',14,'Fontname','Helvetica');
set(gca,'Fontsize',14,'Fontname','Helvetica');
%xlim([8 12]);
saveas(gcf,'CastorOil_Type4_Sp_E_New.fig')

figure
plot(fre/1e9,real(uuuua),'r',fre/1e9,imag(uuuua),'b','Linewidth',1.7);grid
hold on
legend('Re\{\mu\}','Im\{\mu\}');
title('Extracted Permeability','Fontsize',14,'Fontname','Helvetica');
xlabel('Frequency (GHz)','Fontsize',14,'Fontname','Helvetica');
set(gca,'Fontsize',14,'Fontname','Helvetica');
%xlim([8 12]);
saveas(gcf,'CastorOil_Type4_Sp_U_New.fig')

figure
plot(fre/1e9,imag(eeeea)./real(eeeea),'Linewidth',1.7);grid
title('Calculated Dielectric Loss Tangent','Fontsize',14,'Fontname','Helvetica');
xlabel('Frequency (GHz)','Fontsize',14,'Fontname','Helvetica');
set(gca,'Fontsize',14,'Fontname','Helvetica');
%xlim([8 12]);
saveas(gcf,'CastorOil_Type4_Sp_Loss_New.fig')


%Fit the extract ep and mu into a line using Least Squares Polynomial Fit
%Order n=1, meaning a linear fit
c=polyfit((1:numel(eeeea))'/numel(eeeea),real(eeeea),1);
fitted_epr=c(1)*(1:numel(eeeea))/numel(eeeea)+c(2);%fitted real(epsilon)
c=polyfit((1:numel(eeeea))'/numel(eeeea),imag(eeeea),1);
fitted_epi=c(1)*(1:numel(eeeea))/numel(eeeea)+c(2);%fitted imag(epsilon)
c=polyfit((1:numel(uuuua))'/numel(uuuua),real(uuuua),1);
fitted_mur=c(1)*(1:numel(uuuua))/numel(uuuua)+c(2);%fitted real(mu)
c=polyfit((1:numel(uuuua))'/numel(uuuua),imag(uuuua),1);
fitted_mui=c(1)*(1:numel(uuuua))/numel(uuuua)+c(2);%fitted imag(mu)

% Reeeea=real(eeeea);
% save Extracted_Polymer_Real_Epsilon.txt Reeeea /ascii


