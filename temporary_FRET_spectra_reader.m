%% temporary part
S_tan1.spectrum=Spectra{1}';
S_tan1.exp=[3 3 3 3 3 3];
S_tan1.wl=[435 435 505 505 575 575];
S_tan1.emch=[0 1 0 1 0 1];

S_tan2.spectrum=Spectra{2}';
S_tan2.exp=[3 3 3 3 3 3];
S_tan2.wl=[435 435 505 505 575 575];
S_tan2.emch=[0 1 0 1 0 1];

load('Ref_TFP')
load('Ref_YFP')

i=[1 3 5 2 4 6];
figure; plot(1:6, Sd.spectrum(i)./Sd.exp(i),'b', 1:6, Sa.spectrum(i)./Sa.exp(i), 'r', 1:6, S_tan1.spectrum(i)./S_tan1.exp(i),'g',1:6, S_tan2.spectrum(i)./S_tan2.exp(i),'k');
%% calculate FRET+/- spectrum

i0=find(~S_tan1.emch);
i1=find(S_tan1.emch);
i0=[1 3 5];
i1=[2 4 6];

St1_0=(S_tan1.spectrum(i0))';%./S_tan1.exp(i0))';
St1_1=(S_tan1.spectrum(i1))';%./S_tan1.exp(i1))';
St2_0=(S_tan2.spectrum(i0))';%./S_tan2.exp(i0))';
St2_1=(S_tan2.spectrum(i1))';%./S_tan2.exp(i1))';

Sd_0=(Sd.spectrum(i0))';%./Sd.exp(i0))';
Sd_1=(Sd.spectrum(i1))';%./Sd.exp(i1))';
Sa_0=(Sa.spectrum(i0))';%./Sa.exp(i0))';
Sa_1=(Sa.spectrum(i1))';%./Sa.exp(i1))';

Saf_0=(S_a.S.spectrum(i0)./S_a.S.exp(i0))';
Saf_1=(S_a.S.spectrum(i1)./S_a.S.exp(i1))';

% FRET_1==0
% FRET_0 is negative-positive
% St1_0=d1*Sd_0+a1*Sa_0+af1*Saf_0+f1*FRET;
% St1_1=d1*Sd_1+a1*Sa_1+af1*Saf_1;
% St2_0=d2*Sd_0+a2*Sa_0+af2*Saf_0+f2*FRET;
% St2_1=d2*Sd_1+a2*Sa_1+af2*Saf_1;
% assumptions:
% a1/a2=f1/f2; - good for tandem - all donors are good, some acceptors are
% bleached, no free acceptors available
% d1=d2? % first will try not using it
r1=nonnegative_unmix([Sd_1 Sa_1 Saf_1], St1_1');
a1=r1(2); d1=r1(1); af1=r1(3);
r2=nonnegative_unmix([Sd_1 Sa_1 Saf_1], St2_1');
a2=r2(2); d2=r2(1); af2=r2(3);
% f1*FRET=St1_0-d1*Sd_0-a1*Sa_0;
% f2*FRET=St2_0-d2*Sd_0-a2*Sa_0;
F1=St1_0-d1*Sd_0-a1*Sa_0-af1*Saf_0;
F2=St2_0-d2*Sd_0-a2*Sa_0-af2*Saf_0;
%%
Sd=[Sd_0' Sd_1']';
Sa=[Sa_0' Sa_1']';
St1=[St1_0' St1_1']';
St2=[St2_0' St2_1']';
%% load different spectra
Rd1=load('Ref_TFP1');
Rd2=load('Ref_TFP2');
Rd3=load('Ref_TFP3');
Rd4=load('Ref_TFP');
Ra1=load('Ref_YFP1');
Ra2=load('Ref_YFP2');
Ra3=load('Ref_YFP');
%% combine them, check and average
Rd{1}=Rd1.S;
Rd{2}=Rd2.S;
Rd{3}=Rd3.S;
Rd{4}=Rd4.Sd;
Ra{1}=Ra1.S;
Ra{2}=Ra2.S;
Ra{3}=Ra3.Sa;
i0=[1 3 7];
i1=[2 4 8];
figure; plot(i0, Rd{1}.spectrum(i0),'r',i0, Rd{2}.spectrum(i0)*6,'g',i0,Rd{3}.spectrum(i0),'b',[1 3 5],Rd{4}.spectrum([1 3 5]),'k');
figure; plot(i1, Rd{1}.spectrum(i1),'r',i1, Rd{2}.spectrum(i1)*6,'g',i1,Rd{3}.spectrum(i1),'b',[2 4 6],Rd{4}.spectrum([2 4 6]),'k');
i=[1 3 7 2 4 8];
Rmd=mean(([Rd{1}.spectrum(i)', Rd{2}.spectrum(i)',Rd{3}.spectrum(i)', Rd{4}.spectrum([1 3 5 2 4 6])'])');
i=[1 3 5 2 4 6];
Rma=mean(([Ra{1}.spectrum(i)', Ra{2}.spectrum(i)',Ra{3}.spectrum(i)'])');
%%
i0=[1 3 5];
i1=[2 4 6];

Sd_0=Rmd(1:3)';
Sd_1=Rmd(4:6)';
Sa_0=Rma(1:3)';
Sa_1=Rma(4:6)';
St1_0=(S_tan1.spectrum(i0))';%./S_tan1.exp(i0))';
St1_1=(S_tan1.spectrum(i1))';%./S_tan1.exp(i1))';
St2_0=(S_tan2.spectrum(i0))';%./S_tan2.exp(i0))';
St2_1=(S_tan2.spectrum(i1))';%./S_tan2.exp(i1))';

% Saf_1(3)=[];
% Saf_0(3)=[];

r1=nonnegative_unmix([Sd_1 Sa_1 Saf_1], St1_1');
a1=r1(2); d1=r1(1); af1=r1(3);
r2=nonnegative_unmix([Sd_1 Sa_1 Saf_1], St2_1');
a2=r2(2); d2=r2(1); af2=r2(3);
% f1*FRET=St1_0-d1*Sd_0-a1*Sa_0;
% f2*FRET=St2_0-d2*Sd_0-a2*Sa_0;
F1raw=St1_0-d1*Sd_0-a1*Sa_0-af1*Saf_0;
F2raw=St2_0-d2*Sd_0-a2*Sa_0-af2*Saf_0;
%save('tandem_test_04.10.2012.mat')
%%
r1=nonnegative_unmix([Sd_1 Sa_1 Saf_1], S{1}(i1,:)');
a1=r1(:,2); d1=r1(:,1); af1=r1(:,3);
for j=1:length(a1)
    F(j,:)=St1_0'-d1(j)*Sd_0'-a1(j)*Sa_0'-af1(j)*Saf_0';
end
for k=1:3
    Fm(k)=mean(F(11:20,k)./a1(11:20)./d1(11:20)*a1(1)*d1(1),1)/3;  %normalized by exposure
end
Sf.spectrum=Fm;
Sf.spectrum(4:6)=zeros;
Sf.exp=[3 3 3 3 3 3];
Sf.wl=[435 505 575 435 505 575];
Sf.emch=[0 0 0 1 1 1];
save('Ref_FRET.mat','Sf');
