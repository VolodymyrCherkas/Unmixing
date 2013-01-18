%%
for i=1:length(BG)
    Sp(:,i)=S{i}(:,1);
    SpOME{i}.wlconf=OME{i}.wlconf;
    SpOME{i}.expconf=OME{i}.expconf;
    Spn(:,i)=S{i}(:,1)./OME{i}.expconf(:); % normalized by exposure;
end

%%
S_af.mean=mean(Spn')';
S_af.var=var(Spn')';
S_af.OME=SpOME;
S_af.normalized=1;

S_TFP.mean=mean(Spn')';
S_TFP.var=var(Spn')';
S_TFP.OME=SpOME;
S_TFP.normalized=1;
S_TFP.percell=Spn;
% 4 - Glia
% 5 - dead
% 8 - dead
%
S_YFP.mean=mean(Spn')';
S_YFP.var=var(Spn')';
S_YFP.OME=SpOME;
S_YFP.normalized=1;
S_YFP.percell=Spn;

S_TAN.mean=mean(Spn')';
S_TAN.var=var(Spn')';
S_TAN.OME=SpOME;
S_TAN.normalized=1;
S_TAN.percell=Spn;

S_TFP_noaf.mean=S_TFP.mean-S_af.mean/S_af.mean(6)*S_TFP.mean(6);
S_TFP_noaf.var=S_TFP.var;
S_TFP_noaf.OME=S_TFP.OME;
S_TFP_noaf.normalized=1;
for i=1:size(S_TFP.percell,2)
    S_TFP_noaf.percell(:,i)=S_TFP.percell(:,i)-S_af.mean/S_af.mean(6)*S_TFP.percell(6,i);
end
%%
S_YFP_noaf.mean=S_YFP.mean-S_af.mean/S_af.mean(6)*S_YFP.mean(6);
S_YFP_noaf.var=S_YFP.var;
S_YFP_noaf.OME=S_YFP.OME;
S_YFP_noaf.normalized=1;
for i=1:size(S_YFP.percell,2)
    S_YFP_noaf.percell(:,i)=S_YFP.percell(:,i)-S_af.mean/S_af.mean(6)*S_YFP.percell(6,i);
end
%%
% S_TFP_noaf; S_YFP_noaf; S_af; S_FRET - ?
% S_FRET(1:2)~=0; S_FRET(3:6)==0;
S_FRET_minus=-S_TFP_noaf.mean;
S_FRET_minus(3:6)=0;
S_FRET_plus=S_YFP_noaf.mean;
S_FRET_plus(3:6)=0;
%%
i1=ones(1,25)
T=nonnegative_unmix([S_YFP_noaf.mean S_af.mean S_FRET_minus S_FRET_plus], S{1}(:,1)'-S{3}(:,25)'/4);
T=nonnegative_unmix([S_TFP_noaf.mean S_YFP_noaf.mean S_af.mean], S{1}(:,1)'-S{1}(:,25)');

% C_TFP0 + C_YFP0 + C_af0 + FRET0 = S_TAN0
% C_TFPf + C_YFPf + C_aff + FRETf = S_TANf
% FRETf / FRET0 = C_TFPf * C_YFPf / C_TFP0 / C_YFP0


% FRET = FRET_plus + FRET_minus = S_YFP_noaf.mean(1:2) -
% - G*S_TFP_noaf.mean(1:2)? (or 1:4)
% G = const
% (1-G) * C_TFP + C_YFP + C_af + S_YFP_noaf.mean(1:2) = S_TAN
% FRET = k * C_TFP * C_YFP
% k = const

for i=1:size(S_YFP_noaf.percell,2)
    S_YFP_noaf.norm2(:,i)=S_YFP_noaf.percell(:,i)/max(S_YFP_noaf.percell(:,i))*max(S_YFP_noaf.mean);
end
S_YFP_noaf.mean2=mean(S_YFP_noaf.norm2,2)

