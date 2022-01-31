function [missnums] = Wrapper_FailedEval_ver3(foldnam,filepattern2,Maux,Mpri,KDEBKG,OmTyp,btyp,isotyp,...
    srclen,mthres,prodtyp,nprod,seqSpec,minFamSz,tresol,SimUpdateDur)
%
%%
%{
%
foldnam = './outputvars/';

%%

Maux = 0.3
Mpri = 0.3
KDEBKG = 1
OmTyp = 1
btyp = 1
isotyp = 1
srclen = 50
mthres = 0
prodtyp = 0
nprod = 0
tresol = 30
SimUpdateDur = 30
seqSpec = 0;
minFamSz = 0;



filepattern2 = 'Eval_sttm%f_Maux%f_Mpri%f_KDEBKG%f_OmTyp%f_btyp%f_isotyp%f_srclen%f_mthres%f_prodtyp%f_seqSpec%f_minFamSz%f_nprod%f_tresol%f.mat';

[missnums3] = Wrapper_FailedEval_ver3(foldnam,filepattern2,Maux,Mpri,KDEBKG,OmTyp,btyp,isotyp,srclen,mthres,prodtyp,nprod,seqSpec,minFamSz,tresol,SimUpdateDur)


%}
%
%% decide which parameter to use

catnam   =  'Coseismiq';
Regionnam = 'Iceland-Hengill-Active_pval0.1';
Prefix = strcat(catnam,'_',Regionnam);

begsttm = datenum(2020,2,1);

cat0   =  importdata(strcat(Prefix,'.mat'));
modpath = strcat(foldnam,Prefix,'/Eval_updateDur',num2str(SimUpdateDur),'/');

%
endsttm_inv = max(cat0.pricat.datenum);
sttmlist = [begsttm:SimUpdateDur:endsttm_inv];

[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14] = ndgrid(sttmlist,Maux,Mpri,KDEBKG,OmTyp...
    ,btyp,isotyp,srclen,mthres,prodtyp,seqSpec,minFamSz,nprod,tresol);
Varlist = [x1(:),x2(:),x3(:),x4(:),x5(:),x6(:),x7(:),x8(:),x9(:),x10(:),x11(:),x12(:),x13(:),x14(:)];
missnums =  FindFailedEvals(modpath,Varlist,filepattern2);

missnums = array2table(missnums,'VariableNames',{'sttm','Maux','Mpri','KDE_BKG',...
    'OmTyp','btyp','isotyp','srclen','mthres','prodtyp','seqSpec','minFamSz','nprod','tresol'});
end
