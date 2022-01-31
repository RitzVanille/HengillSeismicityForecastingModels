function [missnums] =  FindFailedEvals(modpath,Varlist,filepattern)
% modfoldnam =  'D:/datadump/KDE_BKG/outputvars/Iso_OmoriX_KDEBKG_PPTest_since1990/Inv/';
% modpath = modfoldnam;
% begsttm=datenum(1990,1,1); endsttm=datenum(2020,4,24);
% sttmlist=[begsttm:30:endsttm];
% [X,Y,Z] = meshgrid(sttmlist(1:end-1),1,[0,1,2]);
% Varlist = [X(:),Y(:),Z(:)];
% filepattern = 'Pars_sttm%d_Maux5_Mpri5_KDEBKG%d_OmTyp%d.mat'
% savnam = 'FailInversionList.mat'


fllist  = dir(strcat(modpath,'Eval*.mat'));
if ~isempty(fllist)
    for i = 1:size(fllist)
        compnum(i,:) = sscanf(fllist(i).name, filepattern)';
    end
    [Lia,locb]      = ismember(Varlist,compnum,'rows');
    missnums = Varlist(~Lia,:);
else
    compnum = [];
    missnums = Varlist;
end

%% then find those that are completely missing

check = 1;
end