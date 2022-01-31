function [nmag]=BinMags(mag,bin0,binsz)
% binning is problematic for larger bin size
binedges=[bin0-binsz/2:binsz:11+binsz/2];
[N,edges,bin] = histcounts(mag,binedges);

nmag=bin0+(bin-1)*binsz;
nmag=round(nmag,1);
end