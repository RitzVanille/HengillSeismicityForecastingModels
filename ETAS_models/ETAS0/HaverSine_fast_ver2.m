function [d]=HaverSine_fast_ver2(X,Y)
%tic
X=X';
X=X*pi/180;
Y=Y*pi/180;
coslat1=cos(X(2,:));
coslat2=cos(Y(:,2));

T1=sin(bsxfun(@minus, X(2,:), Y(:,2))/2).^2;
T2=bsxfun(@times, coslat1, coslat2).*sin(bsxfun(@minus, X(1,:), Y(:,1))/2).^2;
a=T1+T2; 
a(a>1)=1;

c=2*atan2(sqrt(a),sqrt(1-a));
d=6371*c;
%toc

end