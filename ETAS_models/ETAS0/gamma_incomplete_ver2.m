function fun=gamma_incomplete_ver2(x,a)
% GAMMA_INCOMPLETE evaluates the upper incomplete gamma function
% (incomplete gamma function of the second kind) $\Gamma(a,x)$
% at non-negative values of the argument. This function extends the
% MATLAB function gammainc to negative values of the parameter a.
%
%  SYNOPSIS:  fun=gamma_incomplete(x,a)
%
%
%  INPUT  x       : function argument
%         a       : parameter
%
%  OUTPUT fun     : a vector of the same length as x; it contains NaN values
%  at places where elements of x are negative.
%
%  REMARK This function extends the MATLAB function gammainc
%  to negative values of the parameter a.
%
%  EXAMPLES
%
%     x=0.01:0.01:8;
%     f=gamma_incomplete(x,1);
%     plot(x,f);
%
%
%     x=0.001:0.001:.1;
%     f=gamma_incomplete(x,-2.3);
%     plot(x,f);
%
%
%     x=0.001:0.001:.2;
%     f=gamma_incomplete(x,-1);
%     plot(x,f);
if size(a,1)==1
    a=a*ones(size(x));
end

if size(x,1)==1
    x=x*ones(size(a));
end

p=(x>=0);
q=(x<0);

r=(a==0);
s=(a>0);
n1=ceil(abs(a));
t=(a+n1)==0 & a<0;
u=(a+n1)>0 & a<0;


fpr = gamma_inc1(x(r & p));
fps = gamma_inc2(a(s & p),x(s & p));
fpt = gamma_inc3(a(t & p),x(t & p));
fpu = gamma_inc4(a(u & p),x(u & p));

fun=zeros(length(x),1);
fun(r & p)=fpr;
fun(s & p)=fps;
fun(t & p)=fpt;
fun(u & p)=fpu;
fun(q)=NaN;
%%
    function res=gamma_inc1(x)
        
        res= expint(x);
        
        
    end
%%
    function res=gamma_inc2(a,x)
        
        res=gamma(a).*gammainc(x,a,'upper');
    end
%%
    function res=gamma_inc3(a,x)
        
        n=-a;
        maxn=max(n);
        pow1=[1:maxn];
        Sold=zeros(size(x));
        for i=1:maxn
            %pow1=[1:n(i)];
            S=((-1).^(pow1(i)-1).*gamma(pow1(i)).*x.^-pow1(i));
            S(n<i)=0;
            Sold=Sold+S;
        end
        
        deno=((-1).^n.*gamma(n+1));
        res=(expint(x)-exp(-x).*Sold)./deno;
        
        
    end

%%
    function res=gamma_inc4(a,x)
        
        n=ceil(abs(a));
        maxn=max(n);
        Sold=zeros(size(x));
        fix=gamma(a);
        prod=ones(length(a),1);
        for i=1:maxn
            prod=prod.*(a+i-1);
            S=x.^(i-1)./prod;
            S(i>n)=0;
            Sold=Sold+S;
        end
        
        mult=fix.*(a+n)./gamma(a+n+1);
        res=(gamma(a+n).*gammainc(x,a+n,'upper').*mult-exp(-x).*x.^a.*Sold);
        
    end
end
