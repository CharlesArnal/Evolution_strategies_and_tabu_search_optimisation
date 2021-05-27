% Computes Rana's function's gradient, expects a line vector x=[x1,...,xn]
% Tedious but straightforward
% In terms of complexity, we consider that one evaluation of RanaFunGradient
% roughly corresponds to n evaluations of Rana's function, where n is the
% ambient dimension
function D=RanaFunGradient(x)  
    n=size(x,2);
    D1=zeros(1,n);     % derivative of the ith term of the sum wrt to x_i
    D2=zeros(1,n);     % derivative of the ith term of the sum wrt to x_i+1
    parfor i=1:n-1
       s1=sign(x(i+1)+x(i)+1);
       s2=sign(-x(i+1)+x(i)-1);
       a=sqrt(abs(x(i+1)+x(i)+1));
       b=sqrt(abs(-x(i+1)+x(i)-1));
       D1(i)= cos(a)*sin(b)...
           +0.5*x(i)*((-1)*sin(a)*sin(b)*s1/a + cos(a)*cos(b)*s2/b) ...
           +0.5*(1+x(i+1))*((-1)*sin(b)*sin(a)*s2/b + cos(b)*cos(a)*s1/a);
       D2(i+1)=0.5*x(i)*(sin(a)*sin(b)*s1/a + cos(a)*cos(b)*(-1)*s2/b)...
           +cos(a)*sin(b)...
           +0.5*(1+x(i+1))*((-1)*sin(b)*sin(a)*(-1)*s2/b + cos(b)*cos(a)*s1/a);
    end
        D=D1+D2;
end
