% Computes Rana's function, expects a line vector x=[x1,...,xn]
% Simple application of the formula
function f=RanaFun(x)  
    n=size(x,2);    
    f=0;
    x_i=x(1:n-1);   % x_i in the formula
    x_ip1=x(2:n);   % x_i+1 in the formula
    v1=sqrt(abs(x_ip1+x_i+1));
    v2=sqrt(abs(x_ip1-x_i+1));
    % vectorized for elegance and clarity
    f=sum(x_i .* cos(v1) .* sin(v2) + (1+x_ip1) .* cos(v2) .* sin(v1));
end
