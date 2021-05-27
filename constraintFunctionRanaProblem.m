
% The simple function corresponding to the constraints |x|<=500
function f=constraintFunction(x)
    if sum(abs(x)>500)==0
        f=1;
    else
        f=-1;
    end
end