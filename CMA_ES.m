% Describe variables
% Provide history too?
% Attention au nombre d'appels - compteur global


% parfor
function [X_parents_history,Y_parents_history,archive]=CMA_ES(Cf,Of,dimension,range,N,...
    scale,numMaxEvaluations,mu,lambda,tolerance,ConvCriterion,ConstraintStrat,minOrMax,archive_size,D1,D2)
    
    if minOrMax=="min"
        Of = @(x) -Of(x);
    end    
   
    
    count=0;
    generation=0;
    numMaxGenerations=ceil(numMaxEvaluations/lambda+1);
    X_parents_history=cell(1,numMaxGenerations);
    Y_parents_history=cell(1,numMaxGenerations);
    archive={zeros(dimension,0),zeros(1,0)};
    Theta=GenerateParametersCMA(dimension,scale);
    
    M=GeneratePopulation(N,range,dimension,Cf);
    [MY,count]=evaluate(M,Of,count);
    [MYsorted,I]=sort(MY,'descend');
    m=M(:,I(1));

    while count<=numMaxEvaluations-lambda | generation==0    % to be modified
        generation=generation+1;
        [X,Y,Theta,m,count,died_out]=MutateCMA(m,Theta,lambda,mu,ConstraintStrat,Of,Cf,generation,scale,count); 
        if died_out
            m=GeneratePopulation(1,range,dimension,Cf);
            Theta=GenerateParametersCMA(dimension,scale);
        end
        Y_parents_history{1,generation}=Y;
        X_parents_history{1,generation}=X;
        archive=update_archive(archive,X,Y,archive_size,D1,D2);
    end
    %to be modified
    count
    if minOrMax=="min"
        for i=1:size(Y_parents_history,2)
            Y_parents_history{1,i}=-Y_parents_history{1,i};
        end
        archive{2}=-archive{2};
    end
    Y_parents_history=Y_parents_history(1,1:generation);
    X_parents_history=X_parents_history(1,1:generation);
end


function par=CMA_par(dimension,w)
    cc=4/(dimension+4);
    mueff=1/sum(w.*w);
    csigma=(mueff+2)/(dimension+mueff+3);
    dsigma=1+2*max(0,sqrt((mueff-1)/(dimension+1))-1)+csigma;
    mucov=mueff;
    ccov=2/(mucov*(dimension+sqrt(2))^2)+(1-1/mucov)*min(1,(2*mueff-1)/((dimension+2)^2+mueff));
    Econstant=sqrt(2)*gamma(dimension/2+0.5)/gamma(dimension/2);
    par={cc,mueff,csigma,dsigma,mucov,ccov,Econstant};
end


function X=GeneratePopulation(N,range,dimension,Cf)
    X=zeros(dimension,N);
    n=0;
    while n<N
        x=2*range*rand(dimension,1)-range;
        if feasible(x,Cf)    % Redundant with the previous line in our special case
            X(:,n+1)=x;
            n=n+1;
        end
    end
end

function Theta=GenerateParametersCMA(dimension,scale) 
        psigma=zeros(dimension,1);
        pc=zeros(dimension,1);
        C=eye(dimension);
        sqrtC=eye(dimension);
        sigma=scale;
        Theta={psigma,pc,C,sqrtC,sigma};
end


% Returns logicals
function [selectedIndices,failure]=ESSelect(Y,mu)
    n=size(Y,2);
    if n==0
        failure=true;
    else
        failure=false;
    end
    [B,I]=sort(Y,'descend');
    selectedIndices=I(1:min(mu,n));
end


function [mkp1]=newMeanCMA(X_parents)
mu=size(X_parents,2);
w=((log(mu+1)-log(1:mu))/sum(log(mu+1)-log(1:mu)))';  %option 2
mkp1=X_parents*w;

end


        
function [X_mutated,Y_mutated,Theta,new_mean,count,died_out]=MutateCMA(old_mean,Theta,lambda,mu,ConstraintStrat,Of,Cf,generation,scale,count)
        died_out=false;
        d=size(old_mean,1);
        new_X=zeros(d,lambda);
        feasible_indices=true(1,lambda);
        parfor i=1:lambda 
            redraw=true;
            while redraw==true
                new_X(:,i)=old_mean+Theta{5}*Theta{4}*normrnd(0,1,d,1);
                redraw = not(feasible(new_X(:,i),Cf));
                if ConstraintStrat=="accept"
                    redraw=false;
                    feasible_indices(1,i)=feasible(new_X(:,i),Cf);
                end
            end
        end
        new_X=new_X(:,feasible_indices);
        [new_Y,count] = evaluate(new_X,Of,count);
        [selectedIndices,failure]=ESSelect(new_Y,mu);
        if failure
            died_out=true;
            "Population had died out"
            new_mean=old_mean;
            X_mutated=zeros(d,0);
            Y_mutated=zeros(d,0);
            return
        end
        X_mutated=new_X(:,selectedIndices);
        Y_mutated=new_Y(:,selectedIndices);
        [new_mean]=newMeanCMA(X_mutated);
        
        real_mu=size(X_mutated,2);
        w=((log(real_mu+1)-log(1:real_mu))/sum(log(real_mu+1)-log(1:real_mu)))';  %option 2

        par=CMA_par(d,w);
        cc=par{1};
        mueff=par{2};
        csigma=par{3};
        dsigma=par{4};
        mucov=par{5};
        ccov=par{6};
        Econstant=par{7};
        
      
        old_psigma=Theta{1};
        old_pc=Theta{2};
        old_C=Theta{3};
        old_sqrtC=Theta{4};
        old_sigma=Theta{5};
        psigma=(1-csigma)*old_psigma+sqrt(csigma*(2-csigma))*...
                    sqrt(mueff)/old_sigma*(old_sqrtC\(new_mean-old_mean));
        if norm(psigma)/sqrt(1-(1-csigma)^(2*generation+2)) < (1.5+1/(d-0.5))*Econstant
               H=1;
        else
               H=0;
        end
        pc=(1-cc)*old_pc+H*sqrt(cc*(2-cc))*sqrt(mueff)/old_sigma*(new_mean-old_mean);
        C=(1-ccov)*old_C+ccov/mucov*(pc*pc')+...
        ccov*(1-1/mucov)/old_sigma^2*(X_mutated-repmat(old_mean,1,size(X_mutated,2)))*diag(w)*(X_mutated-repmat(old_mean,1,size(X_mutated,2)))';                 
        sqrtC=sqrtm(C);
        sigma=min(old_sigma*exp(csigma/dsigma*(norm(psigma)/Econstant -1)),1000);
        Theta={psigma,pc,C,sqrtC,sigma};
end

% Evaluates the function Of on the columns (each corresponds to a solution) of the
% matrix X
% Updates count, the number of evaluations
function [Y,count] = evaluate(X,Of,count)
    n=size(X,2);
    Y=zeros(1,n);
    parfor i=1:size(X,2)
        Y(i)=Of(X(:,i)');
    end
    count = count+n;
end



function feas=feasible(x,Cf)
    feas=true;
    for i=1:size(Cf)
        if Cf{1,i}(x)<0
            feas= false;
        end
    end
end

%to be modified: implémenter un truc plus subtil (avec évolution par
%rapport aux précédentes valeurs
function hasConverged=TestConvergence(Y,tolerance,ConvCriterion)
    hasConverged=false;
    if ConvCriterion=="AbsDiff"
        if abs(max(Y)-min(Y))<tolerance
            hasConverged=true;
            "ES has converged"
        end
    elseif ConvCriterion=="RelDiff"
        if abs(max(Y)-min(Y))<tolerance*abs(mean(Y))
            hasConverged=true;
            "ES has converged"
        end
    else
        "Choose valid convergence criterion"
        return
    end
end
