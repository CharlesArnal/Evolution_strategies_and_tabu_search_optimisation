% Describe variables
% Provide history too?
% Attention au nombre d'appels - compteur global



function [X_parents_history,Y_parents_history,archive,graphdata]=ES(Cf,Of,my_gradient,dimension,range,...
    scale,N,numMaxEvaluation,m1,m2,m3,mu,lambda,tolerance,IDCoeff,StabilizationCoeff,perturbation,...
    veteransNumber,CombinationContVarOp,CombinationStratParOp,CombPairGlob,SelectionStrat,...
    ConstraintStrat,ConvCriterion,endogamy,endogamy_scale,islandsNumber,migrationRate,migrationFrequency,minOrMax,archive_size,D1,D2)
    if minOrMax=="min"
        Of = @(x) -Of(x);
    end    
    
    graphdata=zeros(2,0);
    
    count=0;
    generation=0;
    X_parents_history=cell(1,0);
    Y_parents_history=cell(1,0);
    archive={zeros(dimension,0),zeros(1,0)};
    X=GeneratePopulation(N*islandsNumber,range,dimension,Cf);
    [Y,count] = evaluate(X,Of,count);
    Theta=GenerateStrategyParameters(N*islandsNumber,dimension,scale);
    if IDCoeff ==0
        margin=lambda*islandsNumber;
    else
        margin=lambda*islandsNumber*(1+dimension);
    end
    while count<=numMaxEvaluation-margin &  TestConvergence(Y,tolerance,ConvCriterion) == false      % Corriger
        [selectedIndices,failure]=ESSelect(Y,mu,SelectionStrat,islandsNumber);
        if failure
            "Population had died out"
            return
        end
        X_parents_history{1,generation+1}=X(:,selectedIndices);
        Y_parents_history{1,generation+1}=Y(:,selectedIndices);
        archive=update_archive(archive,X(:,selectedIndices),Y(:,selectedIndices),archive_size,D1,D2);
        if mod(generation,migrationFrequency)==0
            [X_migrated,Theta_migrated]=Migrate(X(:,selectedIndices),Theta(1,selectedIndices),migrationRate);
            [X_children,Theta_children]=Reproduce(X_migrated,Theta_migrated,lambda,CombinationContVarOp,CombinationStratParOp,CombPairGlob,endogamy,endogamy_scale,islandsNumber);
        else
            [X_children,Theta_children]=Reproduce(X(:,selectedIndices),Theta(1,selectedIndices),lambda,CombinationContVarOp,CombinationStratParOp,CombPairGlob,endogamy,endogamy_scale,islandsNumber);
        end
        [X_children,Theta_children,count]=Mutate(X_children,Theta_children,count,m1,m2,m3,IDCoeff,perturbation,my_gradient,ConstraintStrat,Cf,StabilizationCoeff); 
        [Y_children,count] = evaluate(X_children,Of,count);
        if count<=numMaxEvaluation   % to be commented
            [veteransIndices,failure]=ESSelect(Y,veteransNumber,SelectionStrat,islandsNumber);
            [X,Y,Theta]=newGeneration(X,Y,Theta,X_children,Y_children,Theta_children,veteransIndices,veteransNumber,islandsNumber);
            generation=generation+1;
        end
        
        graphdata=[graphdata,[count;-archive{2}(1)]];
    end
    %to be modified
    count
    if minOrMax=="min"
        for i=1:size(Y_parents_history,2)
            Y_parents_history{1,i}=-Y_parents_history{1,i};
        end
        archive{2}=-archive{2};
    end
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

function Theta=GenerateStrategyParameters(N,dimension,scale) 
    Theta=cell(1,N);
    Theta(:)={eye(dimension)*scale};     
end


% Returns logicals
function [selectedIndices,failure]=ESSelect(Y,mu,SelectionStrat,islandsNumber)
    n=size(Y,2);
    k=n/islandsNumber;
    if n==0
        failure=true;
    else
        failure=false;
    end
    selectedIndices=logical(zeros(1,n));
    for j=1:islandsNumber
        if SelectionStrat=="mufirst"
            [B,I]=sort(Y(1+(j-1)*k:j*k),'descend');
            %selectedIndices=I(1:min(mu,n));     % better, automatically ordered
            selectedIndices((j-1)*k+I(1:min(mu,k)))=true;
        else
            "Choose selection strategy"
            return
        end
    end
end

function [X_migrated,Theta_migrated]=Migrate(X,Theta,migrationRate)
    X_migrated=X;
    Theta_migrated=Theta;
    n=size(X,2);
    k=floor(migrationRate*n);    % the number of migrators
    my_perm=randperm(n);
    migrators_indices=my_perm(1:k);     % randomly pick the indices of the solutions that will migrate
    X_migrators=X(migrators_indices);   % the migrators
    Theta_migrators=Theta(migrators_indices);   % the migrators's parameters
    my_perm_2=randperm(k);
    % randomly migrating the migrators
    X_migrated(migrators_indices)=X_migrators(my_perm_2);
    Theta_migrated(migrators_indices)=Theta_migrators(my_perm_2);
end
        

% Careful: potentially less than mu selectedIndices
function [X_children,Theta_children]=Reproduce(X_parents,Theta_parents,lambda,CombinationContVarOp,CombinationStratParOp,CombPairGlob,endogamy,endogamy_scale,islandsNumber)
   k=size(X_parents,2);
   l=k/islandsNumber;
   X_children=zeros(size(X_parents,1),0);
   Theta_children=cell(1,0);
       for t=1:islandsNumber      % reproduction only inside islands (migration occurs before)
           temp_X=zeros(size(X_parents,1),lambda);      % necessary to avoid indexing issues within parfor loop
           temp_Theta=cell(1,lambda);
           parfor i=1:lambda
                if k==1
                    X_children(:,i)=X_parents(:,1);
                else
                        if CombPairGlob=="pair"
                            [a,b]=SelectParents(X_parents(:,1+(t-1)*l:t*l),endogamy,endogamy_scale);
                            temp_X(:,i)= CombControlVariables([X_parents(:,a+(t-1)*l),X_parents(:,b+(t-1)*l)],CombinationContVarOp);
                            temp_Theta{1,i}=CombStratPar({Theta_parents{1,a+(t-1)*l},Theta_parents{1,b+(t-1)*l}},CombinationStratParOp);
                        elseif CombPairGlob=="global"
                            temp_X(:,i)= CombControlVariables(X_parents(:,1+(t-1)*l:t*l),CombinationContVarOp);
                            temp_Theta{1,i}=CombStratPar(Theta_parents(1,1+(t-1)*l:t*l),CombinationStratParOp);
                        end
                end
           end
           X_children=[X_children,temp_X];
           Theta_children=[Theta_children,temp_Theta];
       end
end

function [a,b]=SelectParents(X_parents,endogamy,endogamy_scale)
    k=size(X_parents,2);
    my_perm=randperm(k);
    a=my_perm(1);
    if endogamy==false       
        b=my_perm(2);
    else
        probas=zeros(1,k);
        for t=1:k
            if t~=a
                probas(1,t)=min(0.5,1/(1+ min(9, norm(X_parents(:,a)-X_parents(:,t))/endogamy_scale)));
            end
        end
        b=genBern(probas);
    end

end


function [X_mutated,Theta_mutated,count]=Mutate(X,Theta,count,m1,m2,m3,IDCoeff,perturbation,my_gradient,ConstraintStrat,Cf,StabilizationCoeff)
    n=size(X,2);
    d=size(X,1);
    feasible_indices=true(1,n);
    X_mutated=zeros(d,n);
    Theta_mutated=cell(1,n);
    parfor i=1:n
            redraw=true;
            redrawing_count=0;
            while redraw==true
                if IDCoeff==0
                    if perturbation=="normal"
                        X_mutated(:,i)=X(:,i)+Theta{1,i}*normrnd(0,1,d,1);
                    elseif perturbation=="student"
                        X_mutated(:,i)=X(:,i)+Theta{1,i}*trnd(1,d,1);
                    end
                else
                    if perturbation=="normal"
                        X_mutated(:,i)=X(:,i)+Theta{1,i}*normrnd(0,1,d,1)+IDCoeff*my_gradient(X(:,i)')'*exprnd(1);
                    elseif perturbation=="student"
                        X_mutated(:,i)=X(:,i)+Theta{1,i}*trnd(1,d,1)+IDCoeff*my_gradient(X(:,i)')'*exprnd(1);
                    end  
                    count=count+d;
                end
                redraw = not(feasible(X_mutated(:,i),Cf));
                redrawing_count=redrawing_count+1;
                if ConstraintStrat=="accept"
                    redraw=false;
                    feasible_indices(1,i)=feasible(X_mutated(:,i),Cf);
                end
                if redrawing_count >2^d*10
                    Theta{1,i}=Theta{1,i}*eye(d)*0.1;
                    redrawing_count=0;
                    "covariance rescaled"
                end
            end
            
            
            Theta_mutated{1,i}=rotationPerturbation(d,m3)*Theta{1,i}*variancePerturbation(d,m1,m2,StabilizationCoeff);
        
    end

    X_mutated=X_mutated(:,feasible_indices);

end

function R=rotationPerturbation(d,m3)
    R=eye(d);
    for i=1:d
        for j=i+1:d
            xji=normrnd(0,1)*m3;
            R(i,:)=cos(xji)*R(i,:)-sin(xji)*R(j,:);
            R(j,:)=sin(xji)*R(i,:)+cos(xji)*R(j,:);
        end
    end
end

function D=variancePerturbation(d,m1,m2,StabilizationCoeff)
    D=zeros(d);
    x0=normrnd(0,1);
    for i=1:d
        xl=normrnd(0,1);
        D(i,i)=exp(m1*x0+m2*xl+StabilizationCoeff);
    end
end


function [X,Y,Theta]=newGeneration(X,Y,Theta,X_children,Y_children,Theta_children,veteransIndices,veteransNumber,islandsNumber)
       X_veterans=X(:,veteransIndices);
       Y_veterans=Y(1,veteransIndices);
       lambda=size(X_children,2)/islandsNumber;
       Theta_veterans=Theta(1,veteransIndices);
       X=zeros(size(X_children,1),islandsNumber*(lambda+veteransNumber));
       Y=zeros(1,islandsNumber*(lambda+veteransNumber));
       Theta=cell(1,islandsNumber*(lambda+veteransNumber));
       for i=1:islandsNumber
            X(:,1+(i-1)*(lambda+veteransNumber):i*(lambda+veteransNumber))=[X_veterans(:,1+(i-1)*veteransNumber:i*veteransNumber),X_children(:,1+(i-1)*lambda:i*lambda)];
            Y(1,1+(i-1)*(lambda+veteransNumber):i*(lambda+veteransNumber))=[Y_veterans(1,1+(i-1)*veteransNumber:i*veteransNumber),Y_children(1,1+(i-1)*lambda:i*lambda)];
            Theta(1,1+(i-1)*(lambda+veteransNumber):i*(lambda+veteransNumber))=[Theta_veterans(1,1+(i-1)*veteransNumber:i*veteransNumber),Theta_children(1,1+(i-1)*lambda:i*lambda)];
       end
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

function x= CombControlVariables(X,CombOp)
    if CombOp=="Intermediate"
        x=mean(X,2);
    elseif CombOp=="Discrete"
        [d,k]=size(X);
        x=zeros(d,1);
        for i=1:d
            j=genBern(ones(1,k));
            x(i)=X(i,j);
        end
    end
    
end


function Theta=CombStratPar(THETA,CombOp)
d=size(THETA{1,1},1);
Theta=zeros(d);
k=size(THETA,2);
if CombOp== "Discrete"
    % does not make sense
    for i=1:d
        for j=i:d
            t=genBern(ones(1,k));
            Theta(i,j)=THETA{1,t}(i,j);
            Theta(j,i)=Theta(i,j);
        end 
    end
elseif CombOp=="Intermediate"
    Theta=zeros(d);
    k=size(THETA,2);
    for t=1:k
        Theta=Theta+THETA{1,t}*THETA{1,t}'/k;
    end
    Theta=chol(Theta)';
    %{
    for i=1:size(Theta1,1)
        for j=1:size(Theta1,2)
            if i==j
                Theta(i,i)=(Theta1(i,i) + Theta2(i,i))/2;
            else
                Theta(i,j)=angleAverage(Theta1(i,j),Theta2(i,j));
            end
        end 
    end
    %}
end
end


  
function a=angleAverage(a1,a2)
    r=mod(a2-a1,2*pi);
    if r<=pi
        a=a1+r/2;
    else
       a=a1-(2*pi - r)/2; 
    end
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

function A=TransformationMatrix(T)
    d=size(T,1);
    R=eye(d);
    D=zeros(d);
    for i=1:d
        D(i,i)=T(i,i);
        for j=i+1:d
            R(i,:)=cos(T(i,j))*R(i,:)-sin(T(i,j))*R(j,:);
            R(j,:)=sin(T(i,j))*R(i,:)+cos(T(i,j))*R(j,:);
        end
    end
    A=R*D;
end

% function to sample from an unnormalized generalized Bernoulli distribution

function j = genBern(f)
t = rand()*sum(f);
a = f(1);
j = 1;
while a < t
  j = j+1;
  a = a+f(j);
end
end
