close all;
clear all;

% Evolution Strategy
% ------------------------------

dimension=2;
IDCoeff=0;
scale=0.1;
range=500;
numMaxEvaluation=10000;
StabilizationCoeff=log(1);
veteransNumber=0;
perturbation="normal"; % "normal"/"student"
m1=1/sqrt(2*sqrt(dimension));
m2=1/sqrt(2*dimension);
m3=0.0873;
endogamy=false;
endogamy_scale=200;
islandsNumber=1;
lambda=floor(450/islandsNumber);
mu=floor(90/islandsNumber);
N=floor(1000/islandsNumber);
migrationRate=0.2;
migrationFrequency=1;
CombinationContVarOp="Discrete";
CombinationStratParOp="Intermediate";
CombPairGlob="pair";
minOrMax="min";
ConstraintStrat="redraw";
ConvCriterion="RelDiff";
SelectionStrat="mufirst";
tolerance=0.01;
Cf={@constraintFunctionRanaProblem};
Of=@RanaFun;
my_gradient=@RanaFunGradient;
archive_size=20;
D1=20;
D2=1;

% ---
M=1;



time=zeros(1,M);
best_performance=zeros(1,M);
mean_performance=zeros(1,M);
performance_std_dev=zeros(1,M);
locateBest=zeros(dimension,100);%erase
for j=1:M
NumRuns=1;
Min_run=zeros(1,NumRuns);
tic
for i=1:NumRuns
    rng(i);
    [X_parents_history,Y_parents_history,archive,graphdata]=ES(Cf,Of,my_gradient,dimension,range,...
    scale,N,numMaxEvaluation,m1,m2,m3,mu,lambda,tolerance,IDCoeff,StabilizationCoeff,perturbation,...
    veteransNumber,CombinationContVarOp,CombinationStratParOp,CombPairGlob,SelectionStrat,...
    ConstraintStrat,ConvCriterion,endogamy,endogamy_scale,islandsNumber,migrationRate,migrationFrequency,minOrMax,archive_size,D1,D2);
    Min_run(i)=archive{2}(1);
    locateBest(:,i)=archive{1}(:,1);%erase
    i
end
time(j)=toc
best_performance(j)=min(Min_run)
mean_performance(j)=mean(Min_run)
performance_std_dev(j)=std(Min_run)
[sortedMin,I]=sort(Min_run,"ascend");%erase
best_Y=Min_run(I(1));%erase
best_X=locateBest(:,I(1))%erase
test=RanaFun(best_X')
end




min_value=zeros(1,0);
diff=zeros(1,0);
for i=1:size(X_parents_history,2)
    min_value(1,i)=min(Y_parents_history{1,i});
    diff(1,i)=max(Y_parents_history{1,i})-min(Y_parents_history{1,i});
end
f4=figure;
plot(1:size(X_parents_history,2),min_value);
hold on;
plot(1:size(X_parents_history,2),diff);


f5=figure;
if dimension==2
    u=zeros(201,201);
    for i=1:201
        for j=1:201
            u(j,i)=RanaFun([-505+i*5,-505+j*5]);    % opposé à ce qu'on penserait
        end
    end
    [X_plot,Y_plot]=meshgrid(-500:5:500,-500:5:500);
    mySurface=mesh(X_plot,Y_plot,u)
    hold on;
end

TestX=archive{1};
TestY=archive{2};
plot3(TestX(1,:),TestX(2,:),TestY,'o',"color","red",'MarkerFaceColor','red','MarkerSize',7.5)
hold on;


f6=figure;
if dimension==2
    u=zeros(201,201);
    for i=1:201
        for j=1:201
            u(j,i)=RanaFun([-505+i*5,-505+j*5]);    % opposé à ce qu'on penserait
        end
    end
    [X_plot,Y_plot]=meshgrid(-500:5:500,-500:5:500);
    mesh(X_plot,Y_plot,u)
    hold on;
end

TestX=X_parents_history{1,1};
TestY=Y_parents_history{1,1};
plot3(TestX(1,:),TestX(2,:),TestY,'o',"color","black",'MarkerFaceColor','black','MarkerSize',5)
hold on;

f7=figure;
if dimension==2
    u=zeros(201,201);
    for i=1:201
        for j=1:201
            u(j,i)=RanaFun([-505+i*5,-505+j*5]);    % opposé à ce qu'on penserait
        end
    end
    [X_plot,Y_plot]=meshgrid(-500:5:500,-500:5:500);
    mesh(X_plot,Y_plot,u)
    hold on;
end
TestX=X_parents_history{1,5};
TestY=Y_parents_history{1,5};
plot3(TestX(1,:),TestX(2,:),TestY,'o',"color","black",'MarkerFaceColor','black','MarkerSize',5)

f8=figure;
if dimension==2
    u=zeros(41,41);
    for i=1:41
        for j=1:41
            u(j,i)=RanaFun([-480.025+i*0.025,499.475+j*0.025]);    % opposé à ce qu'on penserait
        end
    end
    [X_plot,Y_plot]=meshgrid(-480:0.025:-479,499.5:0.025:500.5);
    mesh(X_plot,Y_plot,u)
    hold on;
end
TestX=X_parents_history{1,20};
TestY=Y_parents_history{1,20};
plot3(TestX(1,:),TestX(2,:),TestY,'.',"color","black",'MarkerFaceColor','black','MarkerSize',15)


