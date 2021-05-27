close all;
clear all;

dimension=5;
scale=50;
range=500;
N=1000;
numMaxEvaluations=10000;
lambda=200;
mu=25;
minOrMax="min";
tolerance=0.000001;
ContrainStrat="accept";
ConvCriterion="RelDiff";
Cf={@constraintFunctionRanaProblem};
%Cf={@(x) 1};
Of=@RanaFun;
%Of=@(x) norm(x);
archive_size=20;
D1=20;
D2=1;


lambda=[10,30,100,220,450,900];
mu=lambda/2;
M=6;


time=zeros(1,M);
best_performance=zeros(1,M);
mean_performance=zeros(1,M);
performance_std_dev=zeros(1,M);

for j=1:M
NumRuns=10;
Min_run=zeros(1,NumRuns);
tic
for i=1:NumRuns
    rng(i);
    [X_parents_history,Y_parents_history,archive]=CMA_ES(Cf,Of,dimension,range,N,...
    scale,numMaxEvaluations,mu(j),lambda(j),tolerance,ConvCriterion,ContrainStrat,minOrMax,archive_size,D1,D2);
    
    Min_run(i)=archive{2}(1);
    i
end
time(j)=toc
%best_performance(j)=min(Min_run)
mean_performance(j)=mean(Min_run)
performance_std_dev(j)=std(Min_run)
end


%{
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
            u(j,i)=Of([-505+i*5,-505+j*5]);    % opposé à ce qu'on penserait
        end
    end
    [X_plot,Y_plot]=meshgrid(-500:5:500,-500:5:500);
    mesh(X_plot,Y_plot,u)
    hold on;
end

TestX=archive{1};
TestY=archive{2};
plot3(TestX(1,:),TestX(2,:),TestY,'o',"color","red")
hold on;
%}
%{
TestX=X_parents_history{1,1};
TestY=Y_parents_history{1,1};
plot3(TestX(1,:),TestX(2,:),TestY,'x',"color","red")
hold on;

TestX=X_parents_history{1,7};
TestY=Y_parents_history{1,7};
plot3(TestX(1,:),TestX(2,:),TestY,'x',"color","blue")

TestX=X_parents_history{1,15};
TestY=Y_parents_history{1,15};
plot3(TestX(1,:),TestX(2,:),TestY,'x',"color","black")
%}
