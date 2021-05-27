close all;
clear all;

% Tabu Search
% --------------------------

dimension=2;
range=500;
gridRatio=3;
stepsize=200;
stepsize_red_coeff=0.9;
STM_size=7;
MTM_size=5;
intensify_thres=15;
diversify_thres=25;
reduce_thres=30;
numMaxEvaluation=10000;
tolerance=0.000001;
ConvCriterion="temp";
tabu_direction=false;
STM_direction_size=10;
c_wanderlust=1;
concentric=false;
minOrMax="min";
Cf={@constraintFunctionRanaProblem};
Of=@RanaFun;
archive_size=20;
D1=20;
D2=1;



%---------
M=1;

time=zeros(1,M);
best_performance=zeros(1,M);
mean_performance=zeros(1,M);
performance_std_dev=zeros(1,M);
locateBest=zeros(5,100);%erase
for j=1:M
NumRuns=1;
Min_run=zeros(1,NumRuns);
tic
for i=1:NumRuns
    rng(i);
    [X_history,Y_history,archive,graphdata]=Tabu(Cf,Of,dimension,range,gridRatio,stepsize,stepsize_red_coeff,...
        STM_size,tabu_direction,STM_direction_size,c_wanderlust,MTM_size,intensify_thres,diversify_thres,...
        reduce_thres,numMaxEvaluation,concentric,tolerance,minOrMax,archive_size,D1,D2);
        Min_run(i)=archive{2}(1);
        locateBest(:,i)=archive{1}(1);%erase
    i
end
time(j)=toc
best_performance(j)=min(Min_run)
mean_performance(j)=mean(Min_run)
performance_std_dev(j)=std(Min_run)
[sortedMin,I]=sort(Min_run,"ascend");%erase
best_Y=Min_run(I(1));%erase
best_X=locateBest(:,I(1));%erase
end

    
f1=figure;
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
plot3(X_history(1,:),X_history(2,:),Y_history,'o',"color","black",'MarkerFaceColor','black','MarkerSize',5)
hold on;

TestX=archive{1};
TestY=archive{2};
plot3(TestX(1,:),TestX(2,:),TestY,'o',"color","red",'MarkerFaceColor','red','MarkerSize',7.5)
hold on;

 
f2=figure;
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
plot3(X_history(1,1:100),X_history(2,1:100),Y_history(1:100),'o',"color","black",'MarkerFaceColor','black','MarkerSize',5)
hold on;


f3=figure;
if dimension==2
    u=zeros(101,101);
    for i=1:101
        for j=1:101
            u(j,i)=RanaFun([-305.1+i*0.1,494.9+j*0.1]);    % opposé à ce qu'on penserait
        end
    end
    [X_plot,Y_plot]=meshgrid(-305:0.1:-295,495:0.1:505);
    mesh(X_plot,Y_plot,u)
    hold on;
end
testx1=X_history(1,2000:2100);
I=(testx1>-305 & testx1<-295)
testx1=testx1(I);
testx2=X_history(2,2000:2100);
testx2=testx2(I);
testy=Y_history(2000:2100);
testy=testy(I);

plot3(testx1,testx2,testy,'o',"color","black",'MarkerFaceColor','black','MarkerSize',5)
hold on;