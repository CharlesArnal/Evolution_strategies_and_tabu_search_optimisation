
% The main function (parameters are described in the report)
function [X_history,Y_history,archive,graphdata]=Tabu(Cf,Of,dimension,range,gridRatio,...
    stepsize,stepsize_red_coeff,STM_size,tabu_direction,STM_direction_size,c_wanderlust,...
    MTM_size,intensify_thres,diversify_thres,reduce_thres,numMaxEvaluation,concentric,...
    tolerance,minOrMax,archive_size,D1,D2)
    if minOrMax=="min"
        Of = @(x) -Of(x);
    end
    evaluations=0;
    moves=0;
    counter=0;
    X_history=zeros(dimension,0);
    Y_history=zeros(1,0);
    archive={zeros(dimension,0),zeros(1,0)};
    areas=initialize_areas(dimension,gridRatio);
    
    graphdata=zeros(2,0);

    STM=NaN(dimension,0);
    MTM={zeros(1,0),zeros(dimension,0)};
    
    [x,areas]=GeneratePoint(range,dimension,gridRatio,Cf,areas);
    STM=[x];
    y=Of(x);
    evaluations=evaluations+1;
    center=x;
    [moves,areas,X_history,Y_history,MTM,newOptimum,archive]=general_update(x,y,...
        range,dimension,gridRatio,MTM,MTM_size,moves,areas,...
        X_history,Y_history,archive,archive_size,D1,D2);
    
    
    % the main loop
    while evaluations<=numMaxEvaluation-(2*dimension+1) &  stepsize>=tolerance    
        [x,y,STM,evaluations,trapped]=move(x,y,STM,tabu_direction,STM_size,STM_direction_size,c_wanderlust,Of,Cf,stepsize,evaluations,concentric,center);
        [moves,areas,X_history,Y_history,MTM,newOptimum,archive]=general_update(x,y,...
            range,dimension,gridRatio,MTM,MTM_size,moves,areas,...
            X_history,Y_history,archive,archive_size,D1,D2);
        
        if newOptimum
            counter=0;
            center=x;
        else
            counter=counter+1;
        end
        
        if trapped & concentric & counter<intensify_thres
            counter=intensify_thres;
        end
        
        if counter==intensify_thres & feasible(mean(MTM{2},2),Cf)
            x=mean(MTM{2},2);
            STM=[x];
            y=Of(x);
            center=x;
            evaluations=evaluations+1;
            [moves,areas,X_history,Y_history,MTM,newOptimum,archive]=general_update(x,y,...
                range,dimension,gridRatio,MTM,MTM_size,moves,...
                areas,X_history,Y_history,archive,archive_size,D1,D2);
            if newOptimum
                counter=0;
            else
                counter=counter+1;
            end
        end
        if counter==diversify_thres
            [x,areas]=GeneratePoint(range,dimension,gridRatio,Cf,areas);
            STM=[x];
            y=Of(x);
            evaluations=evaluations+1;
            center=x;
            [moves,areas,X_history,Y_history,MTM,newOptimum,archive]=general_update(x,y,...
                range,dimension,gridRatio,MTM,MTM_size,moves,areas,...
                X_history,Y_history,archive,archive_size,D1,D2);
        end
        if counter==reduce_thres
            stepsize=stepsize*stepsize_red_coeff;
            x=MTM{2}(:,end);
            y=MTM{1}(:,end);
            STM=[x];
            center=x;
            [moves,areas,X_history,Y_history,MTM,newOptimum,archive]=general_update(x,y,...
                range,dimension,gridRatio,MTM,MTM_size,moves,areas,...
                X_history,Y_history,archive,archive_size,D1,D2);
            counter=0;
        end
        % specifically to create some graphs
        graphdata=[graphdata,[evaluations;-archive{2}(1)]];
    end
    evaluations
    
    if minOrMax=="min"
        Y_history=-Y_history;
        archive{2}=-archive{2};
    end
end

function [moves,areas,X_history,Y_history,MTM,newOptimum,archive]=general_update(x,y,...
    range,dimension,gridRatio,MTM,MTM_size,moves,areas,X_history,...
    Y_history,archive,archive_size,D1,D2)
    areas=update_areas(x,areas,range,dimension,gridRatio);
    moves=moves+1;
    X_history(:,moves)=x;
    Y_history(:,moves)=y;
    [MTM,newOptimum]=updateMTM(x,y,MTM,MTM_size);
    archive=update_archive(archive,x,y,archive_size,D1,D2);
end



function areas=initialize_areas(dimension,gridRatio)
    numberMoves=zeros(1,gridRatio^dimension);
    indices=zeros(dimension,gridRatio^dimension);
    if gridRatio~=1
        for i=1:gridRatio^dimension
            indices(:,i)=de2bi(i-1,dimension,gridRatio);
        end
    else
       indices=0; 
    end
    areas={numberMoves,indices};
end

% could recover dimension and gridRatio from areas, but makes the code
% messier
% works under the assumption that there are feasible points in each area
% - (it's the case here)
function [x,areas]=GeneratePoint(range,dimension,gridRatio,Cf,areas)
    if sum(areas{1}==0)==0
        areas{1}=zeros(1,gridRatio^dimension);
    end
    candidate_areas=areas{1}(areas{1}==0);
    candidate_areas_indices=areas{2}(:,areas{1}==0);
    acceptable=false;
    while not(acceptable)
        i=randi(numel(candidate_areas));
        index=candidate_areas_indices(:,i);        
        x=(2*range/gridRatio)*(rand(dimension,1)+index)-range;
        acceptable= feasible(x,Cf);    % Redundant with the previous line in our special case
    end
end

function areas=update_areas(x,areas,range,dimension,gridRatio)
    index=floor((x+range)/(2*range/gridRatio));
    if gridRatio~=1;
        i=bi2de(index',dimension,gridRatio)+1;
    else
        i=1;
    end
    areas{1}(i)=areas{1}(i)+1;
end

%to be modified trouver min pas max
function [x,y,STM,evaluations,trapped]=move(x,y,STM,tabu_direction,STM_size,STM_direction_size,c_wanderlust,Of,Cf,stepsize,evaluations,concentric,center)
    d=size(x,1);
    neighb_x=repmat(x,1,2*d)+stepsize*[eye(d),-eye(d)];
    neighb_y=zeros(1,2*d);
    tabu=true(1,2*d);
    parfor i=1:2*d
        if not(concentric)
            tabu(1,i)=not(feasible(neighb_x(:,i),Cf)) | ( ismember(neighb_x(:,i)',STM(:,1:min(size(STM,2),STM_size))','rows') ) ;
        else
            tabu(1,i)=not(feasible(neighb_x(:,i),Cf))| norm(center-neighb_x(:,i))<norm(center-x)
        end
            
        if not(tabu(1,i))
            neighb_y(1,i)=Of(neighb_x(:,i)');
            evaluations=evaluations+1;
        end   
    end
    
    if sum(not(tabu))==0        % This is a rare occurence for concentric==false
        "Point is trapped"      
        STM=[x];
        trapped=true;
        return
    else
        trapped=false;
    end
    
    feasible_x=neighb_x(:,not(tabu));
    feasible_y=neighb_y(1,not(tabu));
    %[sorted_y,I]=sort(feasible_y,'descend');
    [sorted_y,I]=selectBestMove(feasible_y,feasible_x,x,STM,tabu_direction,STM_direction_size,c_wanderlust);
    x_pattern_move=2*feasible_x(:,I(1))-x;
    y=sorted_y(1);
    x=feasible_x(:,I(1));
    STM=update_STM(x,STM,max(STM_size,STM_direction_size));
    if feasible(x_pattern_move,Cf)     % ask that it be not in STM?
        y_pattern_move=Of(x_pattern_move');
        evaluations=evaluations+1;
        if y_pattern_move>y
            y=y_pattern_move;
            x=x_pattern_move;
            STM=update_STM(x,STM,max(STM_size,STM_direction_size));
        end        
    end
end

function STM=update_STM(x,STM,n)
    if size(STM,2)<n
        STM=[x,STM];
    else
        STM=circshift(STM,1,2);
        STM(:,1)=x;
    end
end

function     [sorted_y,I]=selectBestMove(feasible_y,feasible_x,x,STM,tabu_direction,STM_direction_size,c_wanderlust)
    if tabu_direction==false | size(STM,2)<=1
        [sorted_y,I]=sort(feasible_y,'descend');
    else
        past_direction=-(x-mean(STM(:,1:min(end,STM_direction_size)),2))/norm(x-mean(STM(:,1:min(STM_direction_size,end)),2));
        
        n=size(feasible_x,2);
        directional_aversion=zeros(1,n);
        step_size=norm(feasible_x(:,1)-x);
        % to be modified with min instead of max
        for i=1:n
            directional_aversion(1,i)=-((feasible_x(:,i)-x)/step_size)'*past_direction;
        end
        
        modified_y=feasible_y+std(feasible_y,1)  *c_wanderlust* directional_aversion;
        [sorted_modified_y,I]=sort(modified_y,'descend');
        sorted_y=feasible_y(I);
    end
end

% MTM constantly sorted in ascending order
function [MTM,newOptimum]=updateMTM(x,y,MTM,MTM_size)
    if size(MTM{1},2)==0 | y>MTM{1}(end)
        newOptimum=true;
    else
        newOptimum=false;
    end
    if size(MTM{1},2)<MTM_size
        [Y_MTM,I]=sort([y,MTM{1}]);
        temp_X=[x,MTM{2}];
        X_MTM=temp_X(:,I);
        MTM={Y_MTM,X_MTM};
    elseif y>MTM{1}(1,1)
        [Y_MTM,I]=sort([y,MTM{1}(:,2:end)]);
        temp_X=[x,MTM{2}(:,2:end)];
        X_MTM=temp_X(:,I);
        MTM={Y_MTM,X_MTM};
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
