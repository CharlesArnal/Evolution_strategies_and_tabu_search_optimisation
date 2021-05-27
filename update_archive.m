


function arch=update_archive(arch,X,Y,arch_size,D1,D2)
	[Y_arch,J]=sort(arch{2},"descend");
    X_arch=arch{1}(:,J);
    for j=1:size(X,2)
        xj=X(:,j);
        yj=Y(:,j);
        n=size(X_arch,2);
        distances=zeros(n);
        for l=1:n
            distances(l)=norm(xj-X_arch(:,l));
        end
        [sorted_distances,I]=sort(distances,"ascend");
        if n==0 | sorted_distances(1)>D1  % dissimilar enough to the archive
            if n<arch_size  % archive not yet full -> include
                X_arch=[X_arch,xj];
                Y_arch=[Y_arch,yj];
            elseif yj>Y_arch(n)   % better than the worst -> replace it
                X_arch(:,n)=xj;
                Y_arch(n)=yj;
            end
        else
            if yj>Y_arch(1)  % best solution found so far -> replace closest point
                X_arch(:,I(1))=xj;
                Y_arch(I(1))=yj;
            elseif sorted_distances(1)<D2 & yj>Y_arch(I(1))
                % close enough and better than the closest point-> replace it
                X_arch(:,I(1))=xj;
                Y_arch(I(1))=yj;
            end
        end
        % reorder archive
        [Y_arch,J]=sort(Y_arch,"descend");
        X_arch=X_arch(:,J);
    end
    arch={X_arch,Y_arch};
end
