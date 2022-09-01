function [X2] = norm_coeff(X,lb,ub,iopt)

nsize=max(size(X));
X2=zeros(size(X));
indx=zeros(size(X));
db=ub-lb;
val=0;
% Define index for optimization
for i=1:nsize
    if db(i)~=0
        indx(i)=i;
    end
end

if iopt==1 % Raw data -> Normalized data
    for i=1:nsize
        if indx(i)~=0
            k=indx(i);
            X2(k)=(X(k)-lb(k))/(ub(k)-lb(k));
        else
            X2(i)=X(i);
        end
    end
else        % Normalized data -> Raw data
    for i=1:nsize
        if indx(i)~=0
            k=indx(i);
            X2(k)=X(k)*(ub(k)-lb(k))+lb(k);
        else
            X2(i)=X(i);
        end
    end
end
end

