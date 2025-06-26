function neighbours = findneighbours(indx,DT)
        B= DT(any(DT==indx,2),:);
        C= unique(reshape(B,[],1));
        neighbours=C(C~=indx);
end