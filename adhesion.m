function F_adhesive = adhesion(A,i_xcoord,i_ycoord,c_adhesive)
% nargin ('adhesion')

    [vMGFPoints, vMGFcells] = voronoin([A(:,i_xcoord), A(:,i_ycoord)]);
    
    % number_of_cells = size(layer(1).cells, 1)
    % AREA = zeros(length(vMGFcells),1);
    F_adhesive = zeros(length(vMGFcells),1); 

    for ii_cell = 1:length(vMGFcells)

        X = vMGFPoints(vMGFcells{ii_cell},1);
        Y = vMGFPoints(vMGFcells{ii_cell},2);
        X = X(X ~= Inf);
        Y = Y(Y ~= Inf);

          F_adhesive(ii_cell) = c_adhesive*polyarea(X,Y);  %AREA(ii_cell);
    end
        
end