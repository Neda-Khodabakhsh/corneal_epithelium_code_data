function  [diameter_TAC] = diameterfinder(layer)

i_xcoord = 1;
i_ycoord = 2;
i_type   = 3;
v_TAC    = 1; 


[vMGFPoints, vMGFcells] = voronoin([layer(2).cells(:,i_xcoord), layer(2).cells(:,i_ycoord)]);
    
   
    for ii_cell = 1:length(vMGFcells)

        X = vMGFPoints(vMGFcells{ii_cell},1);
        Y = vMGFPoints(vMGFcells{ii_cell},2);
        X = X(X ~= Inf);
        Y = Y(Y ~= Inf); 

        if X < 500
            if -500 < X
                if Y < 500
                  if -500 < Y
        
                      if (X.^2 + Y.^2 < 500^2)
                          if layer(2).cells(ii_cell,i_type)==v_TAC
          
                             diameter_TAC(ii_cell) = 2*(sqrt(polyarea(X,Y)/pi));
                          end
        
                     end 

                 end
        
               end
            end
        end
    end 
  
end

