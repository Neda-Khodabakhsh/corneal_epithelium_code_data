function value = F_pull(i_cell,DT,B,A,i_xcoord,i_ycoord,restlength,c_pull,v_TAC,i_type) 

value = 0;
        
R = 500 ; 

        for j_cell = 1: height(B)

            h = 15;
            
               if (B(j_cell,i_xcoord)-A(i_cell,i_xcoord))^2 + (B(j_cell,i_ycoord)-A(i_cell,i_ycoord))^2 < 300

                     dpush = sqrt((B(j_cell,i_xcoord)-A(i_cell,i_xcoord)).^2 + (B(j_cell,i_ycoord)-A(i_cell,i_ycoord)).^2 + h^2 ); 
                 
                     value = value + (h*c_pull*(dpush- restlength )/dpush);
              
                end
       end
              
end




 