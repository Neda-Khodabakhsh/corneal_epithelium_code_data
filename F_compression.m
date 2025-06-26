function  val = F_compression(A,DT,CC,i_xcoord,i_ycoord,c_compression,restlength)
    

    neighbours = findneighbours(A,DT);            
    distances = sqrt( (CC(neighbours,i_xcoord) - CC(A,i_xcoord)).^2 + (CC(neighbours,i_ycoord) - CC(A,i_ycoord)).^2 );
    new = (restlength - distances);
   
    val = max(c_compression*new);
     
end