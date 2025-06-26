function  val = F_squeeze(A,DT,CC,i_xcoord,i_ycoord,c_squeeze,restlength)
    neighbours = findneighbours(A,DT);            
    distances = sqrt( (CC(neighbours,i_xcoord) - CC(A,i_xcoord)).^2 + (CC(neighbours,i_ycoord) - CC(A,i_ycoord)).^2 );
    
    new = (restlength - distances);
    
    val = (c_squeeze*sum(new));

end