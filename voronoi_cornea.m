function   [layer] = voronoi_cornea

rng(1,"twister");
%--------------------------------------------------------------------------
tfinal = 500;  
plotfigure = 1; 
n_layer  = 2;

% Parameters
n_max    = 3;  %
restlength = 15;   %
lambda_TAC = 0.055; %   
dt = 5;  %5

dist_newcell = 0.05; %
% r_min = 20; 
T_LESC = 48; %
T_TAC = 24;    %

FT_divthreshold = 3; % 1 0.005
FL_divthreshold = 5; % 1 0.001
c_adhesive = 0.005;   %0.005 0.004

c_squeeze = 0.1; %0.1
c_pull = 0.7; %0.7

c_compression = 1; %0.1
c_compressionL = 1;

rshed = 0.002;  %0.002

% Define indices
i_xcoord = 1;
i_ycoord = 2;
i_type   = 3;
i_id     = 4;
i_nd     = 5;
i_clock  = 6;  
i_label =  7;

% Define values of cell types
v_boundary = 0;
v_TAC  = 1;
v_LESC = 2;
v_extruded = -9;
v_shedding = -99;

v_dead = -10;
v_divideLESC = -22;
v_divideTAC = -11;
%--------------------------------------------------------------------------
% Generate disc
%--------------------------------------------------------------------------
R = 500;                  
N_rim = 50;                          
N_interior_basal  = 0;  
N_interior_layern = 0;   

% Generate basal layer
%--------------------------------------------------------------------------
% Points on the rim of basal layer
X_basal = R*cos(2*pi*(0:N_rim-1)'/N_rim);
Y_basal = R*sin(2*pi*(0:N_rim-1)'/N_rim);

% Randomly select interior points of basal layer
n = 0;
while n < N_interior_basal
    new_x = R*(2*rand - 1);
    new_y = R*(2*rand - 1);
    if new_x^2 + new_y^2 <= (R)^2
        n = n + 1;
        X_basal = [X_basal; new_x];
        Y_basal = [Y_basal; new_y];
    end
end
% Set cell types - LESCs on rim of basal layer, all others are TACs
celltypes_basal  = zeros(length(X_basal),1);
celltypes_basal(1:N_rim) = v_LESC;
celltypes_basal((N_rim+1):(N_rim+N_interior_basal)) = v_TAC;
% Set ID for each LESC
ID_basal  = zeros(length(X_basal),1);
ID_basal(1:N_rim) = 1:N_rim;    % All LESCs get an ID 1, 2, 3, ...
% Set numbers of division to 0
ND_basal = zeros(length(X_basal),1);
clock_basal = zeros(length(X_basal),1);

label_basal = zeros(length(X_basal),1);

for i = 1: 50
label_basal(i) = -i;
end

layer(1).cells = [X_basal   Y_basal   celltypes_basal   ID_basal   ND_basal  clock_basal label_basal ];  % basal layer

%--------------------------------------------------------------------------
% Generate upper layers
%--------------------------------------------------------------------------
for ii=2:n_layer

    %Points on the rim of layer n % They are boundary points not LESCs
    X_layern = R*cos(2*pi*(0:N_rim-1)'/N_rim);
    Y_layern = R*sin(2*pi*(0:N_rim-1)'/N_rim);

    % Randomly select interior points of layer n
    n = 0;
    while n < N_interior_layern
        new_x = R*(2*rand - 1);
        new_y = R*(2*rand - 1);
        if new_x^2 + new_y^2 <= (R)^2
            n = n + 1;
            X_layern = [X_layern; new_x];
            Y_layern = [Y_layern; new_y];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    celltypes_layern = zeros(length(X_layern),1);
    celltypes_layern(1:N_rim) = v_boundary;
    celltypes_layern((N_rim+1):(N_rim+N_interior_layern)) = v_TAC;
    ID_layern = zeros(length(X_layern),1);
    ND_layern = zeros(length(X_layern),1);
    clock_layern = zeros(length(X_layern),1);
    
    label_layern = zeros(length(X_layern),1);
    %----------------------------------------------------------------------
    layer(ii).cells = [X_layern  Y_layern  celltypes_layern  ID_layern  ND_layern  clock_layern  label_layern ]; % second layer
   
end

% figure

for t = dt:dt:tfinal
     
    % Find F_adhesive for all cells in layer 1
     F_adhesive = adhesion(layer(1).cells,i_xcoord,i_ycoord,c_adhesive);
      
 for s=1:n_layer  

     XY = [layer(s).cells(:,i_xcoord), layer(s).cells(:,i_ycoord)]; 
     DT = delaunay(XY);
   
     number_of_cells = size(layer(s).cells, 1); 

        for i_cell = 1:number_of_cells

          if layer(s).cells(i_cell,i_type)==v_TAC    

                if (layer(s).cells(i_cell,i_nd) == n_max) && (rand < 0)   %0.008
                    layer(s).cells(i_cell,i_type) = v_dead;  

                elseif (s==2) && (rand < rshed)
                
                    layer(s).cells(i_cell,i_type) = v_shedding; 

                elseif s==1   
                         
                    if F_squeeze(i_cell,DT,layer(s).cells,i_xcoord,i_ycoord,c_squeeze,restlength)+F_pull(i_cell,DT,layer(s+1).cells,layer(s).cells,i_xcoord,i_ycoord,20,c_pull,v_TAC,i_type) > F_adhesive(i_cell)
                           

                     layer(s).cells(i_cell,i_type) = v_extruded;  
                                     
                           
                    elseif (layer(s).cells(i_cell,i_clock) >= T_TAC)&&(layer(s).cells(i_cell,i_nd) < n_max)&&(F_compression(i_cell,DT,layer(s).cells,i_xcoord,i_ycoord,c_compression,restlength) < FT_divthreshold )
                       
                     layer(s).cells(i_cell,i_type) = v_divideTAC;  % Mark TAC as "to divide"

                   end
                end
               
            elseif layer(s).cells(i_cell,i_type)==v_LESC    

                if (layer(s).cells(i_cell,i_clock) >= T_LESC)&&(F_compression(i_cell,DT,layer(s).cells,i_xcoord,i_ycoord,c_compressionL,20) < FL_divthreshold )
                    
                   layer(s).cells(i_cell,i_type) = v_divideLESC;
                  
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
          end
        end
        %------------------------------------------------------------------
        % Remove all dead cells
        indx = (layer(1).cells(:,i_type) == v_dead);
        layer(1).cells(indx,:) = [];
        
        %------------------------------------------------------------------
        % Delaminate all delaminated cells
        indx =  (layer(1).cells(:,i_type) == v_extruded);
        indxshed =  (layer(2).cells(:,i_type) == v_shedding);
        
      
         if s == 1

           layer(s).cells(indx,i_type) = v_TAC;   % Restore cell type to TAC
           layer(s+1).cells = [layer(s+1).cells;    layer(s).cells(indx,:)];
           layer(s).cells(indx,:) = [];  % remove cells from current layer

         end
       
         if s == 2
        % 
            layer(s).cells(indxshed,i_type) = v_TAC;   
            layer(s).cells(indxshed,:) = [];  
        % 
         end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  s==1     % Division only applies to layer 1

            indx_LESC = (layer(s).cells(:,i_type) == v_divideLESC);
            indx_TAC  = (layer(s).cells(:,i_type) == v_divideTAC);
            indx      = indx_LESC|indx_TAC;
            % indx_num  = find(indx);
            indx_num_LESC  = find(indx_LESC);
            indx_num_TAC  = find(indx_TAC);
            
            newTACs_L = zeros(length(indx_num_LESC), 7);
            for n = 1:length(indx_num_LESC)
                i_dividingcell = indx_num_LESC(n);

                newX = 2*(R); newY = 0;   % Set newX and newY to some illegal numbers so that the while loop will run
                while (newX.^2 + newY.^2) > R.^2
                    theta = 2*pi*rand;
                    newX = layer(s).cells(i_dividingcell,i_xcoord) + dist_newcell*cos(theta);
                    newY = layer(s).cells(i_dividingcell,i_ycoord) + dist_newcell*sin(theta);
                    
                end
                %                 X       Y      type        lineage ID                          ND   clock
                newTACs_L(n, :) = [newX    newY    v_TAC    layer(s).cells(i_dividingcell,i_id)    0    0   layer(s).cells(i_dividingcell,i_label)];
 
            end
            layer(s).cells = [layer(s).cells;    newTACs_L];
            %------------------------------------------------------------------------
            newTACs = zeros(length(indx_num_TAC), 7);
            for n = 1:length(indx_num_TAC)
                i_dividingcell = indx_num_TAC(n);

                newX = 2*(R); newY = 0;   % Set newX and newY to some illegal numbers so that the while loop will run
                while (newX.^2 + newY.^2) > R.^2
                    theta = 2*pi*rand;
                    newX = layer(s).cells(i_dividingcell,i_xcoord) + dist_newcell*cos(theta);
                    newY = layer(s).cells(i_dividingcell,i_ycoord) + dist_newcell*sin(theta);
                    
                end    
                %                 X       Y      type        lineage ID                          ND                                      clock
                newTACs(n, :) = [newX    newY    v_TAC    layer(s).cells(i_dividingcell,i_id)    layer(s).cells(i_dividingcell,i_nd)+1     0     layer(s).cells(i_dividingcell,i_label) ];
 
            end
               layer(s).cells = [layer(s).cells;    newTACs];
         %----------------------------------------------------------------------------------------
               layer(s).cells(~indx, i_clock) = layer(s).cells(~indx, i_clock) + dt;

               layer(s).cells(indx_LESC, i_type)   = v_LESC;
               layer(s).cells(indx_TAC, i_type)      = v_TAC;

               layer(s).cells(indx, i_nd)     = layer(s).cells(indx, i_nd) + 1;
               layer(s).cells(indx, i_clock)   = 0;


               layer(s).cells(indx_LESC, i_label)     = 0;
               layer(s).cells(indx_TAC, i_label)      = 0;
        end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Movement
    for k =1:7

        XY = [layer(s).cells(:,i_xcoord), layer(s).cells(:,i_ycoord)]; 
        DT = delaunay(XY);

        indx = (layer(s).cells(:,i_type) == v_TAC);   % Only TACs move
        indx_num = find(indx);
        newpositions = zeros(length(indx_num), 2);
        for n = 1:length(indx_num)
            i_cell = indx_num(n);

            neighbours = findneighbours(i_cell,DT);

            displacementvector = [0 0];
            for num = 1:length(neighbours)
                indx_neighbour = neighbours(num);
                v = [layer(s).cells(indx_neighbour,i_xcoord)  layer(s).cells(indx_neighbour,i_ycoord)] - [layer(s).cells(i_cell,i_xcoord)  layer(s).cells(i_cell,i_ycoord)];
                lengthv = sqrt(v(1)^2 + v(2)^2);
                displacementvector = displacementvector + dt*lambda_TAC*(v/lengthv)*(lengthv - restlength);
            end
            newX = layer(s).cells(i_cell,i_xcoord) + displacementvector(1);  % Update X and Y coordinates
            newY = layer(s).cells(i_cell,i_ycoord) + displacementvector(2);

            % Make sure cells stop at disc boundary
            if newX^2 + newY^2 > R^2
                a = displacementvector(1)^2 + displacementvector(2)^2;
                b = 2*(displacementvector(1)*layer(s).cells(i_cell,i_xcoord) + displacementvector(2)*layer(s).cells(i_cell,i_ycoord));
                c = layer(s).cells(i_cell,i_xcoord)^2 + layer(s).cells(i_cell,i_ycoord)^2 - R^2;
                k = (-b + sqrt(b^2 - 4*a*c))/(2*a);
                newX = layer(s).cells(i_cell,i_xcoord) + k*displacementvector(1);  % Update X and Y coordinates
                newY = layer(s).cells(i_cell,i_ycoord) + k*displacementvector(2);
            end
            newpositions(n,:) = [newX    newY];

        end
        % Update all cell positions
        layer(s).cells(indx,[i_xcoord i_ycoord]) = newpositions;

    end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot
for s=1:n_layer
% daspect([1, 1, 1]);
% daspect([1, 1, 1]);
         subplot(1,2,s)

   plotter(plotfigure,layer(s).cells,i_xcoord,i_ycoord,R,i_id,i_type,v_LESC,t);
end

end

% Plot final Voronoi diagram
   plotfigure = 1;

for s=1:n_layer

% daspect([1, 1, 1]);
% daspect([1, 1, 1]);
         subplot(1,2,s)
        
   plotter(plotfigure,layer(s).cells,i_xcoord,i_ycoord,R,i_id,i_type,v_LESC,t);
end
end






