function plotter(plotfigure,CC,i_xcoord,i_ycoord,R,i_id,i_type,v_LESC,t)

if plotfigure == 1

        [vMGFPoints, vMGFcells] = voronoin(CC(:,[i_xcoord i_ycoord]));
         
        [vMGFx, vMGFy] = voronoi(CC(:,i_xcoord),CC(:,i_ycoord));

       
        hold off
        plot(R*cos(2*pi*(0:0.01:1)), R*sin(2*pi*(0:0.01:1)), 'b', 'LineWidth', 2)
        hold on


        plot(vMGFx,vMGFy,'b','LineWidth', 2)  % Comment this line for empty cornea
        axis(R*[-1 1 -1 1])

        colors = [1   1   1;    % white
            0   0   0;    % black
            0.7 0.7 0.7;  % grey
            1   0   1;    % magenta
            0   1   0;    % green
            0   0   1;    % blue
            0.5 0.5 1];   % light blue

        number_of_cells = size(CC, 1);
        for i_cell = 1:number_of_cells

            currentcolor = colors(mod(CC(i_cell,i_id),7)+1,:);

            X = vMGFPoints(vMGFcells{i_cell},1);
            Y = vMGFPoints(vMGFcells{i_cell},2);
            X = X(X ~= Inf);
            Y = Y(Y ~= Inf);
            patch(X, Y, currentcolor)

            if CC(i_cell,i_id) == 0
                plot(CC(i_cell,i_xcoord), CC(i_cell,i_ycoord),'kx', 'MarkerSize', 6)
            end

        end

        %--------------------------------------------------------------
        % Make pretty boundary around disc
        theta = pi*(0:0.01:0.5);
        X = R*[cos(theta) 1];
        Y = R*[sin(theta) 1];
        patch(X, Y, 'w')

        theta = pi*(0.5:0.01:1);
        X = R*[cos(theta) -1];
        Y = R*[sin(theta) 1];
        patch(X, Y, 'w')

        theta = pi*(1:0.01:1.5);
        X = R*[cos(theta) -1];
        Y = R*[sin(theta) -1];
        patch(X, Y, 'w')

        theta = pi*(1.5:0.01:2);
        X = R*[cos(theta) 1];
        Y = R*[sin(theta) -1];
        patch(X, Y, 'w')

        theta = pi*(0:0.01:2);
        plot(R*cos(theta), R*sin(theta), 'b', 'Linewidth', 2)
        %--------------------------------------------------------------

        for i_cell = 1:number_of_cells
            if CC(i_cell,i_type) == v_LESC
                currentcolor = colors(mod(CC(i_cell,i_id),7)+1,:);
                plot(CC(i_cell,i_xcoord), CC(i_cell,i_ycoord), 'bo', 'Linewidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', currentcolor)
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %--------------------------------------------------------------

        v = axis;
        counter = ['t = ', num2str(t)];

        text(v(1) + 0.8*(v(2) - v(1)), v(3) + 0.93*(v(4) - v(3)), counter,'FontSize', 10 , 'FontWeight','bold')
        
        set(gca, 'Fontsize', 15,'FontWeight','bold')
        set(gca,'linewidth', 1.5,'Fontsize', 15,'FontWeight','bold')
        
         getframe;
         
else
    t
end
end

