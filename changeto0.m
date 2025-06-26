load("layersun.mat")


for i = 51:size(layersun(1).cells,1)
 layersun(1).cells(i,4)= 0;

end

for i = 51:size(layersun(2).cells,1)
 layersun(2).cells(i,4)= 0;

end


for i = 1:50
 layersun(1).cells(i,5)= 0;

end






