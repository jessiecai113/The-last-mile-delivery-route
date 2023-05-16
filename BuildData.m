numOfPoints = 401;                                                           % Number of points（DC+customers）

coordinateOfPoints = rand(numOfPoints, 2) * (100 - 0) + 0;                  % generate coordinate
demandOfPoints = rand(numOfPoints, 1) * (100 - 0) + 1;                      % generate demand


%% save file
fileName = sprintf('./data/vrpData%04d.txt', numOfPoints);
fid = fopen(fileName, 'w');
for i = 1 : numOfPoints
	fprintf(fid,'%.2f %.2f %.1f\n', coordinateOfPoints(i, 1), coordinateOfPoints(i, 2), demandOfPoints(i));
end
fclose(fid);