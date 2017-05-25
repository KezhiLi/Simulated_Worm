addpath('C:\Kezhi\Papers\NIPS2017\codes\drtoolbox\techniques\');

% first manually import data

% remove nan rows
badInds = isnan(featMat2(:, 1));
featMat2(badInds, :) = [];
wormNames2(badInds) = [];

no_dims = 2;
initial_dims = 7;
[uniqueNames, ~, numLabels] = unique(wormNames2);
ydata = tsne(featMat2, numLabels, no_dims, initial_dims);


% plot the points from each compound
figure
for ii = 1:numel(uniqueNames)
    currentInds = strcmp(wormNames2, uniqueNames{ii});
        plot(ydata(currentInds, 1), ydata(currentInds, 2), '.', 'MarkerSize', 50)
     
    %    plot3(ydata(currentInds, 1), ydata(currentInds, 2), ydata(currentInds, 3), '.', 'MarkerSize', 50)
    hold on
end
lgd = legend(uniqueNames{1},uniqueNames{2},uniqueNames{3},uniqueNames{4},uniqueNames{5},uniqueNames{6},uniqueNames{7})
lgd.FontSize = 20;
lgd.FontWeight = 'bold';
hold off