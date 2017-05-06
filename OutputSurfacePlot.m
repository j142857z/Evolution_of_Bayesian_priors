function OutputSurfacePlot(inputFileName, outputFileName, xMin, xMax, yMin, yMax)
GenerateSurfacePlot(inputFileName, xMin, xMax, yMin, yMax)
print('-dpdf',outputFileName)
end

%  function OutputSurfacePlot()
%  GenerateSurfacePlot('fileOfAvgBayesianPopulationByHyperparameterPair.csv')
%  print -dpdf SurfPlot.pdf -r100
%  end

function GenerateSurfacePlot(filename, xMin, xMax, yMin, yMax)

T = csvread(filename)

x = T(:, [1])
y = T(:, [2])
z = T(:, [3])

[qx,qy] = meshgrid(linspace(min(x),max(x)),linspace(min(y),max(y)));
F = TriScatteredInterp(x,y,z);
qz = F(qx,qy);

set(gcf,'renderer','painters');

set(gca, 'FontSize', 20);

h=surf(qx,qy,qz,'EdgeColor','none')

% Create xlabel
%  xlabel('$\alpha$','interpreter','latex', 'FontSize', 20)
xlabel('\alpha_S', 'FontSize', 20)

% Create ylabel
%  ylabel('$\beta$','interpreter','latex', 'FontSize', 20)
ylabel('\beta_S', 'FontSize', 20)

% Create zlabel
zlabel('Percentage of Bayesians in population', 'FontSize', 20);

h=colorbar;
set(h,'fontsize', 20);

axis([xMin, xMax, yMin, yMax])

view([90 90])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 20]);

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [25 20]);

view([0 90])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 20]);

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [25 20]);
end