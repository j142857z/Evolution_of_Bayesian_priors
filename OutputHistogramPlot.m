function OutputHistogramPlot(inputFileName, outputFileName)
GenerateHistogramPlot(inputFileName)
%  print -dpdf HistPlot.pdf -r100
print('-dpdf',outputFileName)
end

function GenerateHistogramPlot(filename)

T = csvread(filename)

x = T(:, [1])

xcenters = 1:100

set(gcf,'renderer','painters');

set(gca, 'FontSize', 20);

hist(x(1:50000), xcenters)

% Create xlabel
xlabel('Bayesian population', 'FontSize', 20)

% Create ylabel
ylabel('Frequency', 'FontSize', 20)

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