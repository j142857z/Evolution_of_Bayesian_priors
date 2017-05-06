#!/bin/bash

./t3 100 100 2 1000 100 1 10

matlab -nodisplay -nosplash -r "OutputSurfacePlot('Output.csv', 'SurfPlot.pdf', 0, 1, 0, 1); exit"
