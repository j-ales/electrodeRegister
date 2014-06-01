function [] = setPlotOptions()
%Setup some default plotting options for a nice looking plot

material dull;
lightangle(240, 30)
lightangle(120, 30)
lightangle(0, 0)
axis normal
axis equal
axis vis3d
axis tight;
axis off
campos(1000*[    1    2.5    1]);
camva(4.5);