%%

color1 = [0.7490 0.7490 0.7490];
color2 = [0 0 0];
color3 = [0.2745 0.1843 1.0000];
color4 = [1.0000 0.2980 0.2980];

%%
N = 3.627288819244943559385972548678215954267142566959751603261 * n;

hold on

plot( log10(Dkm), log10(N), 'Color', color2);
title('Number of Craters per million km^2')
ylabel('Log of Number per million km^2')
xlabel('Log of Diameter (km)')


%% 

slice = log10(Dkm) < 1.1;

dcut = log10(Dkm(slice));
ncut = log10(N(slice));

plot(dcut, ncut, 'Color', color3)

%%

slice2 = log10(Dkm) > 1.6;

dcut2 = log10(Dkm(slice2));
ncut2 = log10(N(slice2));

plot(dcut2, ncut2, 'Color', color4)

%%

func = polyfit(dcut, ncut, 1)

fit = polyval(func, dcut);

plot(dcut, fit, 'Color', color1, 'linestyle', '--')

%%

func2 = polyfit(dcut2, ncut2, 1)

fit2 = polyval(func2, dcut2);

plot(dcut2, fit2, 'Color', color1, 'linestyle', '--')

legend( 'All D', 'D ~ 0 - 1.1', 'D ~ 1.6 - 2.0', 'Fit1', 'Fit2')

%%
ans1 = (polyval(func, log10(5)))
ans2 = (polyval(func, log10(16)))

ans3 = (polyval(func2, log10(5)))
ans4 = (polyval(func2, log10(16)))

scatter(log10(5), ans1, 'MarkerEdgeColor', color3)
scatter(log10(16), ans2, 'MarkerEdgeColor', color3)

scatter(log10(5), ans3, 'MarkerEdgeColor', color4)
scatter(log10(16), ans4, 'MarkerEdgeColor', color4)


