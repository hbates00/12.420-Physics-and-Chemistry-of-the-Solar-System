
%% ---Variable Handling ---

close all

%planet variables
name = catalog{2:end, {'name'}};

mass_p = catalog{2:end, {'mass'}};
mass_p_err_min = catalog{2:end, {'mass_error_min'}};
mass_p_err_max = catalog{2:end, {'mass_error_max'}};

radius_p = catalog{2:end, {'radius'}};
radius_p_err_min = catalog{2:end, {'radius_error_min'}};
radius_p_err_max = catalog{2:end, {'radius_error_max'}};

%semimajor axis variables
axis = catalog{2:end, {'semi_major_axis'}};
axis_err_min = catalog{2:end, {'semi_major_axis_error_min'}};
axis_err_max = catalog{2:end, {'semi_major_axis_error_max'}};

%star variables
mass_s = catalog{2:end, {'star_mass'}};
mass_s_err_min = catalog{2:end, {'star_mass_error_min'}};
mass_s_err_max = catalog{2:end, {'star_mass_error_max'}};

radius_s = catalog{2:end, {'star_radius'}};
radius_s_err_min = catalog{2:end, {'star_radius_error_min'}};
radius_s_err_max = catalog{2:end, {'star_radius_error_max'}};

star_temp = catalog{2:end, {'star_teff'}};
star_temp_err_min = catalog{2:end, {'star_teff_error_min'}};
star_temp_err_max = catalog{2:end, {'star_teff_error_max'}};

for i = linspace(1, length(mass_p), length(mass_p))
    if mass_p_err_min(i) > mass_p(i)
        mass_p_err_min(i) = mass_p(i);
    end
end

for i = linspace(1, length(mass_p), length(mass_p))
    if radius_p_err_min(i) > radius_p(i)
        radius_p_err_min(i) = radius_p(i);
    end
end


%% ---Removing NaNs---

%mass NaNs
indm = [];
for i = linspace(1, length(mass_p), length(mass_p))
    a = isnan(mass_p);
    if a(i) == 1
        indm = [indm, i];
    end
end

indr = [];

%Radius NaNs
for i = linspace(1, length(radius_p), length(radius_p))
    a = isnan(radius_p);
    if a(i) == 1
        indr = [indr, i];
    end
end

ind = [indm, indr];

mass_p(ind) = [];
mass_p_err_min(ind) = [];
mass_p_err_max(ind) = [];

radius_p(ind) = [];
radius_p_err_min(ind) = [];
radius_p_err_max(ind) = [];

mass_p_err_min(isnan(mass_p_err_min)) = 0;
mass_p_err_max(isnan(mass_p_err_max)) = 0;
radius_p_err_min(isnan(radius_p_err_min)) = 0;
radius_p_err_max(isnan(radius_p_err_max)) = 0;

name(ind) = [];

axis(ind) = [];
axis_err_min(ind) = [];
axis_err_max(ind) = [];

%% ---Restricting to Planets with a Mass < 100MEarth---

indp = [];

for i = linspace(1, length(mass_p), length(mass_p))
    if mass_p(i) > 0.315 
        indp = [indp, i];
    end
end

mass_p(indp) = [];
mass_p_err_min(indp) = [];
mass_p_err_max(indp) = [];

radius_p(indp) = [];
radius_p_err_min(indp) = [];
radius_p_err_max(indp) = [];

name(indp) = [];

axis(indp) = [];
axis_err_min(indp) = [];
axis_err_max(indp) = [];

%% ---Earth-Centric Units---

mass_p = 317.828 * mass_p;
mass_p_err_min = 317.828 * mass_p_err_min;
mass_p_err_max = 317.828 * mass_p_err_max;

radius_p = 11.21 * radius_p;
radius_p_err_min = 11.21 * radius_p_err_min;
radius_p_err_max = 11.21 * radius_p_err_max;

%% ---Filtering Out Potential Iron Cores---

bad = [];

k1 = -0.20949;
k2 = 0.0804;
k3 = 0.394;

for i = linspace(1, length(mass_p), length(mass_p))
    testnum = 0;
    
    x1 = mass_p(i) - mass_p_err_min(i);
    x2 = mass_p(i) + mass_p_err_max(i);
    y1 = radius_p(i) - radius_p_err_min(i);
    y2 = radius_p(i) + radius_p_err_max(i);
    
    yt1 = 10.^(k1 + (1./3.0) * log10(x1/5.8) - k2*(x1/5.8).^(k3)) * 2.52;
    yt2 = 10.^(k1 + (1./3.0) * log10(x2/5.8) - k2*(x2/5.8).^(k3)) * 2.52;
    
    if y1 >= yt1
        testnum = testnum + 1;
    else
        testnum = testnum - 1;
    end
    
    if y1 >= yt2
        testnum = testnum + 1;
    else
        testnum = testnum - 1;
    end
    
    if y2 >= yt1
        testnum = testnum + 1;
    else
        testnum = testnum - 1;
    end
    
    if y2 >= yt2
        testnum = testnum + 1;
    else
        testnum = testnum - 1;
    end
    
    if testnum == 4
        bad = [bad, i];
    end
end

mass_p_iron = mass_p;
mass_p_err_min_iron = mass_p_err_min;
mass_p_err_max_iron = mass_p_err_max;

radius_p_iron = radius_p;
radius_p_err_min_iron = radius_p_err_min;
radius_p_err_max_iron = radius_p_err_max;
name_iron = name;

axis_iron = axis;
axis_err_min_iron = axis_err_min;
axis_err_max_iron = axis_err_max;

mass_p_iron(bad) = [];
mass_p_err_min_iron(bad) = [];
mass_p_err_max_iron(bad) = [];

radius_p_iron(bad) = [];
radius_p_err_min_iron(bad) = [];
radius_p_err_max_iron(bad) = [];

name_iron(bad) = [];

axis_iron(bad) = [];
axis_err_min_iron(bad) = [];
axis_err_max_iron(bad) = [];

%% ---Alternative Iron Identification---

bad2 = [];

for i = linspace(1, length(mass_p), length(mass_p))
    x = mass_p(i);
    y = radius_p(i);
    
    t =  10.^(k1 + (1./3.0) * log10(x/5.8) - k2*(x/5.8).^(k3)) * 2.52;
    
    if y >= t;
        bad2 = [bad2, i];
    end
end

mass_p_iron_c = mass_p;
mass_p_err_min_iron_c = mass_p_err_min;
mass_p_err_max_iron_c = mass_p_err_max;

radius_p_iron_c = radius_p;
radius_p_err_min_iron_c = radius_p_err_min;
radius_p_err_max_iron_c = radius_p_err_max;
name_iron_c = name;

axis_iron_c = axis;
axis_err_min_iron_c = axis_err_min;
axis_err_max_iron_c = axis_err_max;

mass_p_iron_c(bad2) = [];
mass_p_err_min_iron_c(bad2) = [];
mass_p_err_max_iron_c(bad2) = [];

radius_p_iron_c(bad2) = [];
radius_p_err_min_iron_c(bad2) = [];
radius_p_err_max_iron_c(bad2) = [];

name_iron_c(bad2) = [];

axis_iron_c(bad2) = [];
axis_err_min_iron_c(bad2) = [];
axis_err_max_iron_c(bad2) = [];
        
%% ---Finding Unreasonable Planets---

byeye = [1, 2, 11];

bad_mass = mass_p_iron_c(byeye);
bad_mass_err_min = mass_p_err_min_iron_c(byeye);
bad_mass_err_max = mass_p_err_max_iron_c(byeye);

bad_radius = radius_p_iron_c(byeye);
bad_radius_err_min = radius_p_err_min_iron_c(byeye);
bad_radius_err_max = radius_p_err_max_iron_c(byeye);


mass_p_iron_c(byeye) = [];
mass_p_err_min_iron_c(byeye) = [];
mass_p_err_max_iron_c(byeye) = [];

radius_p_iron_c(byeye) = [];
radius_p_err_min_iron_c(byeye) = [];
radius_p_err_max_iron_c(byeye) = [];

name_iron_c(byeye) = [];

axis_iron_c(byeye) = [];
axis_err_min_iron_c(byeye) = [];
axis_err_max_iron_c(byeye) = [];

%% ---Extracting Errorbar Groups---

%All planets
full = [];
partial = [];
none = [];

for i = linspace(1, length(mass_p), length(mass_p))
    
    catnum = 0;
    
    if mass_p_err_min(i) == 0 
        catnum = catnum + 1;
    end
    if mass_p_err_max(i) == 0
        catnum = catnum + 1;
    end
    if radius_p_err_max(i) == 0
        catnum = catnum + 1;
    end
    if radius_p_err_min(i) == 0
        catnum = catnum + 1;
    end
    
    if catnum == 4
        none = [none, i];
    elseif catnum == 0
        full = [full, i];
    else
        partial = [partial, i];
    end
end

full_mass_p = mass_p(full);
full_radius_p = radius_p(full);
full_mass_p_err_min = mass_p_err_min(full);
full_mass_p_err_max = mass_p_err_max(full);
full_radius_p_err_min = radius_p_err_min(full);
full_radius_p_err_max = radius_p_err_max(full);

partial_mass_p = mass_p(partial);
partial_radius_p = radius_p(partial);
partial_mass_p_err_min = mass_p_err_min(partial);
partial_mass_p_err_max = mass_p_err_max(partial);
partial_radius_p_err_min = radius_p_err_min(partial);
partial_radius_p_err_max = radius_p_err_max(partial);

nonem = mass_p(none);
noner = radius_p(none);

%Just Iron
full = [];
partial = [];
none = [];

for i = linspace(1, length(mass_p_iron), length(mass_p_iron))
    
    catnum = 0;
    
    if mass_p_err_min_iron(i) == 0 
        catnum = catnum + 1;
    end
    if mass_p_err_max_iron(i) == 0
        catnum = catnum + 1;
    end
    if radius_p_err_max_iron(i) == 0
        catnum = catnum + 1;
    end
    if radius_p_err_min_iron(i) == 0
        catnum = catnum + 1;
    end
    
    if catnum == 4
        none = [none, i];
    elseif catnum == 0
        full = [full, i];
    else
        partial = [partial, i];
    end
end

full_mass_p_iron = mass_p_iron(full);
full_radius_p_iron = radius_p_iron(full);
full_mass_p_err_min_iron = mass_p_err_min_iron(full);
full_mass_p_err_max_iron = mass_p_err_max_iron(full);
full_radius_p_err_min_iron = radius_p_err_min_iron(full);
full_radius_p_err_max_iron = radius_p_err_max_iron(full);

partial_mass_p_iron = mass_p_iron(partial);
partial_radius_p_iron = radius_p_iron(partial);
partial_mass_p_err_min_iron = mass_p_err_min_iron(partial);
partial_mass_p_err_max_iron = mass_p_err_max_iron(partial);
partial_radius_p_err_min_iron = radius_p_err_min_iron(partial);
partial_radius_p_err_max_iron = radius_p_err_max_iron(partial);

nonem_iron = mass_p_iron(none);
noner_iron = radius_p_iron(none);


%Constrained Iron
full2 = [];
partial2 = [];
none2 = [];

for i = linspace(1, length(mass_p_iron_c), length(mass_p_iron_c))
    
    catnum = 0;
    
    if mass_p_err_min_iron_c(i) == 0 
        catnum = catnum + 1;
    end
    if mass_p_err_max_iron_c(i) == 0
        catnum = catnum + 1;
    end
    if radius_p_err_max_iron_c(i) == 0
        catnum = catnum + 1;
    end
    if radius_p_err_min_iron_c(i) == 0
        catnum = catnum + 1;
    end
    
    if catnum == 4
        none2 = [none2, i];
    elseif catnum == 0
        full2 = [full2, i];
    else
        partial2 = [partial2, i];
    end
end

full_mass_p_iron_c = mass_p_iron_c(full2);
full_radius_p_iron_c = radius_p_iron_c(full2);
full_mass_p_err_min_iron_c = mass_p_err_min_iron_c(full2);
full_mass_p_err_max_iron_c = mass_p_err_max_iron_c(full2);
full_radius_p_err_min_iron_c = radius_p_err_min_iron_c(full2);
full_radius_p_err_max_iron_c = radius_p_err_max_iron_c(full2);

partial_mass_p_iron_c = mass_p_iron_c(partial2);
partial_radius_p_iron_c = radius_p_iron_c(partial2);
partial_mass_p_err_min_iron_c = mass_p_err_min_iron_c(partial2);
partial_mass_p_err_max_iron_c = mass_p_err_max_iron_c(partial2);
partial_radius_p_err_min_iron_c = radius_p_err_min_iron_c(partial2);
partial_radius_p_err_max_iron_c = radius_p_err_max_iron_c(partial2);

nonem_iron_c = mass_p_iron_c(none2);
noner_iron_c = radius_p_iron_c(none2);



%% ---All Planets---
hold on

%errorbar(full_mass_p, full_radius_p, full_radius_p_err_min, full_radius_p_err_max, full_mass_p_err_min, full_mass_p_err_max, 'o', 'Color', 'b', 'MarkerFaceColor', 'b')

%errorbar(partial_mass_p, partial_radius_p, partial_radius_p_err_min, partial_radius_p_err_max, partial_mass_p_err_min, partial_mass_p_err_max, 's', 'Color', 'b', 'MarkerFaceColor', 'b')

%plot(nonem, noner, 'o', 'Color', 'b')


%% ---Just Overdense Planets---
redfade = [204, 204, 204]/255.0;

%errorbar(full_mass_p_iron, full_radius_p_iron, full_radius_p_err_min_iron, full_radius_p_err_max_iron, full_mass_p_err_min_iron, full_mass_p_err_max_iron, 'o', 'Color', 'r', 'MarkerFaceColor', 'r')

%errorbar(partial_mass_p_iron, partial_radius_p_iron, partial_radius_p_err_min_iron, partial_radius_p_err_max_iron, partial_mass_p_err_min_iron, partial_mass_p_err_max_iron, 's', 'Color', 'r', 'MarkerFaceColor', 'r')

%plot(nonem_iron, noner_iron, 'o', 'Color', 'r')

%errorbar(bad_mass, bad_radius, bad_radius_err_min, bad_radius_err_max, bad_mass_err_min, bad_mass_err_max, 'o', 'Color', 'r', 'MarkerFaceColor', 'r')
%{
errorbar(full_mass_p_iron_c, full_radius_p_iron_c, full_radius_p_err_min_iron_c, full_radius_p_err_max_iron_c, full_mass_p_err_min_iron_c, full_mass_p_err_max_iron_c, 'o', 'Color', 'r', 'MarkerFaceColor', 'r')

errorbar(partial_mass_p_iron_c, partial_radius_p_iron_c, partial_radius_p_err_min_iron_c, partial_radius_p_err_max_iron_c, partial_mass_p_err_min_iron_c, partial_mass_p_err_max_iron_c, 's', 'Color', 'r', 'MarkerFaceColor', 'r')

plot(nonem_iron_c, noner_iron_c, 'o', 'Color', 'r')


%% ---Plotting Mass/Radius Relationship---

M = linspace(0.1, 1500, 1000);
R = 10.^(k1 + (1./3.0) * log10(M/5.8) - k2*(M/5.8).^(k3)) * 2.52;

semilogx(M, R, 'Color', 'k', 'LineWidth', 3)

%% ---Plotting Names of Superdense Planets---

dx = 0.05;
dy = 0.05;

text(mass_p_iron_c+dx, radius_p_iron_c+dy, name_iron_c)
%}
%% ----Plotting Semimajor Axis---

c = categorical({'Kepler-102 b', 'Kepler-102 c', 'Kepler-106 d', 'Kepler-114 b','Kepler-114 c', 'Kepler-128 b', 'Kepler-128 c', 'Kepler-131 c', 'Kepler-20 e', 'Kepler-20 f', 'Kepler-37 b', 'Kepler-37 c', 'Kepler-406 c', 'Kepler-408 b', 'Kepler-409 b', 'Kepler-42 b', 'Kepler-42 c', 'Kepler-42 d', 'Kepler-62 b', 'Kepler-62 c', 'Kepler-62 e', 'Kepler-62 f', 'Kepler-68 c'});
bar(c, axis_iron_c(1:23), 'w')

a = categorical({'Kepler-102 b', 'Kepler-102 c'});
A = axis_iron_c(1:2);
bar(a, A)
b = categorical({'Kepler-114 b','Kepler-114 c'});
B = axis_iron_c(4:5);
bar(b, B)
c = categorical({'Kepler-128 b', 'Kepler-128 c'});
C = axis_iron_c(6:7);
bar(c,C)
d = categorical({'Kepler-20 e', 'Kepler-20 f'});
D = axis_iron_c(9:10);
bar(d,D)
e = categorical({'Kepler-37 b', 'Kepler-37 c'});
E = axis_iron_c(11:12);
bar(e,E)
f = categorical({'Kepler-42 b', 'Kepler-42 c', 'Kepler-42 d'});
F = axis_iron_c(16:18);
bar(f, F)
g = categorical({'Kepler-62 b', 'Kepler-62 c', 'Kepler-62 e', 'Kepler-62 f'});
G = axis_iron_c(19:22);
bar(g, G)
h = categorical({'Kepler-70 b'});
H = axis_iron_c(24);
bar(h, H)

plot(xlim,[0.2 0.2])
text(.2, .2, "Approximate Expected Upper Bound of Hot Jupiter Semimajor Axis - 0.2au")

title('Semimajor Axis for High-Interest Superdense Exoplanets')
xlabel('Exoplanet Name')
ylabel('Semimajor Axis (au)')


%% ---General Plotting---%

%title('Mass and radius relationship of exoplanets')
%xlabel('Mass (MEarth)')
%ylabel('Radius (REarth)')
%xlim([0.3, 100])
%ylim([0.25, 2.5])

%0.3
%0.25

%line([0.55, 0.55], [0, 100],'Color', 'b', 'LineStyle', '--')
%line([45,45], [0, 100], 'Color', 'm', 'LineStyle', '--')
%line([22,22], [0, 100], 'Color', 'k', 'LineStyle', '--')
%{
text(1, 1, 'Jupiter')
text(1, 1, 'Saturn')
text(1, 1, 'Uranus')

%legend('Fully constrained non-anomalous density planets', 'Partially constrained non-anomalous density planets', 'Unconstrained constrained non-anomalous density planets', 'Fully constrained potentially superdese Planets', 'Partially Constrained Potentially Superdense Planets', 'Unconstrained Potentially Superdense Planets', 'Mass-Radius Relationship for Homogeneous Iron Sphere') 
%legend('Fully Constrained Low-Interest Planets', 'Partially Constrained Low-Interest Planets', 'Unconstrained Low-Interest Planets', 'Fully Constrained High-Interest Planets', 'Partially Constrained High-Interest Planets', 'Unconstrained High-Interest Planets', 'Mass-Radius Relationship for Homogeneous Iron Sphere') 
legend('Fully Constrained High-Interest Planets', 'Partially Constrained High-Interest Planets', 'Unconstrained High-Interest Planets', 'Mass-Radius Relationship for Homogeneous Iron Sphere')%, "Approximate Location of Uranus's Rocky Core", "Approximate Location of Jupiter's Rocky Core", "Approximate Location of Saturn's Rocky Core") 


set(gca,'xscale','log')
%}
