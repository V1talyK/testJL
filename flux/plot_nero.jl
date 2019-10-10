nc = Int64(floor(sqrt(length(xy))));
NeS = (x->Tracker.data(m3(x))[1]).(xy);     NeS = reshape(NeS,nc,nc);
AnS = one_well.(xy,rw);    AnS = reshape(Float64.(AnS),nc,nc)

pl[isnan.(pl)].=0

plt = heatmap(reshape(pl,nc,nc), xscale=0.1, yscale=0.1, xoffset=0, colormap=:inferno);   display(plt);

plt = heatmap(AnS-NeS, xscale=0.01,
                   yscale=0.01, xoffset=0, colormap=:inferno, height = 21, width = 21);   display(plt);
plt = heatmap(AnS, xscale=0.01,
                   yscale=0.01, xoffset=0, colormap=:inferno, height = 21, width = 21);   display(plt);
plt = heatmap(NeS, xscale=0.01,
                   yscale=0.01, xoffset=0, colormap=:inferno, height = 50, width = 50);    display(plt);

plt = heatmap(PyS, xscale=0.1, yscale=0.1, xoffset=0, colormap=:inferno);     display(plt);

plt = heatmap(clamp.(NeS,0,1).-AnS, xscale=0.01, yscale=0.01, xoffset=0, colormap=:inferno); display(plt)

BoS = AxyW.(xy);    BoS = reshape(BoS,nc,nc)

x_ax = 1:nc;
nh = Int64.(floor(nc/2)+1)
P = AnS[nh,:]
plt = lineplot(x_ax, NeS[:,nh], title = "Press", name = "Ner", xlabel = "x", ylabel = "P",width = 50,xlim = [0, 50]);
#lineplot!(plt,x_ax,NeS[nh,:], name = "NeS");
lineplot!(plt,x_ax,P, name = "Analitic");
#lineplot!(plt,x_ax,BoS[nh,:],name = "AxyW");
#lineplot!(plt,x_ax,CoS[nh,:],name = "Cos");
display(plt)

plt = scatterplot(NeS[nh,:],P, xlabel = "Analitic", ylabel = "Nero" )

plt = heatmap(reshape(Pnum[:,1],50,50), xscale=0.01,
                   yscale=0.01, xoffset=0, colormap=:inferno, height = 50, width = 50);   display(plt);

plt = heatmap(reshape(Pnum[:,1],50,50).-clamp.(NeS,60,100), xscale=0.01,
                  yscale=0.01, xoffset=0, colormap=:inferno, height = 50, width = 50);   display(plt);

plt = lineplot(0:0.01:1, fsip.(0:0.01:1));
display(plt)

plt = lineplot(1:50, reshape(kh,50,50)[10,:], title = "Press", name = "Ner", xlabel = "x", ylabel = "P",width = 50,xlim = [0, 50]);
