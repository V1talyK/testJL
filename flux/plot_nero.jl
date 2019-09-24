nc = Int64(floor(sqrt(length(xy))));
NeS = (x->Tracker.data(m3(x))[1]).(xy);     NeS = reshape(NeS,nc,nc);
AnS = one_well.(xy,rw);    AnS = reshape(Float64.(AnS),nc,nc)



plt = heatmap(BoS, xscale=0.1, yscale=0.1, xoffset=0, colormap=:inferno);   display(plt);

plt = heatmap(AnS-NeS, xscale=0.01,
                   yscale=0.01, xoffset=0, colormap=:inferno, height = 21, width = 21);   display(plt);
plt = heatmap(AnS, xscale=0.01,
                   yscale=0.01, xoffset=0, colormap=:inferno, height = 21, width = 21);   display(plt);
plt = heatmap(NeS, xscale=0.01,
                   yscale=0.01, xoffset=0, colormap=:inferno, height = 21, width = 21);    display(plt);

plt = heatmap(PyS, xscale=0.1, yscale=0.1, xoffset=0, colormap=:inferno);     display(plt);

plt = heatmap(clamp.(NeS,0,1).-AnS, xscale=0.01, yscale=0.01, xoffset=0, colormap=:inferno); display(plt)

NeS = (x->Tracker.data(m3(x))[1]).(xy);     NeS = reshape(NeS,nc,nc);
BoS = AxyW.(xy);    BoS = reshape(BoS,nc,nc)

x_ax = 1:nc;
nh = Int64.(floor(nc/2)+1)
P = AnS[nh,:]
plt = lineplot(x_ax, NeS[nh,:], title = "Press", name = "Ner", xlabel = "x", ylabel = "P",width = 21,xlim = [0, 21]);
#lineplot!(plt,x_ax,NeS[nh,:], name = "NeS");
lineplot!(plt,x_ax,P, name = "Analitic");
#lineplot!(plt,x_ax,BoS[nh,:],name = "AxyW");
#lineplot!(plt,x_ax,CoS[nh,:],name = "Cos");
display(plt)

plt = scatterplot(NeS[nh,:],P, xlabel = "Analitic", ylabel = "Nero" )

plt = heatmap(reshape(kh,21,21), xscale=0.01,
                   yscale=0.01, xoffset=0, colormap=:inferno, height = 21, width = 21);   display(plt);
