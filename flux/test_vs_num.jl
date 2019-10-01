using UnicodePlots
t=1
plt = heatmap(reshape(P[:,1],50,50), xscale=0.01,
             yscale=0.01, xoffset=0, colormap=:inferno, height = 50, width = 50);   display(plt);

pt1 = reshape(P[:,1],50,50)
plt = lineplot(1:50, pt1[:,10], title = "Num vs. Nero", name = "10", xlabel = "x", ylabel = "P",
                width = 50,xlim = [0, 50],ylim = [80, 90]);
lineplot!(plt,1:50,pt1[:,25], name = "25");
display(plt)


funkh, dkx, dky = funKH(xa,xb,ya,yb)


kh = map(x->funkh(x[1]/1000,x[2]/1000),xy)
plt = heatmap(kh, xscale=0.01,
             yscale=0.01, xoffset=0, colormap=:inferno, height = 50, width = 50);   display(plt);
kh = kh
kh = kh.+0.1
