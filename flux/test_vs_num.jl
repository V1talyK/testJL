using UnicodePlots
t=1
plt = heatmap(reshape(P[:,1],50,50), xscale=0.01,
             yscale=0.01, xoffset=0, colormap=:inferno, height = 50, width = 50);   display(plt);

pt1 = reshape(P[:,1],50,50)
plt = lineplot(1:50, pt1[:,10], title = "Num vs. Nero", name = "10", xlabel = "x", ylabel = "P",
                width = 50,xlim = [0, 50],ylim = [80, 100]);
lineplot!(plt,1:50,pt1[:,25], name = "25");
display(plt)


funkh, dkx, dky = funKH(xa,xb,ya,yb)


kh = map(x->funkh(x[1]/1000,x[2]/1000),xy)
kh = kh
kh = kh.+0.1

plt = heatmap(reshape(kh,50,50), xscale=0.01,
             yscale=0.01, xoffset=0, colormap=:inferno, height = 50, width = 50);   display(plt);
Base.source_path()

using JLD
JLD.save(joinpath(rdata,"PPP1.jld"),"P",P)
Pnum = JLD.load("/home/lik/Документы//testJL/flux/PPP1.jld")
Pnum = Pnum["P"]

mean(abs.(NeS[:]*100 .-Pnum[:,1]))
nh = 25
x_ax = 10:20:1000
plt = lineplot(x_ax, NeS[:,nh], title = "Press", name = "Ner", xlabel = "x", ylabel = "P",width = 50,xlim = [0, 50]);
lineplot!(plt,x_ax,reshape(Pnum[:,1],50,50)[:,nh], name = "Num");
display(plt)

z0 = (x->fRBF(x')[1]).(xy)
z1 = funB0.(xy)
zp = (x->prod(pw .-R(x,wxy,pw,rw))).(xy)
zp = (x->prod(1 .-fromWell(x,wxy,rw))).(xy)
zpp = (x->Tracker.data(pde_trialA(x,m1(x)))[1]).(xy)

plt = heatmap(reshape(zp.*z1,50,50), xscale=0.01,
             yscale=0.01, xoffset=0, colormap=:inferno, height = 50, width = 50);   display(plt);
