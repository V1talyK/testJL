Ra = 0:0.1:1;
Rb = 0:0.1:10;
R1 = 1;
R2 = 1;
Rv = zeros(length(Ra)*length(Rb))
k=0
for v in Iterators.product(Ra,Rb)
    global k+=1
    #println(v)
    Rv[k] = 2*(R1+R2+v[1])*(v[2]+v[1])/(R1+R2+v[1]+v[2]+v[1])
end
Rm = reshape(Rv,length(Ra),length(Rb))
plt = heatmap(Rm, xscale=0.1,
                  yscale=0.1, xoffset=0, colormap=:inferno);   display(plt);
