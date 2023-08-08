using UnicodePlots, Term

panel(p; kw...) = begin
  p.margin[] = p.padding[] = 0  # make plots more compact
  Panel(string(p; color=true); style="hidden", fit=true, kw...)
end

function plot_P_lr(PM, P_lr, nw)
    plt = Vector(undef, Int64(ceil(nw)))
    for (k,v) in enumerate(1:nw)

        plt[v] = lineplot(PM[v,:], ylim = [floor(minimum(PM)),ceil(maximum(PM))],
                name = "fact", ylabel = "P", title = "скв. $(v)")
        #for i in v[2:end]
            lineplot!(plt[k], P_lr[v,:], name = "calc")
        #end
    end
    grid(panel.(plt); layout=(3, nothing)) |> print
end


function plot_P_lr_ex(PM, P_lr, P_ex, nw)
    plt = Vector(undef, Int64(ceil(nw)))
    for (k,v) in enumerate(1:nw)
        x_ex = size(PM,2) - size(P_ex, 2) +1 : size(PM,2)
        plt[v] = lineplot(PM[v,:], ylim = [floor(minimum(PM)),ceil(maximum(PM))],
                name = "fact", ylabel = "P", title = "скв. $(v)")
        #for i in v[2:end]
            lineplot!(plt[k], P_lr[v,:], name = "calc")
            lineplot!(plt[k], x_ex, P_ex[v,:], name = "exam")
        #end
    end
    grid(panel.(plt); layout=(3, nothing)) |> print
end
