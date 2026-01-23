using Plots
XX = zeros(32, 16)
YY = zeros(32)

for (k,v) in enumerate(17:32)
    if k==1
        XX[v,k] = -1
        XX[v,k+1] = 1
    elseif k==16
        XX[v,k] = -1
        XX[v,k-1] = 1
    else
        XX[v,k] = -2
        XX[v,k+1] = 1
        XX[v,k-1] = 1
    end
end

for (k,v) in enumerate(17:32)
    if k==1
        XX[v,k] = -2
        XX[v,k+1] = 1.5
        XX[v,k+2] = 0.5
    elseif k==2
        XX[v,k-1] = 1.5
        XX[v,k] = -3.5
        XX[v,k+1] = 1.5
        XX[v,k+2] = 0.5
    elseif k==15
        XX[v,k-2] = 0.5
        XX[v,k-1] = 1.5
        XX[v,k] = -3.5
        XX[v,k+1] = 1.5
    elseif k==16
        XX[v,k-2] = 0.5
        XX[v,k-1] = 1.5
        XX[v,k] = -2
    else
        XX[v,k-2] = 0.5
        XX[v,k-1] = 1.5
        XX[v,k] = -4
        XX[v,k+1] = 1.5
        XX[v,k+2] = 0.5

    end
end


XX[1:8,4] .= 1
XX[9:16,11] .= 1

YY[1:8] .= 1
YY[9:16] .= 2

XX[20,:].=0
XX[27,:].=0
XX[23,7] = -1.1
XX[23,8] = 0.1

XX[24,7] = 0.1
XX[24,8] = -1.1

aa = (XX'*XX)\(XX'*YY)

plot(aa, lw = 2)
1+1