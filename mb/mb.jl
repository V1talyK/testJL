function mb(q, Pa; lam = 1)
    nt = length(q)
    P = zeros(length(q))
    dP_dx = zeros(length(q))
    dP_du = zeros(length(q))
    dP_db = zeros(length(q))
    dP_dl = zeros(length(q))
    bet = 150
    dt = 30
    P0 = Pa;

    function mbt(x,u)
        x = (lam*Pa+bet/dt*x - u)/(bet/dt+lam)
        return x
    end

    for t=1:nt
        P[t] = mbt(P0, q[t])
        dP_dx[t] = bet/dt/(bet/dt+lam)
        dP_du[t] = 1 ./ (bet/dt+lam)
        dP_db[t] = ((1/dt*P0)*(bet/dt+lam) - (lam*Pa+bet/dt*P0 - q[t])*(1/dt))/(bet/dt+lam)^2
        dP_dl[t] = ((Pa)*(bet/dt+lam) - (lam*Pa+bet/dt*P0 - q[t]))/(bet/dt+lam)^2
        P0 = P[t]
    end
    return P, dP_dx, dP_du, dP_db, dP_dl, mbt
end
