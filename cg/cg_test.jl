using LinearAlgebra
a = rand(3)
norm(a)
x = zeros(2,10)
df = zeros(2,10)
w = zeros(10)
S = similar(x)

function fun(x)
    2*x[1]^2 + 2*x[2]^2 + 2*x[1]*x[2] + 20*x[1] + 10*x[2] + 10
end

function dfun(x)
    [4*x[1]+2*x[2]+20,2*x[1]+4*x[2]+10]
end

function make_dt(x,df)
    chi = 4*df[1]*x[1]+4*df[2]*x[2]+2*x[1]*df[2]+2*df[1]*x[2]+20*df[1]+10*df[2]
    zna = 4*df[1]^2+4*df[2]^2+4*df[1]*df[2]
    chi/zna
end

ε=0.1;
d = zeros(2,10)

f = fun(x[:,1])
df[:,1] = dfun(x[:,1])
norm(df)<ε
t = make_dt(x[:,1],df[:,1])
x[:,2] = x[:,1]-t.*df[:,1]

for i=1:9
    #f = fun(x[:,i])
    df[:,i] = dfun(x[:,i])
    norm(df[:,i])<ε
    if i>1
        b = (norm(df[:,i])/norm(df[:,i-1]))^2
        d[:,i] = df[:,i]+b.*df[:,i-1]
        t = make_dt(x[:,i],d[:,i])
        x[:,i+1] = x[:,i]-t.*d[:,i]
    else
        t = make_dt(x[:,i],df[:,i])
        x[:,i+1] = x[:,i]-t.*df[:,i]
    end

end

A = [-2 1 1 0; 1 -2 0 1; 1 0 -2 1; 0 1 1 -2];
W = ones(4)
W1 =  diagm(0=>W)
AW = A-diagm(0=>W)
pwf = [1 0 0 0.8]'

function sim(pw)
    b = -W.*pw
    p = AW\b
    q = vec(W.*(p-pw))
end

function fun_dJ_dpw(qf, qc)
    dp_dpw = -AW\W1
    dq_dpw = W1*dp_dpw-W1
    dJ_dpw = -2*vec((qf-qc)'*dq_dpw)
end

function fun_dJ_dt(dpw,qf, qc)
    dp_dt = -AW\(W1*dpw)
    dq_dt = W1*dp_dt+W1*dpw
    dJ_dt = -2*(qf-qc)'*dq_dt
end

function make_dt(qf,qc,W1,AW,dpw)
    dp_dt = -AW\(W1*dpw)
    chi=sum((qf.-qc).*W1*(dp_dt+dpw))
    zna=sum(W1.*dpw.*W1*(dp_dt+dpw))
    chi/zna
end

pwc = [1 0 0 0.5]'

qf = sim(pwf)
qc = sim(pwc)

dpw = fun_dJ_dpw(qf, qc)

dpw = zeros(4,10)
d = zeros(4,10)
x = zeros(4,10)
for i=1:9
    #f = fun(x[:,i])
    qc = sim(x[:,i])
    dpw[:,i] = fun_dJ_dpw(qf, qc)
    norm(dpw[:,i])<ε
    if i>1
        b = (norm(dpw[:,i])/norm(dpw[:,i-1]))^2
        d[:,i] = dpw[:,i]+b.*dpw[:,i-1]
        t = make_dt(qf,qc,W1,AW,d[:,i])
        x[:,i+1] = x[:,i]-t.*d[:,i]
    else
        t = make_dt(qf,qc,W1,AW,dpw[:,i])
        x[:,i+1] = x[:,i]-t.*dpw[:,i]
    end

end
