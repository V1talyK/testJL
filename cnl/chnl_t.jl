x = zeros(length(b),8)
ch = Channel{Any}(4)
ts = Threads.@spawn begin
    thr_id = Threads.threadid()
    for item in ch
        CL=cholesky(mA);
        x[:,thr_id] = CL\items;
        println(extrema(view(x,:,thr_id)))
    end
end


put!(ch,b)
for i=1:100
    put!(ch,rand(length(b)))
end

x = zeros(length(b),8)
CL8=[cholesky(mA) for i=1:8];
CL=cholesky(mA)

function foo(x,CL8)
    @sync for i=1:100
        Threads.@spawn begin
            thr_id = Threads.threadid()
            @views cholesky!(CL8[thr_id],mA)
            @views x[:,thr_id] = CL8[thr_id]\b;
            #x[:,thr_id] = CL\b;
        end
    end
end

@time foo(x,CL8)
@btime foo($x,$CL8)


ch1 = Channel{Any}(3)
k=0
f2() = f1(ch1)
ts = @task begin
    for item in ch1
        thr_id = Threads.threadid()
        println(thr_id)
        println(item)
        k+=1
    end
end

ts = @task f2()
schedule(ts)
yield()
put!(ch1,2)
istaskstarted(ts)
istaskdone(ts)
