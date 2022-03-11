
ch = Channel{Any}(4)
@async begin
    for item in ch
        x = mA\item
        println(extrema(x))
    end
end


put!(ch,b)
for i=1:100
    put!(ch,rand(length(b)))
end
