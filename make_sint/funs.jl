function writeToFile(name,D)
    r=dirname(Base.source_path())
    rfile = joinpath(r,"$name.tsv"); #Синтетика 2
    OUT = joinpath(rfile)
    ioW = Base.open(OUT,"w");
    writeToPipe(ioW, D...)
    close(ioW)
end

function writeToPipe(io, args...)
    k = 0

    newline = '\n'
    delim = '\t'
    @inbounds while k < length(args[1])
        k += 1
        ncols = length(args)
        @inbounds for i = 1:ncols
            v = getStr(args[i][k])
            write(io, "$v")
            write(io, ifelse(i == ncols, newline, delim))
        end

    end
    close(io)
    return true
end


@inline getStr(v::Number) = ifelse(isnan(v), "\\N", string(v))
#@inline getVal(v::Array{T,1}) where T <: Number = (length(v) > 0) ? "["*join(map(x -> string(x), v),", ")*"]" : "\\N"
#@inline getVal(v::Array{T,1}) where T <: AbstractString = (length(v) > 0) ? "["*join(map(x -> "'$x'", v),", ")*"]" : "\\N"
@inline getStr(v::Tuple) = string(v)
@inline getStr(v::Date) = "$v"
@inline getStr(v::String) = "$v"
@inline getStr(v::Array{Float64,1}) = "$v"
@inline getStr(v::Array{<:Signed,1}) = "$v"
#@inline getVal(v::Nothing) = "\\N"
function mlin(ts)
    return collect((ts.-ts[1])./(ts[end]-ts[1]))
end
