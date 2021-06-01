 using LoopVectorization

function mydot(a, b)
      s = 0.0
      @inbounds @simd for i ∈ eachindex(a,b)
          s += a[i]*b[i]
      end
      s
 end

function mydotavx(a, b)
      s = 0.0
      @turbo for i ∈ eachindex(a,b)
          s += a[i]*b[i]
      end
      s
end


a = rand(1024*10); b = rand(1024*10); c = similar(a)
b = view(b,rand(1:length(b),length(b)))

@btime mydot($a, $b)
@btime mydotavx($a, $b)
@btime dot($a, $b)

function mycopy(a,b)
    @turbo thread = 3 for i ∈ eachindex(a,b)
        a[i] = b[i]
    end
end

@btime mycopy($a,$b)
@btime copyto!($a,$b)

function vmul!(c,a,b)
    @turbo thread = 4 for i ∈ eachindex(c,a,b)
        c[i] = a[i]*b[i]
    end
end

@btime vmul!($c,$a,$b)
@btime $c .= $a.*$b;
