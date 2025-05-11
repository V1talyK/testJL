using LinearAlgebra
A = rand(100, 100)
A = A+A'
b = rand(100)
M = rand(4,100)
x = A\b
y0 = M*x

z = A\M'
y1 = z'*b

y1-y0


using Mmap
io = open("/media/lik/Том N/1/mmap.bin", "w+");
io2 = open("/pool0/media/mmap.bin", "w+");
B = mmap(io, Array{Float64,2}, (25000,30000));
D = mmap(io2, Array{Float64,2}, (25000,30000));
C = rand(25000, 30000);

B[3, 4] = 1;
Mmap.sync!(B);
close(io);
open("mmap.bin", "r+");
C = mmap(io, BitArray, (25,30000));
close(io)
close(io2)
rm("/media/lik/Том N/1/mmap.bin")
rm("/pool0/media/mmap.bin")

@time B .= B.+B;
@time C .= C.+C;
@time D .= D.+D;

@time Threads.@threads for i = axes(B)[2]
    C[:,i] .= view(B,:,i);
end