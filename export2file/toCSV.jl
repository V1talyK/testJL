outfile = string(rsrc,"/qw.txt")
open(outfile, "w") do f
  for i in qw
      for j in i
          print(f, round(j,digits = 3), " ")
      end
      println(f,"")
  end
end # the fi
