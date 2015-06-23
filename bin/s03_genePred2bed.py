import sys
f_in = open( sys.argv[1],"r" )
for line in f_in:
   line = line.strip('\n')
   f    = line.split()
   f[5] = ",".join( f[5].split(",")[:-1] )
   f[6] = ",".join( f[6].split(",")[:-1] )
   print "\t".join(f)
f_in.close()