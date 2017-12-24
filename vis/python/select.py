

file = open("list", "r")

lines=file.readlines()

file.close()

size=len(lines)

partlines=[]

for i in range(0,size,15):
  partlines.append(lines[i])

outfile=open("listnew","w")

for line in partlines:
   outfile.writelines(line)

outfile.close()
