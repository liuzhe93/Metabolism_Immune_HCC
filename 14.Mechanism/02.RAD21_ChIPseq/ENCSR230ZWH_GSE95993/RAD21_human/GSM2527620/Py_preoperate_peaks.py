import os,sys

for i in range(len(sys.argv)):
    "sys.argv[%d] = %s" % (i, sys.argv[i])

if len(sys.argv) != 3:
    print('''python Py.py input output''')
    exit(0)

fin = open(sys.argv[1], "r")
fout = open(sys.argv[2], "w")

for line in fin:
    line = line.strip()
    if(line.startswith("chr")):
        fout.write(line + "\n")
    else:
        continue

fin.close()
fout.close()

