f = open("../.env")
f_src = open("../src/batch/batch.py", 'a')
const_list = f.read().split()
for c in const_list:
    par = c.split(sep="=")[0]
    # print(par, type(par))
    f_src.write(par+" = float(os.environ.get(\""+par+"\"))\n")
# print(f.read().split())
f.close()
f_src.close()
