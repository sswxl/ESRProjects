filename='2007_train.txt'
f=open(filename,'r')
allfile=f.readlines()

transfile=[]
for i in range(len(allfile)):
    transfile.append(allfile[i].replace(r"C:\Users\WHY\Desktop\yolo3-pytorch-master/",""))
print(transfile)
f.close()

f=open(filename,'w')
for i in range(len(transfile)):
    f.write(transfile[i])
f.close()