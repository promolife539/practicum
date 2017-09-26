N = 8
listN = []
for i in range(0,8) :
    listN.append(0)

count=0
listk = [2,5,7,1,8,5,4,5,3,6,3,4,1]
for i in range(0, len(listk)) :
    listN[listk[i]-1]+=1
for i in range(0,8) :
    if listN[i] != 0 :
        count += 1

print(count)