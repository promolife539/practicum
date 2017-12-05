import random
N = 8
listN = []
for i in range(0,8) :
    listN.append(0)

count=0

for i in range(0,8):
    listk[i] = random.randint(1,8)
print(listk)
for i in range(0, len(listk)) :
    listN[listk[i]-1]+=1
for i in range(0,8) :
    if listN[i] != 0 :
        count += 1

print(count)