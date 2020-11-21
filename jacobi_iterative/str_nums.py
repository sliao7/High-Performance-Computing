from decimal import Decimal
filename = "5000A16A1A.txt"
f = open(filename,"r")
average = 0
for i in range(10):
	average += float(f.readline().strip())
average /= 10
print(average)