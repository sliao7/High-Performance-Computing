#generate pbs files to be used by submit.sh to run in the terminal
ns = [100, 500, 1000, 5000, 10000]
ps = [1, 4, 9, 16, 25, 36, 49, 64]
nodes = [1,1,1,1,1,2,2,3]
ppn = [1,4,9,16,25,18,25,22]
diff = [0.1, 0.5, 0.9]


for n in ns:
	for i,p in enumerate(ps):
		for d in diff:
				instance = str(n)+"A"+str(p)+"A"+str(int(10*d))+"A"
				fileName = instance +".pbs"
				f = open(fileName,"w")
				s = '''#PBS -N PA3
#PBS -l nodes='''+str(nodes[i])+''':ppn='''+str(ppn[i])+'''
#PBS -l mem=2gb 
#PBS -l walltime=30:00
#PBS -q coc-ice-multi

cd /nv/coc-ice/sliao7/pa3-template
module load gcc mvapich2/2.2
mpirun -np '''+ str(p)+''' ./jacobi -n ''' + str(n)+ ''' -d '''+str(d)+''' >> '''+instance+'''.txt
mpirun -np '''+ str(p)+''' ./jacobi -n ''' + str(n)+ ''' -d '''+str(d)+''' >> '''+instance+'''.txt
mpirun -np '''+ str(p)+''' ./jacobi -n ''' + str(n)+ ''' -d '''+str(d)+''' >> '''+instance+'''.txt
mpirun -np '''+ str(p)+''' ./jacobi -n ''' + str(n)+ ''' -d '''+str(d)+''' >> '''+instance+'''.txt
mpirun -np '''+ str(p)+''' ./jacobi -n ''' + str(n)+ ''' -d '''+str(d)+''' >> '''+instance+'''.txt
mpirun -np '''+ str(p)+''' ./jacobi -n ''' + str(n)+ ''' -d '''+str(d)+''' >> '''+instance+'''.txt
mpirun -np '''+ str(p)+''' ./jacobi -n ''' + str(n)+ ''' -d '''+str(d)+''' >> '''+instance+'''.txt
mpirun -np '''+ str(p)+''' ./jacobi -n ''' + str(n)+ ''' -d '''+str(d)+''' >> '''+instance+'''.txt
mpirun -np '''+ str(p)+''' ./jacobi -n ''' + str(n)+ ''' -d '''+str(d)+''' >> '''+instance+'''.txt
mpirun -np '''+ str(p)+''' ./jacobi -n ''' + str(n)+ ''' -d '''+str(d)+''' >> '''+instance+'''.txt
'''
				f.write(s)
				f.close()

