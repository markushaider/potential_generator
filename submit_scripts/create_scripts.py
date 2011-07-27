#!/usr/bin/python

clusterNumbers = ['600', '601', '602', '603', '604', '605', '607', '608', '609', '610']
resolutions = ['32','64','128','256']
#clusterNumbers = ['600']
#resolutions = ['32']

for res in resolutions:
	for cNumber in clusterNumbers:
		file = open('submit'+cNumber+'_'+res+'.sh',"w")
		file.write('#!/bin/bash\n')
		file.write('../Potcomo_'+res+ \
	' /RAID3/markus/clusterdata/potentials/snapshots/'+cNumber+ \
	'/ /RAID3/markus/clusterdata/potentials/'+res+'/'+cNumber+'new/\n')
		file.close()

for res in resolutions:
	for cNumber in clusterNumbers:
		file = open('start_'+cNumber+'_'+res+'.sh',"w")
		file.write('#!/bin/bash\n')
		file.write('\n')
		file.write('#$ -q intel.q\n')
		file.write('#$ -o output_'+cNumber+'_'+res+'.out\n')
		file.write('#$ -cwd\n')
		file.write('#$ -j yes\n')
		file.write('./submit'+cNumber+'_'+res+'.sh\n')
		file.close()


#submit script
file = open('all_submit.sh','w')
file.write('#!/bin/bash\n')
for res in resolutions:
	for cNumber in clusterNumbers:
		file.write('qsub start_'+cNumber+'_'+res+'.sh\n')
		file.write('sleep 0.5\n')
file.close()
