#!/usr/bin/python

import os

clusterNumbers = ['600', '601', '602', '603', '604', '605', '607', '608', '609', '610']
resolutions = ['32','64','128','256']
#clusterNumbers = ['600', '601']
#resolutions = ['32','64']


#for res in resolutions:
#	path='/RAID3/markus/clusterdata/potentials/'+res+'/'
#	os.mkdir( path);
#	print 'created '+path+'\n'

for res in resolutions:
	for cNumber in clusterNumbers:
		path='/RAID3/markus/clusterdata/potentials/'+res+'/'+cNumber+'new/'
		os.mkdir( path);
		print 'created '+path+'\n'

