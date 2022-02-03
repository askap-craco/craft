#!/usr/bin/env python

from pyclustering.cluster.dbscan import dbscan
import numpy as np
import sys

# Index values for various candidate fields
sn  = 0     # S/N
t   = 1     # sampno
s   = 2     # secs from file start
w   = 3     # boxcar width
idt = 4     # number of samples difference between top and bottom of frequency range
dm  = 5     # DM
bno = 6     # Beam number
mjd = 7     # Modified Julian Date

def _main():
	if len(sys.argv) < 4:
		print("USAGE: ./run_DBSCAN.py <.cand file> <eps> <nmin>")
	else:
		fname = sys.argv[1]
		cands = open_file(fname)


		abbr_cands = np.array([ (cand[t], cand[w], cand[dm]) for cand in cands ])

		eps = float(sys.argv[2])
		nmin = int(sys.argv[3])

		DBSCAN = dbscan(abbr_cands, eps, nmin, True)
		DBSCAN.process()
		cluster_idxs = DBSCAN.get_clusters()

		clusters = [ cands[idxs] for idxs in cluster_idxs ]

		#get the brightest cand in each cluster to represent the cluster
		cluster_reps = [ cluster[np.argmax(cluster[:,sn])] for cluster in clusters ]

		#filter the reps by DM and width
		cluster_reps = [ rep for rep in cluster_reps if rep[dm] >= 100 and rep[w] <= 10 ]

		write_cands(fname, cluster_reps)

def open_file(fname, sort_idx=sn, ret_nparr=True):
	cands = []
	with open(fname, 'r') as f:
		for i, line in enumerate(f):
			if line[0] != '#' and len(line) > 5:
				# In case the file has more columns than we need, trim the extras off the end
				new_cand = list(map(float, line.split()[0:mjd+1]))
				cands.append(new_cand)
				
				# Sometimes there's no mjd field by default, so we need to add it
				while len(new_cand) <= mjd:
					new_cand.append(0.0)

	cands.sort(key=lambda x: x[sort_idx])

	if ret_nparr:
		return np.array(cands)
	else:
		return cands

def write_cands(fname, cands, suffix='.dbs'):
	header = 'S/N, sampno, secs from file start, boxcar, idt, dm, beamno, mjd'
	intf = '%d'
	floatf = '%0.2f'
	formats = (floatf, intf, floatf, intf, intf, floatf, intf, '%0.15f')
	#npcands = np.array(cands)
	np.savetxt(fname+suffix, cands, fmt=formats, header=header)

if __name__ == '__main__':
	_main()
