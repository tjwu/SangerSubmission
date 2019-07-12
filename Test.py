import sys, os
import datetime

batchID = sys.argv[1]
site = sys.argv[2]
resultsDir = "/hgsccl_software/devel/TJ/SangerSubmission_allfiles" 

longDateString = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
vipFileName = "%s/%s.%s.vip.txt" % (resultsDir, batchID, longDateString)
command = "/hgsccl/codified/bin/vips \"%s-%s\" > %s" % (batchID, site, vipFileName)
print command
os.system(command)

