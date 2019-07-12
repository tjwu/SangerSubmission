import sys, os, csv
import sqlite3 


F = csv.reader(open("CodexExportForSanger.csv", "r"), delimiter = ",")

n = 0 

for i in F:
    n += 1
    if n == 1:
        continue
    batch_site, site, externalID, category, gene, chr, pos, ref, var = i[0:9]
    batchID = batch_site[:5]

    DB="/hgsccl/databases/flowcellStatus"
    DB_done_list = []
    db = sqlite3.connect(DB)
    cursor = db.cursor()
    cursor.execute(''' select sampleFlowcellLaneBarcode from samples s join flowcells f on f.id = s.flowcellID where samplename=? and f.batchID=?''', (externalID, batchID))
    data = cursor.fetchall()
    try:
        LaneBarcode = '-'.join(str(data[0][0]).split('-')[-2:])
    except IndexError:
        print externalID, batchID
        continue
    First = batch_site + '-' + LaneBarcode
    TempList= [ First, chr, pos, ref , var, "Heterozygous", "yer", "no", "no" ]
    print '\t'.join(TempList)


    
print n 
