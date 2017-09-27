# By Tsung-Jung Wu
# python Sanger_temp.py 170105_D00341_0512_AH7LYGBCXY

import sys
import os
import xml.etree.ElementTree as ET
import csv
import datetime
import logging
import glob
import re
import sqlite3
import shutil


# logging
logging.basicConfig(filename='SangerPush_Run.log',format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s', datefmt='%Y/%m/%d %H:%M:%S',level=logging.DEBUG)

dateString = datetime.date.today()
dateTime = datetime.datetime.now().strftime("%Y-%b-%d %H:%M:%S")
PreviousdateString = datetime.date.today()- datetime.timedelta(days=1)

#FC = sys.argv[1]
# PATHs
INSTRUMENT_DIR="/hgsccl/next-gen/Illumina/Instruments/"
FCD_path="../../Data/Intensities/BaseCalls/FCDefinition.xml"
python="/hgsccl/codified/bin/python2.7" 
#Codex="/hgsccl/codified/bin/exportCodex"
Codex="/hgsccl/codified/bin/exportCodex.manual"
FC_Need_Logging=[]


def Read_FC_database():
    # Get FCs which have done Mercury analysis
    DB="/hgsccl/databases/flowcellStatus"
    DB_done_list = []
    db = sqlite3.connect(DB)
    cursor = db.cursor()
    # FC has finished time should mean it has done Mercury analysis. 
    cursor.execute('''SELECT flowcellName FROM flowcells WHERE whenFinished is not null AND id >=22''')
    for row in cursor:
            DB_done_list.append(row[0])
    return DB_done_list

# After get FCs which have done Mercury, get FCs which have done Sanger.
# Diff Mercury finished FC and Sanger finished FC. The rest is the 
def Get_New_FC():
    #SangerDB="/hgsccl_software/prod/SangerSubmission_allfiles/SangerRecord_db"
    SangerDB="/hgsccl/databases/SangerRecord_db"
    Sanger_FC_Lab_list = []
    Sanger_FC_done_list = []
    sangerdb = sqlite3.connect(SangerDB)
    sangerfccursor = sangerdb.cursor()
    sangerfccursor.execute('''SELECT flowcellName, SamplingLab FROM Sanger_FC''')

    for sangerfcrow in sangerfccursor:
        FC,Lab = sangerfcrow
        FC_Lab='_'.join(sangerfcrow)
        Sanger_FC_done_list.append(FC)

    #NewFC = list(set(DB_done_list) - set(Sanger_FC_done_list))
    NewFC=DB_done_list
    #NewFC.remove("160328_D00341_0445_AHHNT5BCXX")
    #NewFC.remove("160330_D00222_0435_BHHGWFBCXX")
    #NewFC.remove("160403_D00222_0436_AHMF2CBCXX")
    return NewFC

# Run throug new FCs and get FCD info. This step gets BatchID and Lablist.
def Read_FCDef_Get_Lab(FC):
    # Get FC dir. 
    FC_DIR= glob.glob(INSTRUMENT_DIR + "*/*%s/Results/Project_*%s/" % (FC, FC))[0]
    os.chdir(FC_DIR)
    
    # Get FCD path ready to read it. 
    Read_XML = glob.glob('%s/*/*%s/Data/Intensities/BaseCalls/FCDefinition.xml' % (INSTRUMENT_DIR, FC))
    #print Read_XML, "<-reading xml"
    document=ET.parse(Read_XML[0])
    root = document.getroot()

    # Start reading it. 
    ele = root.find('LaneBarcodeInfo')
    BatchID=root.attrib['BatchId']
    Lablist = []
    for i in ele.findall('LaneBarcode'):
        SampleID = i.attrib['ID']
        Lab = i.attrib['sampling_lab']
        ##############################################################
        ## Get information for each sample # check if sample is IDMB?
        ##############################################################
        ##if SampleID.find("IDMB") >=0:
        ##    print "This is a target pannel sample"

        # Each lab should have one file. So, get all labs for later file reading purpose.
        #Lab = re.sub("\s+", '_', Lab)
        #Lab.replace("'","\\'").replace(" ", "\ ")
        if Lab.find(" ") >=0:
            continue
        else:
            Lab = re.sub('[^A-Za-z0-9]+', '-', Lab)

        if Lab in Lablist:
            continue
        else:
            Lablist.append(Lab)
    return [BatchID,Lablist]

# Each lab should generate a tsv file in codexExport directory. 
# Add them all into a list. Later, read through it and get samples which needs sanger sequencing. 
def Get_new_codex_files(FC, BatchID, Lablist):
    NewFilelist=[]
    #print FC, BatchID, Lablist
    for site in Lablist:
        codexFileName = "/hgsccl/codified/codexExport/%s-%s.%s.tsv" % (BatchID, site, dateString)
        #codexFileName = "/hgsccl/codified/webreports/%s-%s.%s.tsv" % (BatchID, site, dateString)
        In = BatchID + "-" + site
        CMD = python + " " + Codex + " " + In
        print CMD
        os.system(CMD)
        NewFilelist.append(codexFileName)

        Current_FC_Lab = FC+"_"+site
        #print Current_FC_Lab
        SangerFClist=[ FC, BatchID, dateString, site ]
        # Add FC and sampling lab to database. 
        # CREATE TABLE Sanger_FC(id INTEGER PRIMARY KEY, flowcellName TEXT, BatchID TEXT, Finish_Time BLOB, SamplingLab TEXT);

        #SangerDB="/hgsccl_software/prod/SangerSubmission_allfiles/SangerRecord_db"
        SangerDB="/hgsccl/databases/SangerRecord_db"
        Sanger_FC_Lab_list = []
        sangerdb = sqlite3.connect(SangerDB)
        sangerfccursor = sangerdb.cursor()
        sangerfccursor.execute('''SELECT flowcellName, SamplingLab FROM Sanger_FC WHERE FINISH_Time is not null''')

        for sangerfcrow in sangerfccursor:
            FC_Lab='_'.join(sangerfcrow)
            Sanger_FC_Lab_list.append(FC_Lab)

        if Current_FC_Lab not in Sanger_FC_Lab_list:
            sangerfccursor.execute('''INSERT INTO Sanger_FC( flowcellName, BatchID, Finish_Time, SamplingLab ) VALUES(?,?,?,?) ''', (SangerFClist))
            sangerdb.commit()
    return NewFilelist

# Reading through each codexgenerated file. 
# Get those samples and positions which need sanger sequencing. Put them all into one csv file. Mail this file to receivers.  
def Get_sample_need_sanger(FC, BatchID, NewFilelist):
    BAR = FC[20:]
    n=0
    Sanger_File_List=[]
    SampleVariantlist=[]
    for j in NewFilelist:
        OF= csv.reader(open(j,"r"), delimiter="\t")
        for k in OF:
            # Skip header line in case there are multiple FCs in one day. 
            if k[0] == "sample":
                Title=k
                continue
            else:
                if k[6] == "yes":
                    n+=1
                    Sanger_File_Name = "/hgsccl/SangerSubmission_allfiles/Sanger_Sequencing_%s_%s_%s.csv" % (FC, BatchID, dateString)
                    Wsanger=csv.writer(open(Sanger_File_Name,"a"), delimiter= ",", lineterminator="\n")
                    if n == 1:
                        Wsanger.writerow(Title)
                    Wsanger.writerow(k)
                    # Add samples which need sanger into a list for later use. Commnad it out since Diviay wants to get new ones only. TJ Sep 21 2017.
                    #SampleVariantlist.append(k[0:3]+[k[5]])
                    #CREATE TABLE Sanger_Samples(id INTEGER PRIMARY KEY, FCLBC TEXT, Codex_SampleID TEXT,Chr INTEGER, Pos INTEGER, Ref TEXT, Var TEXT, Zygosity TEXT, Report TEXT, Comment TEXT)
                    #Sample,chr pos ref var zygosity    report  confirmedWithSanger confirmedWithoutSanger  comment
                    #IR038-Columbia-1-IDMB52 13  48947542    A   T   Heterozygous    no  no  no
                    #Sanger_File_List.append(Sanger_File_Name)
                    Lane, IDMB = k[0].split('-')[-2:]
                    FCLBC = '-'.join([BAR, Lane, IDMB])
                    insert_Sample_list = [FCLBC]+k[0:7]+[dateTime]    
                    DB_key='_'.join([FCLBC]+k[1:3])

                    #SangerDB="/hgsccl_software/prod/SangerSubmission_allfiles/SangerRecord_db"
                    SangerDB="/hgsccl/databases/SangerRecord_db"
                    Sanger_Samples_list = []
                    sangerdb = sqlite3.connect(SangerDB)
                    sangercursor = sangerdb.cursor()
                    sangercursor.execute('''SELECT FCLBC, Chr, Pos FROM Sanger_Samples''')
                    for sangerrow in sangercursor:
                        DB_FCLBC, DB_Chr, DB_Pos = sangerrow
                        DB_value_key='_'.join([DB_FCLBC, str(DB_Chr), str(DB_Pos)])
                        Sanger_Samples_list.append(DB_value_key)
                    # Whatever in database, it should have done sanger sequencing.
                    # Whatever not in database, treat it as new and put into list for later use. 
                    if DB_key not in Sanger_Samples_list:
                        sangercursor.execute('''INSERT INTO Sanger_Samples(FCLBC, Codex_SampleID, Chr, Pos, Ref, Var, Zygosity, Report, Date_Time ) VALUES(?,?,?,?,?,?,?,?,?) ''', (insert_Sample_list))
                        sangerdb.commit()
                        # Including Zygosity here since later output needs this value
                        SampleVariantlist.append(k[0:3]+[k[5]])
                        if FC in FC_Need_Logging:
                            continue
                        else:
                            FC_Need_Logging.append(FC)
    if n >=1:
        mailString='mailx -a %s -s mailing_file_%s tsungjuw@bcm.edu  <<-EOF' % (Sanger_File_Name, BatchID)
        #os.system(mailString)
    return SampleVariantlist

# Check if all positions have vcf. If so, get the vcf. If not, run neptune
# Will add to here. 
# Coming soon when I have time.

def Get_Neptune_vcf(FC, SampleVariantlist):
    BAR = FC[20:]
    Count=0
    ref_count=0
    #Vcf_Missing_Count=0
    #Vcf_Missing_FileName="/hgsccl/SangerSubmission_allfiles/VCF_missing_Variant_list_%s.csv" % dateString
    #Vcf_Missing_File=csv.writer(open((Vcf_Missing_FileName), "a"), delimiter=',', lineterminator='\n')
    for l in SampleVariantlist:
        Count+=1
        if len(l[0].split("-")) > 4:
            One, Two, Three, Four, Five = l[0].split("-")
            BATCH = One
            Lab = Two+"-"+Three
            Lane = Four
            IDMB = Five
        else:
            BATCH, Lab, Lane, IDMB = l[0].split("-")

            
        WriteVCF=csv.writer(open(("/hgsccl/SangerSubmission_allfiles/Sanger_Sequencing_%s_%s_%s.vcf" % (FC,BATCH,dateString)), "a"), delimiter="\t")

        chr, pos, zygosity = l[1:]
        try:
            SampleBWA = glob.glob(INSTRUMENT_DIR + "*/*%s/Results/Project_*%s/Sample_*-%s-%s/BWA*" % (FC,FC,Lane,IDMB))[0]
        except:
            print "Assume sample is still running early analysis steps"
            continue
        GetID_cmd= "grep INDIVIDUAL_ID %s | cut -d= -f2" % SampleBWA
        externalID=os.popen(GetID_cmd).read().strip()
        LaneBarcode=Lane+"-"+IDMB
        try:
            SampleNeptuneVCF = max(glob.iglob(INSTRUMENT_DIR + "*/*%s/Results/Project_*%s/Sample_*-%s-%s/neptune/*vcf" % (FC,FC,Lane,IDMB)),key=os.path.getctime)
        except ValueError:
            print "Assume sample is still gneratingvcf"
            continue
        ReadVCF = csv.reader(open(SampleNeptuneVCF,"r"),delimiter="\t")
        for m in ReadVCF:
            if Count == 1:
                if m[0].find("#") >=0:
                    WriteVCF.writerow(m)
                else:
                    if ref_count == 0:
                        ref_count+=1
                        if m[0]==chr and m[1] ==pos:
                            newINFO = m[7]+";externalID(INDIVIDUAL_ID)=%s;LaneBarcode=%s;Batch-Lab=%s;Zygosity=%s" % (externalID,LaneBarcode,BATCH+"-"+Lab, zygosity)
                            m[7]= newINFO
                            WriteVCF.writerow(["#Reference=/stornext/snfs0/hgsc-refs/Illumina/bwa_references/h/hg19/original/hg19.fa"])
                            WriteVCF.writerow(m)
                        else:
                            WriteVCF.writerow(["#Reference=/stornext/snfs0/hgsc-refs/Illumina/bwa_references/h/hg19/original/hg19.fa"])
                    else:
                        if m[0]==chr and m[1] ==pos:
                            newINFO = m[7]+";externalID(INDIVIDUAL_ID)=%s;LaneBarcode=%s;Batch-Lab=%s;Zygosity=%s" % (externalID,LaneBarcode,BATCH+"-"+Lab, zygosity)
                            #newINFO = m[7]+";externalID(INDIVIDUAL_ID)=%s;LaneBarcode=%s" % (externalID,LaneBarcode)
                            m[7]= newINFO
                            WriteVCF.writerow(m)
            else:
                if m[0]==chr and m[1]==pos:
                    newINFO = m[7]+";externalID(INDIVIDUAL_ID)=%s;LaneBarcode=%s;Batch-Lab=%s;Zygosity=%s" % (externalID,LaneBarcode,BATCH+"-"+Lab, zygosity)
                    #newINFO = m[7]+";externalID(INDIVIDUAL_ID)=%s;LaneBarcode=%s" % (externalID,LaneBarcode)
                    m[7]= newINFO
                    WriteVCF.writerow(m)
    # If
    if ref_count ==0:
        ref_count+=1
        try:
            WriteVCF.writerow(["#Reference=/stornext/snfs0/hgsc-refs/Illumina/bwa_references/h/hg19/original/hg19.fa"])
        except UnboundLocalError:
            pass
    #if Vcf_Missing_Count > 0:
    #    mailString='mailx -a %s -s mailing_file_%s tsungjuw@bcm.edu  <<-EOF' % (Sanger_File_Name, BatchID)
        #os.system(mailString)

IGNFC=open("/hgsccl/SangerSubmission_allfiles/IgnoreFC.txt","r")
Exceptionlist=[]
for i in IGNFC.readlines():
    i = i.strip()
    Exceptionlist.append(i)

# Moveing old files to old files' directory to avoid confusion.

#Sanger_prior_dir = "/hgsccl/SangerSubmission_allfiles/Prior_%s_files" % dateString
#if not os.path.isdir(Sanger_prior_dir):
#    os.mkdir(Sanger_prior_dir)
#Sanger_prior_files = glob.glob("/hgsccl/SangerSubmission_allfiles/Sanger_Sequencing*.*")
#
#for file in Sanger_prior_files:
#    dest_dir = Sanger_prior_dir
#    shutil.move(file, dest_dir )
#
#Codex_prior_dir = "/hgsccl/codified/codexExport/Prior_%s_files" %dateString
#if not os.path.isdir(Codex_prior_dir):
#    os.mkdir(Codex_prior_dir)
#Codex_prior_files = glob.glob("/hgsccl/codified/codexExport/*.tsv")
#
#for file in Codex_prior_files:
#    dest_dir = Codex_prior_dir
#    shutil.move(file, dest_dir )
#

def Archive_prior_file(Target_dir, Target_format):
    Prior_dir = Target_dir + "Prior_%s_files" % dateString
    if not os.path.isdir(Prior_dir):
        os.mkdir(Prior_dir)
    Prior_files = glob.glob(Target_dir+Target_format)

    for file in Prior_files:
        dest_dir = Prior_dir
        filename=file.split("/")[-1]
        # dest_dir+"/"+filename
        if os.path.isfile(dest_dir+"/"+filename):
            os.rename(dest_dir+"/"+filename, dest_dir+"/old_"+filename)
            shutil.move(file, dest_dir)
        else:
            shutil.move(file, dest_dir)
        
Archive_prior_file("/hgsccl/codified/codexExport/","*.tsv")
Archive_prior_file("/hgsccl/SangerSubmission_allfiles/","Sanger_Sequencing*.*")



DB_done_list=Read_FC_database()
NewFC_list=Get_New_FC()
NewFC_list.sort()
#NewFC_list=["161110_D00341_0503_AH3N3TBCXY"]
#NewFC_list=["170127_D00341_0514_AHFG37BCXY"]
#NewFC_list=["170301_D00222_0498_BHG72GBCXY"]
for FC in NewFC_list:
    #logging.info(FC)
    if FC in Exceptionlist:
        continue
    BatchID, Lablist = Read_FCDef_Get_Lab(FC)
    NewFilelist = Get_new_codex_files(FC, BatchID, Lablist)
    SampleVariantlist = Get_sample_need_sanger(FC, BatchID, NewFilelist)
    Get_Neptune_vcf(FC,SampleVariantlist)

if len(FC_Need_Logging) == 0:
    logging.info("")
else:
    For_logging=','.join(FC_Need_Logging)
    logging.info(For_logging)


#Copy all files from SangerSubmission_allfiles to SangerSubmission directory Diviya wants directory looks like this. 
def Move_VCF():
    Allfiles=glob.glob("/hgsccl/SangerSubmission_allfiles/*")
    Target_Dir="/hgsccl/SangerSubmission/"
    for i in Allfiles:
        if i.endswith("vcf"):
            shutil.copy2(i, Target_Dir)
Move_VCF()


#SampleVariantlist=Get_sample_need_sanger("170330_D00222_0504_AHG3W3BCXY","IR038",["/hgsccl/codified/codexExport/IR038-Columbia.2017-06-22.tsv","/hgsccl/codified/codexExport/IR038-Vanderbilt.2017-06-22.tsv","/hgsccl/codified/codexExport/IR038-CHOP.2017-06-22.tsv"])
#Get_Neptune_vcf("170330_D00222_0504_AHG3W3BCXY",SampleVariantlist)
