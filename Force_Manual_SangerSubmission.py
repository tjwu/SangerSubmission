# Christie's request or from DM directly.
# This script runs on manually generated list which looks like codified list. 
# Then, read through it and got into each sample's directory to get neptune vcf.
# Finally, output each batch's vcf which includes some more info to /hgsccl/SangerSubmission.
# Divya's script is able to read this directory and get vcf for later sanger checking.
# python SangerSubmission.py
# By Tsung-Jung Wu

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
RefLine="#Reference=/stornext/snfs0/hgsc-refs/Illumina/bwa_references/h/hg19/original/hg19.fa"
All_Files = "/hgsccl/SangerSubmission_allfiles/"

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
    NewFC = DB_done_list
    return NewFC

# Run throug new FCs and get FCD info. This step gets BatchID and Lablist.
def Read_FCDef_Get_Lab(FC):
    # Get FC dir. # Add exception here to handle no results directory 
    try:
        FC_DIR= glob.glob(INSTRUMENT_DIR + "*/*%s/Results/Project_*%s/" % (FC, FC))[0]
    except IndexError:
        return False
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

# Codex command modify to output all samples need to report at once. Modify function to handle this change. # Apr 5 2018
def RunCodexCMD():
    codexFileName =  "%s/Codex_Report_Out_%s.tsv" % (All_Files, dateString)
    CMD = python + " " + Codex + " > " + codexFileName
    print CMD                                                                                 
    #os.system(CMD)
    codexFileName = "/hgsccl_software/devel/TJ/SangerSubmission/FakeCodexFile.tsv" 
    return codexFileName

# Add them all into a list. Later, read through it and get samples which needs sanger sequencing. 
def Get_new_codex_files(FC, BatchID, Lablist):
    NewFilelist=[]

    #print FC, BatchID, Lablist
    for site in Lablist:
        #codexFileName = "/hgsccl/codified/codexExport/%s-%s.%s.tsv" % (BatchID, site, dateString)
        #In = BatchID + "-" + site
        #CMD = python + " " + Codex + " " + In
        #print CMD
        #os.system(CMD)
        #NewFilelist.append(codexFileName)

        Current_FC_Lab = FC+"_"+site
        SangerFClist=[ FC, BatchID, dateString, site ]

        # Add FC and sampling lab to database. 
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
            #sangerdb.commit()
    return NewFilelist

def Get_new_vip_files(BatchID, Lablist, BAR):
    NewVIPlist = []
    FC_Dir = glob.glob("%s/*/*%s/Results/Project_*%s/" %  (INSTRUMENT_DIR, BAR, BAR))[0]
    try:
        FC_Dir_Vip = max(glob.glob("%s/*/*%s/Results/Project_*%s/%s*.vip.txt" %  (INSTRUMENT_DIR, BAR, BAR, BatchID)), key=os.path.getctime)
        NewVIPlist.append(FC_Dir_Vip)
    except ValueError:
        resultsDir = "/hgsccl/SangerSubmission_allfiles/"
        for site in Lablist:
            vipFileName = "%s/%s.%s.vip.txt" % (resultsDir, BatchID, dateString)
            CMD = "/hgsccl/codified/bin/vips \"%s-%s\" > %s" % (BatchID, site, vipFileName)
            #os.system(CMD)
            NewVIPlist.append(vipFileName)
            shutil.copy2(vipFileName, FC_Dir)
    return NewVIPlist

def Get_VIP():
    NewVIPlist = []
    resultsDir = "/hgsccl/SangerSubmission_allfiles/"
    vipFileName = "%s/All-%s-vip.txt" % (resultsDir, dateString)
    CMD = "/hgsccl/codified/bin/vips > %s" % vipFileName
    #os.system(CMD)
    NewVIPlist.append(vipFileName) 
    return NewVIPlist

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
            elif len(k) < 9: # Add here to avoid unexpected line show up.
                continue
            elif k[0].find(BatchID) < 0:
                continue
            else:
                if k[6] == "yes" and k[8] == "no":
                    n+=1
                    Sanger_File_Name = "/hgsccl/SangerSubmission_allfiles/Sanger_Sequencing_%s_%s_%s.csv" % (FC, BatchID, dateString)
                    Wsanger=csv.writer(open(Sanger_File_Name,"a"), delimiter= ",", lineterminator="\n")
                    if n == 1:
                        Wsanger.writerow(Title)
                    Wsanger.writerow(k)
                    # Add samples which need sanger into a list for later use. Commnad it out since Diviya wants to get new ones only. TJ Sep 21 2017.
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
                        #sangerdb.commit()
                        # Including Zygosity here since later output needs this value
                        SampleVariantlist.append(k[0:6])
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

def Check_missing(InRow):
    MissingValue = 0
    MissingList = []
    Terms = ["Gene=", "Transcript_Used=", "Location=", "Nucleotide_corrected=", "AminoAcid="]
    for t in Terms:
        if InRow[7].find(t)<0:
            InRow[7] = InRow[7] + ( ";%s" % t)
            MissingList.append(t)
    return InRow        
            #InRow[7] = InRow[7]+";%s=N/A" % t
    #if len(MissingList) > 0:
    #    return MissingList
    #else:
    #    return "Pass"

def FindAndWrite(m, externalID, LaneBarcode, BATCH, Lab, zygosity, WriteVCF):
    #if m[7].find("Location")<0:
    #    m[7] = m[7].replace("Nucleotide_corrected","Location=Exon N/A;Nucleotide_corrected")
    #elif m[7].find("Transcript_Used")<0:
    #    m[7] = m[7].replace(";Location=",";Transcript_Used=N/A;Location=")
    #MSV = Check_missing(m)
    #if MSV != "Pass":
    #    print LaneBarcode, " missing values", MSV
    m = Check_missing(m)
    VR, RR, DP = m[9].split(":")[1:4]
    VRDP_ratio = round(float(VR)/float(DP), 2)
    if VRDP_ratio < 0.3: # When allele frequency is below 0.3, put in "Possible Mosaic" for the zygosity. Eric and Yunyun's email.
        zygosity = "Possible Mosaic"
    newINFO = m[7]+";externalID(INDIVIDUAL_ID)=%s;LaneBarcode=%s;Batch-Lab=%s;Zygosity=%s" % (externalID,LaneBarcode,BATCH+"-"+Lab, zygosity)
    m[7]= newINFO
    WriteVCF.writerow(m)

def RescueAndWrite(row, externalID, LaneBarcode, BATCH, Lab, zygosity, ReadDepth, SampleName, WriteVCF):
    row = row[:-2] + ReadDepth
    #if row[7].find("Location")<0:
    #    row[7] = row[7].replace("Nucleotide_corrected","Location=Exon N/A;Nucleotide_corrected")
    #elif row[7].find("Transcript_Used")<0:
    #    row[7] = row[7].replace(";Location=",";Transcript_Used=N/A;Location=")
    #MSV = Check_missing(row)
    #if MSV != "Pass" :
    #    print SampleName, MSV
    row = Check_missing(row)
    GOF, GQ = row[9].split(":")[2:4]
    VRDP_ratio = round(float(GOF)/float(GQ), 2)
    if VRDP_ratio < 0.3:
        zygosity = "Possible Mosaic"
    newINFO = row[7]+";sampleID=%s;externalID(INDIVIDUAL_ID)=%s;LaneBarcode=%s;Batch-Lab=%s;Zygosity=%s" % (SampleName, externalID, LaneBarcode, BATCH+"-"+Lab, zygosity)
    row[7]= newINFO
    WriteVCF.writerow(row)

def ReadSampleVCF(BAR, Lane, IDMB, chr, pos, INSTRUMENT_DIR):
    SampleVCFList = glob.glob("%s/*/*%s/Results/Project_*%s/Sample_%s-%s-%s/*/*s.vcf" %  (INSTRUMENT_DIR, BAR, BAR, BAR, Lane, IDMB))
    for vcf in SampleVCFList:
        Rvcf = csv.reader(open(vcf, "r"), delimiter= "\t")
        for RV in Rvcf:
            if RV[0] == chr and RV[1] == pos:
                Def, Val = RV[8:]
                ReadDepth = [Def, Val]
    if not ReadDepth:
        print "Create empty list for later FCs can run"
        ReadDepth = []
    return ReadDepth

def GetSampleName (BAR, Lane, IDMB, INSTRUMENT_DIR):
    BWAConfig = glob.glob("%s/*/*%s/Results/Project_*%s/Sample_%s-%s-%s/BWAConfigParams.txt" %  (INSTRUMENT_DIR, BAR, BAR, BAR, Lane, IDMB))[0]
    ReadBWA = csv.reader(open(BWAConfig, "r"), delimiter="\t")
    #ReadBWA = open(BWAConfig, "r")
    for i in ReadBWA:
        if i[0].find("SAMPLE_NAME") >=0:
            SampleName = i[0][12:].strip()
    #print SampleName
    return SampleName

def InsertRefLine(OutputVCF, HeaderList):
    f = open(OutputVCF, "r")
    contents = f.readlines()
    f.close()

    newfile = HeaderList + contents

    f = open(OutputVCF, "w")
    newfile = "".join(newfile)
    f.write(newfile)
    f.close()

def Get_Neptune_vcf(FC, SampleVariantlist, NewVIPlist):
    print FC
    BAR = FC[20:]
    Count=0
    HeaderList=[]
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

        VCF_for_Write = ("/hgsccl/SangerSubmission_allfiles/Sanger_Sequencing_%s_%s_%s.vcf" % (FC,BATCH,dateString))

        WriteVCF = csv.writer(open(VCF_for_Write , "a"), delimiter="\t")

        chr, pos, ref, var, zygosity = l[1:]
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
        VCFline = 0
        VCFcontent = 0
        DoneFinding = 0
        for m in ReadVCF:
            VCFcontent += 1
            VCFline += 1
            if m[0].find("#") >= 0:
                if len(HeaderList) < 3 and VCFline < 3:
                    HeaderList.append(m[0] + "\r\n")
                elif VCFline == 3 and len(HeaderList) < 3:
                    m[9] = ""
                    newm = '\t'.join(m)
                    HeaderList.append(newm +"\r\n")
                    HeaderList.append(RefLine + "\r\n")
            else:
                #if m[0] == chr and m[1] == pos and m[3] == ref and m[4] == var:
                if m[0] == chr and m[1] == pos and m[3] ==ref and m[4] ==var:
                    FindAndWrite(m, externalID, LaneBarcode, BATCH, Lab, zygosity, WriteVCF)
                    DoneFinding += 1
                #else:
                    #missingCount += 1
                    #FCLBC = "%s-%s-%s" % (BAR, Lane, IDMB)
                    #ReadDepth = ReadSampleVCF(BAR, Lane, IDMB, chr, pos, INSTRUMENT_DIR)
                    #for vl in NewVIPlist:
                    #    ReadVL = csv.reader(open(vl, "r"), delimiter = "\t")
                    #    for v in ReadVL:
                    #        if v[0] == chr and v[1] == pos and v[3] == ref and v[4] == var:
                    #            RescueAndWrite(v, externalID, LaneBarcode, BATCH, Lab, zygosity, ReadDepth, WriteVCF)
        #if VCFcontent <= 3 or missingCount >0:
        if  DoneFinding == 0:
            FCLBC = "%s-%s-%s" % (BAR, Lane, IDMB)
            #print "%s neptune VCF in %s %s use VIP DB annotation" % (FCLBC, BATCH, Lab)
            
            ReadDepth = ReadSampleVCF(BAR, Lane, IDMB, chr, pos, INSTRUMENT_DIR)
            SampleName = GetSampleName (BAR, Lane, IDMB, INSTRUMENT_DIR)

            for vl in NewVIPlist:
                ReadVL = csv.reader(open(vl, "r"), delimiter = "\t")
                for v in ReadVL:
                    #if v[0] == chr and v[1] == pos and m[3] == ref and m[4] == var:
                    if v[0] == chr and v[1] == pos and v[3] == ref and v[4] == var:
                        #print v[0:5], "second"
                        RescueAndWrite(v, externalID, LaneBarcode, BATCH, Lab, zygosity, ReadDepth, SampleName, WriteVCF)
    try:
        InsertRefLine(VCF_for_Write, HeaderList)
    except UnboundLocalError:
        pass

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
#Archive_prior_file("/hgsccl/SangerSubmission_allfiles/","Codex_*.*")
#Archive_prior_file("/hgsccl/SangerSubmission_allfiles/","All-*.*")


DB_done_list=Read_FC_database()
NewFC_list=Get_New_FC()
NewFC_list.sort()
AllFilelist = [RunCodexCMD()]
NewVIPlist = Get_VIP()

#NewFC_list=["170615_D00222_0510_BHLHLVBCXY","170818_D00222_0516_AHMYKVBCXY"]
#NewFC_list = ['171228_D00341_0578_BH5GM2BCX2'] 
#NewFC_list = ['171020_D00222_0528_BHYFC2BCXY'] 
#NewFC_list = ['170803_D00341_0538_AHL3L3BCXY'] 
#NewFC_list = ['171011_D00341_0549_BHYGMMBCXY'] 
#NewFC_list = ['170809_D00341_0539_AHMWYNBCXY'] 
   
for FC in NewFC_list:
    #logging.info(FC)
    BAR = FC[20:]
    if FC in Exceptionlist:
        continue
    # Handel error here if no results directory
    try:
        BatchID, Lablist = Read_FCDef_Get_Lab(FC)
    except TypeError:
        continue
    #NewFilelist = Get_new_codex_files(FC, BatchID, Lablist) 
    Get_new_codex_files(FC, BatchID, Lablist)

    #NewVIPlist = Get_new_vip_files(BatchID, Lablist, BAR)
    
    SampleVariantlist = Get_sample_need_sanger(FC, BatchID, AllFilelist)
    Get_Neptune_vcf(FC,SampleVariantlist, NewVIPlist)

if len(FC_Need_Logging) == 0:
    logging.info("")
else:
    For_logging=','.join(FC_Need_Logging)
    logging.info(For_logging)

# VCF might be empty due to empty neptune vcf
def Remove_Empty_VCF():
    Allfiles = glob.glob("/hgsccl/SangerSubmission_allfiles/*vcf")
    print "Checking file content"
    for i in Allfiles:
        Final_Check=csv.reader(open(i, "r"),delimiter="\t")
        FileCheckRow=0
        for j in Final_Check:
            FileCheckRow+=1
        # 4 rows means this file only has title. Remove vcf but keep csv to avoid confusion
        if FileCheckRow == 4:
            os.remove(i)

#Copy all files from SangerSubmission_allfiles to SangerSubmission directory Diviya wants directory looks like this. 
def Move_VCF():
    Allfiles=glob.glob("/hgsccl/SangerSubmission_allfiles/*")
    Target_Dir="/hgsccl/SangerSubmission/"
    for i in Allfiles:
        if i.endswith("vcf"):
            shutil.copy2(i, Target_Dir)
    
Remove_Empty_VCF()
Move_VCF()

#SampleVariantlist=Get_sample_need_sanger("170330_D00222_0504_AHG3W3BCXY","IR038",["/hgsccl/codified/codexExport/IR038-Columbia.2017-06-22.tsv","/hgsccl/codified/codexExport/IR038-Vanderbilt.2017-06-22.tsv","/hgsccl/codified/codexExport/IR038-CHOP.2017-06-22.tsv"])
#Get_Neptune_vcf("170330_D00222_0504_AHG3W3BCXY",SampleVariantlist)
