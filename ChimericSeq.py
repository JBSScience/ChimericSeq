# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 15:46:04 2015

@author: JBSUSER
"""

import subprocess,os,math,re,sys,zipfile,datetime,multiprocessing,webbrowser,threading,traceback,time,csv
import platform
import difflib
from tkinter import *
from tkinter import filedialog,ttk
import tkinter as tki
from threading import Thread,Lock
from Bio.Seq import Seq
from bisect import bisect
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

 
class Core:
    workingDirectory=os.getcwd().replace("\\","/")   #Base directory = D://ChimericSeq for example
    BUILDING=False
    doGTF=True
    addPathsBool=True
    outDir=os.getcwd().replace("\\","/")
    Bowtie2Folder=workingDirectory+'/bowtie2-2.2.9'
    ViralRefFolder=workingDirectory+'/Viral_Reference'
    HostRefFolder=workingDirectory+'/Host_Reference'
    ViralRefFa='None Selected'
    HostRefFa='None Selected'
    HostGTF='None Selected'
    readLocation=''
    readBaseName=''
    distanceThreshold=10000
    readSimilarity=95
    ntStretchCount=8
    ntUnidentifiedCount=10
    readDir=''
    useBasic=False
    GeneData=None
    FocusData=None
    ChromList=None
    folderMode=False
    Cores=[]
    coreOption=2
    logger=[]
    CSminLen=10
    STOP=False
    doReadClean=False 
    YESNO=None
    paired=False
    loadAlignments=False #change to false in final code
    overlapRequested=True
    Na=115
    HostPrefix='humanRef'
    viralPrefix='viralRef'
    GO=False
    autoMode=False #change to false in final code
    savedReads=[] # used by multiple files
    savedReadsHlength=[] # host length
    savedReadsMapData=[] # used by multiple files
    savedReadsFilename=""
    savedReadsFilenameFinal=""
    savedReadsLogFilename=""
    savedReadsDir=""
    multipleFileFinalReads=[]

    processLargeFile=False # batch mode
    processLargeFileDir="" # directory for files
    processingLargeFile = False
    processingLargeFileComplete = True
    processingSaveData=False
    processingLargeFileLock = threading.Lock()
    processingLargeFileCompleteLock = threading.Lock()
    processingSaveDataLock=threading.Lock()
    inProcesingLoop=False # processing one file
    inProcesingDirLoop=False # start to proces files
    processLargeFileIndex=0 # >=0 for index, -1 for all files processed
    processLargeFileCheckGTF=None
    processLargeFileInterrupt=None
    processLargeFileFirstRead=True
    splitInfoFilename="_split_info.dat"
    splitFileCount=0
    splitFile1=''
    splitFile2=''
    splitStartTime=0
    splitEndTime=0
    splitSelectedFiles=[] # specified in file
    splitSelectedFilesUI=[] # specified in UI
    # trim3
    Trim_max=70
    Trim_min=0
    Trim_avg=0
    Trim=Trim_avg
    # trim5
    Trim2_max = 70
    Trim2_min = 0
    Trim2_avg = 0
    Trim2=Trim2_avg
    
    otempFilter=70
    olenFilter=24
    htempFilter=25
    hlenFilter=10
    vtempFilter=25
    vlenFilter=10
    runmode=0 #0=integration mode, 1=translocation mode
    preprocessReadLength=10
    preprocessUnknownPercent=100
    microhomologystringency_selection=0
    microhomologystringency_min=80  # percentage
    microhomologystringency_avg=90
    microhomologystringency_max=100
    microhomologystringency=microhomologystringency_avg
    microhomologystringency1_min=15 # count
    microhomologystringency1_avg=20
    microhomologystringency1_max=35
    microhomologystringency1=microhomologystringency1_avg
    platformType = platform.system()

    def __init__(self):
        os.chdir(self.workingDirectory)  #change dir to working directory (ChimericSeq Root folder)
        if self.platformType == "Windows":                
            if self.addPathsBool:
                self.addPaths()
        for i in range (1,multiprocessing.cpu_count()+1,1):  #grab the number of CPUs and make list 2 through X
            self.Cores.append(i)
        self.lock=Lock()
        self.event=threading.Event()
        self.createInterface()
        
    def createInterface(self):
        self.topLevelTK=Tk()
        self.topLevelTK.withdraw()
        root=tki.Toplevel(self.topLevelTK)
        root.protocol('WM_DELETE_WINDOW',self.topLevelTK.destroy)
        self.interface=Interface(root,self)
        root.geometry("1200x555")
        self.app=self.interface
        self.app.printToLog('Logging Session')
        self.app.printToLog('Date Time: '+str(datetime.datetime.now()))
        self.app.printToLog('Number of CPU cores: '+str(multiprocessing.cpu_count()))
        self.printDefaults()
        root.mainloop()

    def resetLargeFileProcess(self):
        self.savedReads=[] # used by multiple files
        self.savedReadsHlength=[] # host length
        self.savedReadsMapData=[] # used by multiple files
        self.savedReadsFilename=""
        self.savedReadsFilenameFinal=""
        self.savedReadsLogFilename=""
        self.savedReadsDir=""
        self.multipleFileFinalReads=[]

        self.processLargeFile=False # batch mode
        self.processLargeFileDir="" # directory for files
        self.processingLargeFile = False
        self.processingLargeFileComplete = True
        self.processingSaveData=False
        self.processingLargeFileLock = threading.Lock()
        self.processingLargeFileCompleteLock = threading.Lock()
        self.processingSaveDataLock=threading.Lock()
        self.inProcesingLoop=False # processing one file
        self.inProcesingDirLoop=False # start to proces files
        self.processLargeFileIndex=0 # >=0 for index, -1 for all files processed
        self.processLargeFileCheckGTF=None
        self.processLargeFileInterrupt=None
        self.processLargeFileFirstRead=True
        self.splitInfoFilename="_split_info.dat"
        self.splitFileCount=0
        self.splitFile1=''
        self.splitFile2=''
        self.splitStartTime=0
        self.splitEndTime=0
        self.splitSelectedFiles=[]
        self.splitSelectedFilesUI=[]
        
    def hello(self):
        print('hello')
        self.printToLog('hola')
        #self.printToLog('hello')
    
    def loadGTF(self):
        if self.HostGTF!='None Selected':
            start=time.time()
            self.printToLog('Building Gene Database...')
            f=open(self.HostGTF,'r')
            focusData=[] #focusData
            chromList=[] #chromList
            genes=[]
            for line in f:
                test=line.split("\t")
                if test[0] not in chromList:
                    chromList.append(test[0])
                    focusData.append([])
                    genes.append([])
                    i=-1
                else:
                    i=chromList.index(test[0])
                try:
                    focusData[i].append(test)
                    if test[2]=='gene':
                        genes[i].append(test)
                except:
                    pass
            f.close()
            self.FocusData=focusData
            self.GeneData=genes
            self.ChromList=chromList
            stop=time.time()
            runtime=stop-start
            self.printToLog('Build complete at '+str(runtime)+ 'seconds')
            self.printToLog('======================================================================')
            self.GO=True
            self.event.set()
        else:
            if self.processLargeFile==False or  self.processLargeFileCheckGTF == None:
                self.printToLog("No selected Host GTF file to build on")
                self.printToLog('Gene info will not be available. Continue?')
                self.printToLog('======================================================================')
                self.app.yesbttn.configure(command=lambda: self.modGo(True,'set'))
                self.app.nobttn.configure(command=lambda: self.modGo(False,'set'))
            if self.processLargeFile==True:
                if self.processLargeFileCheckGTF == None:
                    self.processLargeFileCheckGTF = True
                    self.processLargeFileInterrupt = False
                else:
                    self.modGo(True,'set')
            #ask user yes or no buttons
    
    def modGo(self,arg,command):
        self.GO=arg
        if self.GO==False:
            self.printToLog('Run Cancelled')

            if self.processLargeFile==True:
                self.printToLog('Stop processing multiple files')
                self.changeProcessingLargeFileStatus(False)
                self.processLargeFileIndex = -1
                self.splitSelectedFiles=[]
                self.splitSelectedFilesUI=[]
                self.processLargeFileCheckGTF = None
                self.processLargeFileInterrupt = True
                self.changeProcessingLargeFileCompleteStatus(True)
            
        if command=='set':
            self.event.set()
            
    def installer(self,install,OS):
        os.chdir(self.workingDirectory)
        if 'blast' in install:
            os.chdir(self.workingDirectory)
            subprocess.call(['ncbi-blast-2.2.31+-win64.exe'])
        if 'bowtie' in install:
            z=zipfile.ZipFile('bowtie2-2.2.5-mingw-win64.zip','r')
            z.extractall(self.Bowtie2Folder)

    def printToLog(self,message):
        self.app.printToLog(message)
        
    def printDefaults(self):
        self.printToLog('Working Directory Set To:\n'+self.workingDirectory)
        self.printToLog('Bowtie2 Folder is set to:\n'+self.Bowtie2Folder)
        self.printToLog('Viral Reference File is set to:\n'+self.ViralRefFa)
        self.printToLog('Bowtie2 viral index prefix: '+self.viralPrefix)
        self.printToLog('Host Reference File is set to:\n'+self.HostRefFa)
        self.printToLog('Bowtie2 host index prefix: '+self.HostPrefix)
        self.printToLog('CPU option: '+str(self.coreOption)+' CPUs engaged')
        self.printToLog('Load Prexisting Alignments: '+str(self.loadAlignments))
        self.printToLog('Clipped Sequence Min Length: '+str(self.CSminLen))
        self.printToLog('Salt concentration: '+str(self.Na)+'mM')
        self.printToLog('Prompting: '+str(not self.autoMode))
        
        
    def setWorkingDir(self):
        self.workingDirectory=filedialog.askdirectory(title='Select Working Directory')
        try:
            os.chdir(self.workingDirectory)
            self.printToLog('Working Directory changed to: '+self.workingDirectory)
        except OSError:
            self.printToLog('No New Working Directory Selected')
        

        
    def getReads(self):
        #self.printToLog("checkbutton value = "+str(self.app.selectDirVar.get()))
        if self.app.selectDirVar.get() == 1:
            self.processLargeFile = True
        else:
            self.processLargeFile = False
        if self.folderMode==False and self.processLargeFile == False:
            if self.platformType == "Windows":
                temp=filedialog.askopenfilenames(title='Select Reads',initialdir=self.workingDirectory,filetypes=(("FastQ Files", "*.fq;*.fastq;*.txt"),("Fasta Files", "*.fa;*.fasta;*.txt"),("All files", "*.*")))
            else:
                temp=filedialog.askopenfilenames(title='Select Reads',initialdir=self.workingDirectory)

            
            if temp:
                if len(temp)>2:
                    self.printToLog('Max selection of 2 reads permitted (forward/reverse)')
                    self.printToLog('If you would like to run more, please use directory mode')
                    return 'break'
                else:
                    self.readLocation=temp
                    if len(self.readLocation) ==2 :
                        if (self.readLocation[0].endswith('.fq') | self.readLocation[0].endswith('.fastq')):
                            if (self.readLocation[1].endswith('.fq') | self.readLocation[1].endswith('.fastq')):
                                self.paired=True
                            else:
                                self.printToLog('Cannot mix file types in mate pair (Fasta with FastQ)')
                                return 'break'
                        elif(self.readLocation[0].endswith('.fa') | self.readLocation[0].endswith('.fasta')):
                            self.printToLog('Cannot have mate pair of Fasta file type')
                            return 'break'
                    else:
                        self.paired=False
                    self.readBaseName=os.path.basename(self.readLocation[0])
                    self.readBaseName=os.path.splitext(self.readBaseName)[0]
                    self.readBaseName=self.readBaseName[0:len(self.readBaseName)-1]
                    self.readDir=os.path.dirname(self.readLocation[0])
                    self.printToLog('Reads Selected:')
                    for i in range (0,len(self.readLocation)):
                        self.printToLog(str(self.readLocation[i]))
                    return self.readLocation
            else:
                self.printToLog('No reads selected')
                return 'break'
        else:
            if (self.getProcessingLargeFileStatus() == False) and (self.processLargeFileIndex == 0):
                temp=filedialog.askdirectory(title='Select Reads Directory',initialdir=self.workingDirectory+'/Reads')
                if temp:
                    self.readDir=temp
                    self.processLargeFileDir = temp
                    self.printToLog('Directory Selected:\n'+self.readDir)
                    self.paired=True
                    try:
                        filename=self.processLargeFileDir+"/"+self.splitInfoFilename
                        #self.printToLog("split info file = "+filename)
                        try:
                            splitInfoFile=open(filename, "r")
                        except:
                            self.printToLog('\n------- Error!!! ----------------------------')
                            self.printToLog('can not locate split information file '+filename+".")
                            self.printToLog('---------------------------------------------\n')
                            return 'break'
                        fileCount=splitInfoFile.readline().replace("\n","")
                        self.splitFileCount=int(fileCount)
                        file1=splitInfoFile.readline().replace("\n","")
                        file2=splitInfoFile.readline().replace("\n","")
                        self.splitFile1=self.processLargeFileDir+"/"+file1
                        self.splitFile2=self.processLargeFileDir+"/"+file2
                        # read selected files list if any
                        try:
                            selectedFiles = splitInfoFile.readline().replace("\n","")
                            selectedFiles = selectedFiles.replace(" ", "")
                            if selectedFiles =="": #blank line
                                self.splitSelectedFiles=[]
                            else:
                                self.splitSelectedFiles = selectedFiles.split(",")   
                        except:
                            self.splitSelectedFiles=[]
                        self.printToLog('FSS == self.splitSelectedFilesUI len'+str(len(self.splitSelectedFilesUI)))
                        if len(self.splitSelectedFilesUI) > 0:
                            self.splitSelectedFiles = self.splitSelectedFilesUI
                        splitInfoFile.close()
                        self.processLargeFileIndex += 1   #increase index
                        # check selected files specified
                        if len(self.splitSelectedFiles) == 0: # no selected files specified
                            file1=self.processLargeFileDir+"/"+file1.replace("%", str(self.processLargeFileIndex))
                            file2=self.processLargeFileDir+"/"+file2.replace("%", str(self.processLargeFileIndex))
                        else: # selected files specified
                            while (self.processLargeFileIndex <= self.splitFileCount) and (str(self.processLargeFileIndex) not in self.splitSelectedFiles):
                                self.processLargeFileIndex += 1
                            if self.processLargeFileIndex > self.splitFileCount:
                                self.processLargeFileIndex = 1
                            file1=self.processLargeFileDir+"/"+file1.replace("%", str(self.processLargeFileIndex))
                            file2=self.processLargeFileDir+"/"+file2.replace("%", str(self.processLargeFileIndex))
                        #self.printToLog("====== "+file1+"  "+file2)
                        if not(os.path.exists(file1)):
                            self.readLocation=''
                            self.printToLog("File not found: "+file1)
                            self.processLargeFileIndex = -1
                            self.processLargeFileCheckGTF=None
                            self.splitSelectedFiles=[]
                            self.splitSelectedFilesUI=[]
                            self.changeProcessingLargeFileCompleteStatus(True)
                            return 'break'
                        if not(os.path.exists(file2)):
                            self.readLocation=''
                            self.printToLog("File not found: "+file2)
                            self.processLargeFileIndex = -1
                            self.processLargeFileCheckGTF=None
                            self.splitSelectedFiles=[]
                            self.splitSelectedFilesUI=[]
                            self.changeProcessingLargeFileCompleteStatus(True)
                            return 'break'
                        self.readLocation=[file1, file2]
                        self.readBaseName=os.path.basename(self.readLocation[0])
                        self.readBaseName=os.path.splitext(self.readBaseName)[0]
                        self.readBaseName=self.readBaseName[0:len(self.readBaseName)-1]
                        self.readDir=os.path.dirname(self.readLocation[0])
                        self.updateRunDirectory()
                        self.inProcesingDirLoop = True
                        return self.readDir
                    except:
                        return 'break'
                else:
                    self.printToLog('No directory selected')
                    return 'break'
            else:
                self.printToLog('find next file to be processed........')
                if self.processLargeFileIndex == -1:
                    self.splitFile1=''
                    self.splitFile2=''
                    self.printToLog('all files processed')
                    self.changeProcessingLargeFileCompleteStatus(True)
                    return 'break'
                self.processLargeFileIndex += 1   #increase index
                if self.processLargeFileIndex > self.splitFileCount:
                    self.readLocation=''
                    self.processLargeFileIndex = -1
                    self.processLargeFileCheckGTF=None
                    self.splitSelectedFiles=[]
                    self.splitSelectedFilesUI=[]
                    self.changeProcessingLargeFileCompleteStatus(True)
                    return 'break'
                
                # check selected files specified
                if len(self.splitSelectedFiles) == 0: # no selected files specifies
                    file1=self.splitFile1.replace("%", str(self.processLargeFileIndex))
                    file2=self.splitFile2.replace("%", str(self.processLargeFileIndex))
                else: # selected files specified
                    while (self.processLargeFileIndex <= self.splitFileCount) and (str(self.processLargeFileIndex) not in self.splitSelectedFiles):
                        self.processLargeFileIndex += 1
                    if self.processLargeFileIndex > self.splitFileCount:
                        self.readLocation=''
                        self.processLargeFileIndex = -1
                        self.processLargeFileCheckGTF=None
                        self.splitSelectedFiles=[]
                        self.splitSelectedFilesUI=[]
                        self.changeProcessingLargeFileCompleteStatus(True)
                        return 'break'
                    file1=self.splitFile1.replace("%", str(self.processLargeFileIndex))
                    file2=self.splitFile2.replace("%", str(self.processLargeFileIndex))
                    
                self.printToLog('get next file........')
                if not(os.path.exists(file1)):
                    self.readLocation=''
                    self.printToLog("File not found: "+file1)
                    self.processLargeFileIndex = -1
                    self.processLargeFileCheckGTF=None
                    self.splitSelectedFiles=[]
                    self.splitSelectedFilesUI=[]
                    self.changeProcessingLargeFileCompleteStatus(True)
                    return 'break'
                if not(os.path.exists(file2)):
                    self.readLocation=''
                    self.printToLog("File not found: "+file2)
                    self.processLargeFileIndex = -1
                    self.processLargeFileCheckGTF=None
                    self.splitSelectedFiles=[]
                    self.splitSelectedFilesUI=[]
                    self.changeProcessingLargeFileCompleteStatus(True)
                    return 'break'
                self.readBaseName=os.path.basename(self.readLocation[0])
                self.readBaseName=os.path.splitext(self.readBaseName)[0]
                self.readBaseName=self.readBaseName[0:len(self.readBaseName)-1]
                self.readDir=os.path.dirname(self.readLocation[0])
                self.readLocation=[file1, file2]
                self.changeProcessingLargeFileStatus(True)
                return self.readDir

    def checkBowtie2FolderExists(self):
        self.printToLog('Verifying Bowtie2 Integrity.....................')
        GO=False
        if self.platformType == "Windows":
            if os.path.exists(self.Bowtie2Folder) and os.path.isfile(self.Bowtie2Folder+'/bowtie2.bat'):
                GO=True
        else:
            if os.path.exists(self.Bowtie2Folder) and os.path.isfile(self.Bowtie2Folder+'/bowtie2'):
                GO=True
        if GO:
            self.printToLog('Bowtie2 Integrity Verified')
        else:
            self.printToLog('Please check Bowtie2 Folder is chosen correctly and is not corrupted')
        return GO
    
    def checkBowtieIndexs(self):
        self.printToLog('Checking Bowtie2 Indexes...')
        stopv=False
        stoph=False
        self.viralPrefix=str(self.app.vp.get(1.0,'end-1c'))
        self.HostPrefix=str(self.app.hp.get(1.0,'end-1c'))
        for i in range (1,7):
            if i<5:
                name=self.ViralRefFolder+'/'+self.viralPrefix+'.'+str(i)+'.bt2'
                name2=self.HostRefFolder+'/'+self.HostPrefix+'.'+str(i)+'.bt2'
            else:
                name=self.ViralRefFolder+'/'+self.viralPrefix+'.rev.'+str(i-4)+'.bt2'
                name2=self.HostRefFolder+'/'+self.HostPrefix+'.rev.'+str(i-4)+'.bt2'
            if not os.path.isfile(name):
                stopv=True
            if not os.path.isfile(name2):
                stoph=True
        if stopv:
            self.printToLog('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            self.printToLog('Bowtie 2 index for viral not found, please create Bowtie2 index or specify correct path')
            self.printToLog('Checked in: '+self.ViralRefFolder+'\nFor: '+self.viralPrefix)
            self.printToLog('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        if stoph:
            self.printToLog('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            self.printToLog('Bowtie 2 index for host not found, please create Bowtie2 index or specify correct path')
            self.printToLog('Checked in: '+self.HostRefFolder+'\nFor: '+self.HostPrefix)
            self.printToLog('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        if stoph==False and stopv==False:
            self.printToLog('Bowtie2 indexes existences verified')
            return True
        else:
            return False
                
    
 

                
    def askbuildBowtie2ViralIndex(self,prefix):
        if self.ViralRefFa!='None Selected':
            self.printToLog('Building Index for file:\n'+self.ViralRefFa)
            self.printToLog("In:\n"+self.ViralRefFolder)  
            self.printToLog('As:\n'+prefix)
            self.printToLog('Please confirm build...(Y/N) button')
            self.app.yesbttn.configure(command=lambda: self.buildBowtie2ViralIndex(prefix))
            self.app.nobttn.configure(command=self.cancel_buildBowtie2Index)
        else:
            self.printToLog('No selected viral reference fasta file to build on')

    def buildBowtie2ViralIndex(self,prefix):    
        def target():
            try:
                self.app.update()
                os.chdir(str(self.ViralRefFolder))
                if self.platformType == "Windows":                
                    temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2-build.bat',self.ViralRefFa,prefix],stdout=subprocess.PIPE,bufsize=1)#,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
                else:
                    temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2-build',self.ViralRefFa,prefix],stdout=subprocess.PIPE,bufsize=1)#,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)    
                self.app.update()
                for line in iter(temp.stdout.readline,b''):
                    self.printToLog(line.decode('utf-8'))
                temp.stdout.close()
                self.printToLog('======================================================================\nBuild Complete')
            except:
                err=traceback.format_exc()
                status=-1
        thread=Thread(target=target)
        thread.start()
        os.chdir(self.workingDirectory)
        self.app.resetYesNoButtons()
        
    def cancel_buildBowtie2Index(self):
        self.printToLog('Bowtie2 Index build canceled!')
        self.app.resetYesNoButtons()
        
    def saveSettings(self):
        try: 
            out=open(self.workingDirectory+'/config.txt','w')
            go=True
        except:
            self.printToLog('Cannot save data, file currently in use')
            go=False
        if go:
            settings=[]
            settings.append('BT2D,'+self.Bowtie2Folder)
            settings.append('VRFa,'+self.ViralRefFa)
            settings.append('VID,'+self.ViralRefFolder)
            settings.append('VP,'+self.viralPrefix)
            settings.append('HRFa,'+self.HostRefFa)
            settings.append('HID,'+self.HostRefFolder)
            settings.append('HP,'+self.HostPrefix)
            settings.append('GTF,'+str(self.HostGTF))
            settings.append('OutDir,'+str(self.outDir))
            settings.append('TC,'+str(self.coreOption))
            settings.append('MinLen,'+str(self.CSminLen))
            settings.append('Na,'+str(self.Na))
            settings.append('Auto,'+str(self.autoMode))
            settings.append('Load,'+str(self.loadAlignments))
            settings.append('DoGene,'+str(self.doGTF))
            settings.append('Trim5,'+str(self.Trim))
            settings.append('Trim3,'+str(self.Trim2))
            settings.append('OverlapTmMax,'+str(self.otempFilter))
            settings.append('OverlapLenMax,'+str(self.olenFilter))
            settings.append('HostTmMin,'+str(self.htempFilter))
            settings.append('HostLenMin,'+str(self.hlenFilter))
            settings.append('ViralTmMin,'+str(self.vtempFilter))
            settings.append('ViralLenMin,'+str(self.vlenFilter))
            settings.append('FilterBasic,'+str(self.useBasic))
            settings.append('GeneThresh,'+str(self.distanceThreshold))
            settings.append('Microhomology_sel,'+str(self.microhomologystringency_selection))
            if (self.microhomologystringency==""):
                settings.append('Microhomology,'+str(self.microhomologystringency_avg))
            else:
                settings.append('Microhomology,'+str(self.microhomologystringency))
            if (self.microhomologystringency1==""):
                settings.append('Microhomology,'+str(self.microhomologystringency1_avg))
            else:
                settings.append('Microhomology,'+str(self.microhomologystringency1))
            for item in settings:
                out.write(item+'\n')
            out.close()
        self.printToLog('Settings saved successfully!')

    def splitLargeFile(self):
        self.file1=None
        self.file2=None
        def printInfo(info):
            self.splitLargeFileOutText.configure(state=NORMAL)
            self.splitLargeFileOutText.insert(END, info+"\n")
            self.splitLargeFileOutText.see(END)
            self.splitLargeFileOutText.configure(state=DISABLED)
            self.splitLargeFileOutText.update()
            return

        def startSplitFiles():
            start_time=time.time()
            file1=self.file1
            file2=self.file2
            def batch_iterator(iterator, batch_size):
                entry = True  # Make sure we loop once
                while entry:
                    batch = []
                    count = 0
                    while len(batch) < batch_size:
                        try:
                            entry = iterator.__next__()
                            count +=1
                            if (count % 200000) == 0:
                                printInfo("reads processed "+str(count)+" of "+str(batch_size))
                                #time.sleep(1)
                        except StopIteration:
                            entry = None
                        if entry is None:
                            # End of file
                            printInfo("reads processed "+str(count)+" of "+str(batch_size))
                            break
                        batch.append(entry)
                    if batch:
                        yield batch
            self.splitLargeFileBttn2.configure(state=DISABLED) # prevent spliting again in the middle of process
            printInfo("Splitting Files...")            
            printInfo(file1)
            if (file2 != None):
                printInfo(file2)
        
            readsCount = 1000000
            #readsCount = 40000
            dirName, pathName=os.path.split(file1)
            if pathName.lower().endswith("fastq"):
                newDirName=pathName[:-7]+"_split"
                a=str(pathName[:-7])
                a0=str(pathName[:-7])
            else:
                newDirName=pathName[:-4]+"_split"
                a=str(pathName[:-4])
                a0=str(pathName[:-4])

            newDirNameFull=dirName+"/"+newDirName

            try:
                os.mkdir(newDirNameFull)
            except FileExistsError:
                if os.path.exists(newDirNameFull):
                    printInfo("path exist already")
                else:
                    printInfo("create failed")
                    return

            """    
            record_iter = SeqIO.parse(open(file1),"fastq")
            printInfo('Reading '+file1)
            file_count1=0
            for i, batch in enumerate(batch_iterator(record_iter, readsCount)):
                b="_"+str(i+1)
                c='_1.fastq'
                filename=newDirNameFull+"/"+a+b+c
                printInfo(filename)
                handle = open(filename, "w")
                count = SeqIO.write(batch, handle, "fastq")
                handle.close()
                file_count1 +=1
                printInfo("Wrote %i records to %s" % (count, filename))
            b="_%"
            filename1_format=a0+b+c
            
            file_count2=0
            if (file2 != None):
                printInfo(file2)
                record_iter = SeqIO.parse(file2,"fastq")
                for i, batch in enumerate(batch_iterator(record_iter, readsCount)):
                    b="_"+str(i+1)
                    c='_2.fastq'
                    filename=newDirNameFull+"/"+a+b+c
                    printInfo(filename)
                    handle = open(filename, "w")
                    count = SeqIO.write(batch, handle, "fastq")
                    handle.close()
                    file_count2 +=1
                    printInfo("Wrote %i records to %s" % (count, filename))
            self.splitLargeFileBttn2.configure(state=DISABLED)
            end_time=time.time()
            b="_%"
            filename2_format=a0+b+c
            """
            showReadsCount = readsCount / 5 
            printInfo('Reading '+file1)
            file_count1=1
            b="_"+str(file_count1)
            c='_1.fastq'
            filename=newDirNameFull+"/"+a+b+c
            printInfo("creating file: "+filename)
            handle = open(filename, "w")
            cnt = 0
            for (title, sequence, quality) in FastqGeneralIterator(open(file1)) :
                handle.write("@"+title+"\n")
                handle.write(sequence+"\n+\n")
                handle.write(quality+"\n")
                cnt +=1
                if (cnt % showReadsCount) == 0: # display the progres
                    printInfo("reads processed "+str(cnt)+" of "+str(readsCount))
                if (cnt % readsCount) == 0: # one small file done
                    handle.close()
                    printInfo("Wrote %i records to %s" % (cnt, filename))
                    file_count1 += 1
                    b="_"+str(file_count1)
                    c='_1.fastq'
                    filename=newDirNameFull+"/"+a+b+c
                    printInfo("creating file: "+filename)
                    handle = open(filename, "w")
                    cnt = 0
            handle.close()
            printInfo("Wrote %i records to %s" % (cnt, filename))
            b="_%"
            filename1_format=a0+b+c

            printInfo('Reading '+file2)
            file_count2=1
            b="_"+str(file_count2)
            c='_2.fastq'
            filename=newDirNameFull+"/"+a+b+c
            printInfo("creating file: "+filename)
            handle = open(filename, "w")
            cnt = 0
            for (title, sequence, quality) in FastqGeneralIterator(open(file2)) :
                handle.write("@"+title+"\n")
                handle.write(sequence+"\n+\n")
                handle.write(quality+"\n")
                cnt +=1
                if (cnt % showReadsCount) == 0: # display the progres
                    printInfo("reads processed "+str(cnt)+" of "+str(readsCount))
                if (cnt % readsCount) == 0: # one small file done
                    handle.close()
                    printInfo("Wrote %i records to %s" % (cnt, filename))
                    file_count2 += 1
                    b="_"+str(file_count2)
                    c='_2.fastq'
                    filename=newDirNameFull+"/"+a+b+c
                    printInfo("creating file: "+filename)
                    handle = open(filename, "w")
                    cnt = 0
            handle.close()
            printInfo("Wrote %i records to %s" % (cnt, filename))
            b="_%"
            filename2_format=a0+b+c


            end_time=time.time()
                        
            printInfo("=====================================================================")
            printInfo("filename1 format = "+filename1_format)
            printInfo("filename2 format = "+filename2_format)
            out=open(newDirNameFull+"/"+self.splitInfoFilename,"w")
            out.write(str(file_count1)+"\n")
            out.write(filename1_format+"\n")
            out.write(filename2_format+"\n")
            out.close()
            printInfo(" new folder '"+newDirNameFull+"' with "+str(file_count1+file_count2)+" files created")
            printInfo(" time used = "+str(end_time-start_time)+" seconds")
            printInfo("========================    Split Complete   ========================")
            printInfo("=====================================================================")    
            return
        
        def readLargeFilles():
            if self.platformType == "Windows":
                temp=filedialog.askopenfilenames(title='Select Reads',initialdir=self.workingDirectory,filetypes=(("FastQ Files", "*.fq;*.fastq;*.txt"),("Fasta Files", "*.fa;*.fasta;*.txt"),("All files", "*.*")))
            else:
                temp=filedialog.askopenfilenames(title='Select Reads',initialdir=self.workingDirectory)
            if len(temp)>2:
                printInfo('Max selection of 2 reads permitted (forward/reverse)')
                return

            if (not temp[0].lower().endswith("fastq")) and (not temp[0].lower().endswith("fq")):
                printInfo(temp[0]+' not in *.fastq format')
                return

            if len(temp)==2:
                if (not temp[1].lower().endswith("fastq")) and (not temp[1].lower().endswith("fq")):
                    printInfo(temp[1]+' not in *.fastq format')
                    return
                    
            files=temp[0]
            self.file1=temp[0]
            if len(temp)==2:
                files=files+"\n"+temp[1]
                self.file2=temp[1]
            self.app.changeLabel(self.splitLargeFileFastqBox, files)
            self.splitLargeFileBttn2.configure(state=NORMAL)
            return

        def closeSplitLargeFileWin():
            self.splitLargeFileWin.destroy()
            return
        self.printToLog('splitLargeFile!')
        self.splitLargeFileWin=tki.Toplevel(self.app.master)
        self.splitLargeFileWin.title('Split large paired fastq files')
        self.splitLargeFileWin.geometry('620x620')   
        self.splitLargeFileWin.protocol('WM_DELETE_WINDOW',closeSplitLargeFileWin)
        self.labelSplit=Text(self.splitLargeFileWin,height=1,width=20,state=DISABLED,relief=FLAT,background=self.app.defaultbg)
        self.app.changeLabel(self.labelSplit,'Fastq File(s):')
        self.labelSplit.grid(row=0,column=0,columnspan=1,sticky=W)
        
        self.splitLargeFileFrame1=Frame(self.splitLargeFileWin)
        self.splitLargeFileFrame1.grid(row=1,column=0,sticky='nw',columnspan=2)         
        self.splitLargeFileBttn1=Button(self.splitLargeFileFrame1,text = '. . .', command=readLargeFilles)
        self.splitLargeFileBttn1.grid(row=0,column=0,sticky='',ipadx=3,ipady=5)
        self.splitLargeFileFastqBox=Text(self.splitLargeFileFrame1,height=2,width=70,state=DISABLED,background=self.app.defaultbg)
        self.splitLargeFileFastqBox.grid(row=0,column=1,sticky=W,columnspan=1)
        self.splitLargeFileScrollbarFile1=ttk.Scrollbar(self.splitLargeFileFrame1,command=self.splitLargeFileFastqBox.yview)
        self.splitLargeFileFastqBox.config(yscrollcommand=self.splitLargeFileScrollbarFile1.set)
        self.splitLargeFileScrollbarFile1.grid(row=0,column=2,sticky='ns')

        self.splitLargeFileBttn2=Button(self.splitLargeFileFrame1,text = ' Split Files ', state = 'disabled', command=startSplitFiles)
        self.splitLargeFileBttn2.grid(row=1,column=1,sticky='',ipadx=3,ipady=5)

        self.splitLargeFileOutText=Text(self.splitLargeFileFrame1,height=30,width=70,state=DISABLED, background=self.app.defaultbg)
        self.splitLargeFileOutText.grid(row=2,column=0,columnspan=2,sticky='wens')
        self.splitLargeFileScrollbarFile2=ttk.Scrollbar(self.splitLargeFileFrame1,command=self.splitLargeFileOutText.yview)
        self.splitLargeFileOutText.config(yscrollcommand=self.splitLargeFileScrollbarFile2.set)
        self.splitLargeFileScrollbarFile2.grid(row=2,column=2,sticky='ns')
        """
        self.splitLargeFileScrollbarFile3=ttk.Scrollbar(self.splitLargeFileFrame1,command=self.splitLargeFileOutText.xview, orient=HORIZONTAL)
        self.splitLargeFileOutText.config(yscrollcommand=self.splitLargeFileScrollbarFile3.set)
        self.splitLargeFileScrollbarFile3.grid(row=3,column=0,columnspan=2,sticky='ew')
        """
        self.splitLargeFileWin.grab_set()
        return
    
    def loadSettings(self):
        if self.platformType == "Windows":
            temp=filedialog.askopenfilename(title='Select Config File',initialdir=self.workingDirectory,filetypes=(("text", "*.txt"),("All files", "*.*")))
        else:
            temp=filedialog.askopenfilename(title='Select Config File',initialdir=self.workingDirectory)  
        if temp:
            out=open(temp,mode='r')
            for line in out:
                a=line.replace("\n","").split(',')
                if a[0]=='BT2D':
                    self.Bowtie2Folder=a[1]
                    self.app.changeLabel(self.app.bt2Text,a[1])
                elif a[0]=='VRFa':
                    self.ViralRefFa=a[1]
                    self.app.changeLabel(self.app.vrf,a[1])
                elif a[0]=='VID':
                    self.ViralRefFolder=a[1]
                    self.app.changeLabel(self.app.vdd,a[1])
                elif a[0]=='VP':
                    self.viralPrefix=a[1]
                    self.app.changeLabelFree(self.app.vp,a[1])
                elif a[0]=='HRFa':
                    self.HostRefFa=a[1]
                    self.app.changeLabel(self.app.hrf,a[1])
                elif a[0]=='HID':
                    self.HostRefFolder=a[1]
                    self.app.changeLabel(self.app.hdd,a[1])
                elif a[0]=='HP':
                    self.HostPrefix=a[1]
                    self.app.changeLabelFree(self.app.hp,a[1])
                elif a[0]=='GTF':
                    self.HostGTF=a[1]
                    self.app.changeLabel(self.app.hgtf,a[1])
                elif a[0]=='OD':
                    self.workingDirectory=a[1]
                    self.app.changeLabel(self.app.out,a[1])
                elif a[0]=='TC':
                    self.coreOption=int(a[1])
                    self.app.cpuIntVar.set(int(a[1]))
                elif a[0]=='MinLen':
                    self.CSminLen=int(a[1])
                    self.app.changeLabelFree(self.app.clipText,a[1])
                elif a[0]=='Na':
                    self.Na=int(a[1])
                    self.app.changeLabelFree(self.app.saltText,a[1])
                elif a[0]=='Auto':
                    self.autoMode=bool(a[1])
                    if bool(a[1]):
                        self.app.autoModeVar.set(1)
                elif a[0]=='Load':
                    self.loadAlignments=bool(a[1])
                    if bool(a[1]):
                        self.app.loadAlignmentBoolVar.set(1)
                """
                elif a[0]=='Microhomology_sel':
                    self.microhomologystringency_selection = int(a[1])
                elif a[0]=='Microhomology':
                    self.microhomologystringency = int(a[1])
                elif a[0]=='Microhomology1':
                    self.microhomologystringency1 = int(a[1])
                """  
            out.close()
        self.printToLog('Config file loaded')
        self.printDefaults()
        
        
        
    def askbuildBowtie2HostIndex(self,prefix): 
        if self.HostRefFa!='None Selected':
            self.printToLog('Building Index for file:\n'+self.HostRefFa)
            self.printToLog("In:\n"+self.HostRefFolder)  
            self.printToLog('As:\n'+prefix)
            self.printToLog('Please confirm build...(Y/N) button')
            self.app.yesbttn.configure(command=lambda: self.buildBowtie2HostIndex(prefix))
            self.app.nobttn.configure(command=self.cancel_buildBowtie2Index)
        else:
            self.printToLog('No selected host reference fasta file to build on')
            
    def buildBowtie2HostIndex(self,prefix): 
        def target():
            try:
                self.app.update()
                os.chdir(str(self.HostRefFolder))
                if self.platformType == "Windows":
                    temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2-build.bat',self.HostRefFa,prefix],stdout=subprocess.PIPE,bufsize=1)#,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
                else:
                    temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2-build',self.HostRefFa,prefix],stdout=subprocess.PIPE,bufsize=1)#,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)                
                self.app.update()
                for line in iter(temp.stdout.readline,b''):
                    self.printToLog(line.decode('utf-8'))
                temp.stdout.close()
                self.printToLog('======================================================================\nBuild Complete')
            except:
                err=traceback.format_exc()
                status=-1     
                self.printToLog('Error Occured\n'+err)
        thread=Thread(target=target)
        thread.start()
        os.chdir(self.workingDirectory)
        self.app.resetYesNoButtons()

    def run1(self, x):
        self.printToLog('Starting read preprocessing with button - 3')
        return

    def changeProcessingLargeFileStatus(self, status):
        self.processingLargeFileLock.acquire()
        self.processingLargeFile = status
        self.processingLargeFileLock.release()
        return

    def getProcessingLargeFileStatus(self):
        self.processingLargeFileLock.acquire()
        status = self.processingLargeFile
        self.processingLargeFileLock.release()
        return status

    def changeProcessingLargeFileCompleteStatus(self, status):
        self.processingLargeFileCompleteLock.acquire()
        self.processingLargeFileComplete=status
        self.processingLargeFileCompleteLock.release()
        return

    def getProcessingLargeFileCompleteStatus(self):
        self.processingLargeFileCompleteLock.acquire()
        status=self.processingLargeFileComplete
        self.processingLargeFileCompleteLock.release()
        return status
    
    def changeProcessingSaveDataStatus(self, status):
        self.processingSaveDataLock.acquire()
        self.processingSaveData=status
        self.processingSaveDataLock.release()
        return
    
    def getProcessingSaveDataStatus(self):
        self.processingSaveDataLock.acquire()
        status=self.processingSaveData
        self.processingSaveDataLock.release()
        return status
    
    def run(self):
        def target():
            start=time.time()
            if self.doReadClean:
                self.printToLog('Starting read preprocessing')
                try:
                    self.preprocessReadLength=int(self.app.readminlength.get(1.0,'end-1c'))
                    self.preprocessUnknownPercent=int(self.app.maxunknown.get(1.0,'end-1c'))
                except:
                    pass
                for read in self.readLocation:
                    self.sequence_cleaner(read,self.preprocessReadLength,self.preprocessUnknownPercent)
                self.event.wait()
                self.event.clear()
            """
            else:
                self.printToLog('Skipping read preprocessing')
            """
            if self.doGTF:
                self.loadGTF()
                self.event.wait()
                self.event.clear()
            else:
                self.GO=True
            if self.GO:
                self.askviralAlignment()
                self.event.wait()
                self.event.clear()
                if self.GO:
                    self.getCS()
                    if self.GO:
                        self.askHostAlignment()
                        self.event.wait()
                        self.event.clear()
                        if self.GO:
                            self.printToLog('Building Internal Tables...')
                            self.mapper()
                            self.printToLog('Tables Filled!')
                            stop=time.time()
                            runtime=stop-start
                            self.printToLog('Check for unique reads...')
                            self.countUnique()
                            self.printToLog('======================================================================\nRun Completed in '+str(runtime)+' seconds')
                            self.GO=False
                            # check for mictohomology & similarity
                            self.printToLog('\nCheck for microhomology & high similarity reads...')
                            self.app.listRedraw()
                            self.app.selectDirButton.configure(state=NORMAL)
                            self.app.bttnStart.configure(state=NORMAL)                            
                            if self.processLargeFile:
                                self.changeProcessingSaveDataStatus(True)
                                self.saveData()
                                self.outDir='None Selected'
                                self.changeProcessingLargeFileStatus(False)
        try:
            self.distanceThreshold=int(self.app.DistanceText.get(1.0,'end-1c'))
        except:
            self.distanceThreshold=10000
        if not self.processLargeFile:
            temp=self.checkReads()
            stop=False
            if temp==0:
                stop=True
                self.printToLog('Aborting: No Read(s) Selected')
            if stop==False:
                temp=self.checkBowtie2FolderExists()
                if temp:
                    temp=self.checkBowtieIndexs()
                    if temp:
                        self.printToLog('Initiating Run\n======================================================================')
                        self.app.selectDirButton.configure(state=DISABLED)
                        self.app.bttnStart.configure(state=DISABLED)
                        self.event.clear()
                        thread=Thread(target=target)
                        thread.start()
        else:
            self.autoMode = True
            if self.getProcessingLargeFileStatus()==False:
                self.changeProcessingLargeFileCompleteStatus(False)
                self.printToLog('Processing multiple files...........')
                if (self.processLargeFileFirstRead):
                    self.splitStartTime=time.time()
                    self.app.selectDirButton.configure(state=DISABLED)
                    self.app.bttnStart.configure(state=DISABLED)
                    temp=self.checkReads()
                    self.processLargeFileFirstRead=False
                else:
                    if self.getProcessingSaveDataStatus() == True:
                        self.printToLog('waiting for save complete...')
                        tc = threading.Timer(5, self.run)
                        tc.start()
                        return
                    temp=self.getReads()
                    if temp=='break':
                        self.printToLog('No more files to be processed......')
                        self.changeProcessingLargeFileStatus(False)
                        self.changeProcessingLargeFileCompleteStatus(True)
                        self.app.selectDirButton.configure(state=NORMAL)
                        self.app.bttnStart.configure(state=NORMAL)
                        if self.processLargeFileInterrupt == None:
                            self.saveMultipleFinalReads()
                        self.splitEndTime=time.time()
                        self.printToLog('\n======================================================================\n')
                        self.printToLog('Run Multiple Files Completed in '+str(self.splitEndTime - self.splitStartTime)+' seconds \n')
                        if (self.savedReadsLogFilename != ""):
                            try: 
                                self.printToLog('Final log will be saved to file:\n'+self.savedReadsLogFilename)
                                out2=open(self.savedReadsLogFilename, "w")
                                for line in self.logger:
                                    out2.write(line+'\n')
                                out2.close()
                                self.savedReadsLogFilename = ""
                            except:
                                self.printToLog('Cannot save log to '+self.savedReadsLogFilename)
                        return
                    temp=self.checkReads()                    
                stop=False
                if temp==0:
                    stop=True
                    self.printToLog('Aborting: No Directory Selected')
                    self.changeProcessingLargeFileStatus(False)
                    self.changeProcessingLargeFileCompleteStatus(True)
                    self.app.selectDirButton.configure(state=NORMAL)
                    return
                self.printToLog(str(temp)+" == "+self.processLargeFileDir)
                if stop==False:
                    temp=self.checkBowtie2FolderExists()
                    if temp:
                        temp=self.checkBowtieIndexs()
                        if temp:
                            self.printToLog('Initiating Run\n======================================================================')
                            if (len(self.readLocation) < 2):
                                self.printToLog('Can not find paired input files!')
                                stop=True
                                return
                            self.printToLog('File1 = '+self.readLocation[0])
                            self.printToLog('File2 = '+self.readLocation[1])
                            self.readBaseName=os.path.basename(self.readLocation[0])
                            self.readBaseName=os.path.splitext(self.readBaseName)[0]
                            self.readBaseName=self.readBaseName[0:len(self.readBaseName)-1]
                            self.readDir=os.path.dirname(self.readLocation[0])
                            self.outDir=os.path.dirname(self.readLocation[0])
                            self.runDirectory=os.path.dirname(self.readLocation[0])
                            #self.updateRunDirectory()
                            self.event.clear()
                            thread=Thread(target=target)
                            thread.start()
                self.printToLog('start threading ...........')
                self.changeProcessingLargeFileStatus(True)

            if self.getProcessingLargeFileCompleteStatus() == False:
                if self.GO:
                    self.printToLog('processing multiple files...')
                tc = threading.Timer(15, self.run)
                tc.start()
                return
            else:    
                self.printToLog('processing multiple files complete...')
                return
            
                
        
    def checkReads(self):
        if self.processLargeFile == True:
            if self.processLargeFileDir == "":
                return 0
            else:
                return 1
        if not self.folderMode:
            return len(self.readLocation)
    
    def addPaths(self):
        if self.platformType == "Windows":
            try:
                os.environ['PATH']+=os.pathsep+self.workingDirectory+'/Support/Python'
                os.environ['PATH']+=os.pathsep+self.workingDirectory+'/Support/Perl'
                os.environ['PATH']+=os.pathsep+self.workingDirectory+'/Support/Perl/perl'
                os.environ['PATH']+=os.pathsep+self.workingDirectory+'/Support/Perl/perl/bin'
                os.environ['PATH']+=os.pathsep+self.workingDirectory+'/Support/Perl/perl/lib'

                    
            except KeyError:
                os.environ['PATH']=os.pathsep+self.workingDirectory+'/Support/Python'
                os.environ['PATH']+=os.pathsep+self.workingDirectory+'/Support/Perl'
                os.environ['PATH']+=os.pathsep+self.workingDirectory+'/Support/Perl/perl'
                os.environ['PATH']+=os.pathsep+self.workingDirectory+'/Support/Perl/perl/bin'
                os.environ['PATH']+=os.pathsep+self.workingDirectory+'/Support/Perl/perl/lib'

        else:
            try: 
                os.environ['PATH']+=os.pathsep+self.workingDirectory+'/tcl/tk8.5.2/Python'
            except KeyError:
                os.environ['PATH']=os.pathsep+self.workingDirectory+'/tcl/tk8.5.2/Python'
            
    def askviralAlignment(self):
        stop=False
        forward=-1
        reverse=-1
        mode=None
        count=0
        if not self.folderMode:      #if in single alignment mode (using only one read or one read pair, not batch to folder)
            if len(self.readLocation) ==2 :  #if you have exactly two reads selected
                mode='paired'
                mode2='fq'
                for i in range (0,2):
                    if (self.readLocation[i].endswith('.fq') | self.readLocation[i].endswith('.fastq')):  #check if each read ends with .fq or .fastq
                        if self.readLocation[i].endswith('1.fq')|self.readLocation[i].endswith('1.fastq'):
                            forward=i
                            count+=1
                    else:
                        stop=True
                        self.printToLog("Read "+self.readLocation[i]+' does not have correct file extenstion (expecting .fq or .fastq)') #print error if not ending in .fq or .fastq
                if count!=1:  #makes sure a forward read is identified 
                    stop=True
                    self.printToLog('Could not designate forward read, naming convention set to\nForward=name1.fq/.fastq Reverse=name2.fq/.fastq') 
                else:
                    if forward==1:
                        reverse=0
                    else:
                        reverse=1
                    if self.doReadClean: #change for final code (delete)
                        temp=self.readLocation[forward].split('/')
                        temp[-1]='cleaned_'+temp[-1]
                        s='/'
                        temp2=s.join(temp)
                        self.printToLog('Forward read: '+temp2)
                        temp=self.readLocation[reverse].split('/')
                        temp[-1]='cleaned_'+temp[-1]
                        s='/'
                        temp2=s.join(temp)
                        self.printToLog('Reverse read: '+temp2)
                    else:
                        self.printToLog('Forward read: '+self.readLocation[forward])
                        self.printToLog('Reverse read: '+self.readLocation[reverse])
            elif len(self.readLocation)==1:  #if exactly 1 read is selected
                if (self.readLocation[0].endswith('.fq') | self.readLocation[0].endswith('.fastq')):  #check if read ends with .fq or .fastq
                    mode='single'
                    mode2='fq'
                elif (self.readLocation[0].endswith('.fasta')| self.readLocation[0].endswith('.fa')):
                    mode='single'
                    mode2='fa'
                else:
                    stop=True
                    self.printToLog("Read "+self.readLocation+' does not have correct file extenstion (expecting .fq or .fastq)') #print error if not ending in .fq or .fastq
            else:
                stop=True
                self.printToLog('No Read(s) Selected')
            if not stop:
                changeForFinal=True # change for final code (delete) ### change to True
                if self.doReadClean and mode=='single' and changeForFinal :
                    test=[]
                    for i in range (0,len(self.readLocation)):
                        temp=self.readLocation[i].split('/')
                        temp[-1]='cleaned_'+temp[-1]
                        s='/'
                        temp2=s.join(temp)
                        test.append(temp2)
                    self.readLocation=test
                self.printToLog('Alignment check complete!')
                self.printToLog('Beginning Alignment of Reads to Viral Reference: '+os.path.basename(self.ViralRefFa))
                if self.autoMode:
                    self.viralAlignment(mode,mode2,forward,reverse)
                else:
                    self.printToLog('Please confirm alignment...(Y/N) button')
                    self.app.yesbttn.configure(command=lambda: self.viralAlignment(mode,mode2,forward,reverse))
                    self.app.nobttn.configure(command=self.cancel_Alignment)
                
                
    def cancel_Alignment(self):
        self.printToLog('Alignment canceled!')
        self.app.resetYesNoButtons()
        self.GO=False
        self.event.set()
        
    def updateRunDirectory(self):
        if self.outDir=='None Selected':
            self.runDirectory=self.readDir+'/'+self.readBaseName
        else:
            self.runDirectory=self.outDir+'/'+self.readBaseName
    
    def viralAlignment(self,mode,mode2,forward,reverse):
        stop=False
        if mode=='single' or mode=='paired':
            self.mode=mode
        try:
            self.Trim=(self.app.TrimmingText.get(1.0,'end-1c'))
            self.Trim2=(self.app.TrimmingText2.get(1.0,'end-1c'))
            self.printToLog('Trimming '+self.Trim+' bases off 5\' end')
            self.printToLog('Trimming '+self.Trim2+' bases off 3\' end')
        except:
            self.printToLog('Trimming not entered correctly in config')
            err=traceback.format_exc()
            self.printToLog(err)
            stop=True
        if not stop:
            self.printToLog('Starting Viral Alignment...')
            os.chdir(self.ViralRefFolder)
            try:
                self.app.update()
                self.updateRunDirectory()
                if not os.path.exists(self.runDirectory):
                    os.makedirs(self.runDirectory)
                self.printToLog('Aligning...')
                os.chdir(self.ViralRefFolder)
                outfile=self.runDirectory+'/'+self.readBaseName+'_viralAlign.sam'
                if os.path.isfile(outfile) and self.loadAlignments:
                    self.printToLog('Loaded Prexisting Alignment:\n'+outfile)
                else:
                    if mode=='paired':
                        if self.platformType == "Windows":
                            temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2.bat','--trim5',str(self.Trim),'--trim3',str(self.Trim2),'--very-sensitive-local','--threads',str(self.coreOption),'-x',str(self.viralPrefix),'-1',self.readLocation[forward],'-2',self.readLocation[reverse],'-S',outfile],stderr=subprocess.PIPE)
                        else:
                            temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2','--trim5',str(self.Trim),'--trim3',str(self.Trim2),'--very-sensitive-local','--threads',str(self.coreOption),'-x',str(self.viralPrefix),'-1',self.readLocation[forward],'-2',self.readLocation[reverse],'-S',outfile],stderr=subprocess.PIPE)
                        #temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2.bat','--trim5',str(self.Trim),'--trim3',str(self.Trim2),'-D 20 -R 3 -N 0 -L 10 -i S,1,0.30','--threads',str(self.coreOption),'-x',str(self.viralPrefix),'-1',self.readLocation[forward],'-2',self.readLocation[reverse],'-S',outfile],stderr=subprocess.PIPE)
                    elif mode=='single':
                        #temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2.bat','--local','-D','20','-R','3','-N','1','-L','12','-i','S,1,0','--trim5',str(self.Trim),'--trim3',str(self.Trim2),'--threads',str(self.coreOption),'-x',str(self.viralPrefix),'-U',self.readLocation[0],'-S',outfile],stderr=subprocess.PIPE)
                        if self.platformType == "Windows":
                            temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2.bat','--trim5',str(self.Trim),'--trim3',str(self.Trim2),'--very-sensitive-local','--threads',str(self.coreOption),'-x',str(self.viralPrefix),'-U',self.readLocation[0],'-S',outfile],stderr=subprocess.PIPE)
                        else:
                            temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2','--trim5',str(self.Trim),'--trim3',str(self.Trim2),'--very-sensitive-local','--threads',str(self.coreOption),'-x',str(self.viralPrefix),'-U',self.readLocation[0],'-S',outfile],stderr=subprocess.PIPE)
                    for line in iter(temp.stderr.readline,b''):
                        self.printToLog(line.decode('utf-8'))
                    temp.stderr.close()
                self.printToLog('======================================================================\nAlignment Complete')
                self.GO=True
                self.event.set()
            except:
                err=traceback.format_exc()
                self.printToLog('Error Occured\n'+err)  
                self.GO=False
                self.event.set()
            os.chdir(self.workingDirectory)
            self.app.resetYesNoButtons()           

    def askHostAlignment(self):
        self.updateRunDirectory()
        self.CSLocation=self.runDirectory+'/'+self.readBaseName+'_CSList.txt'
        self.CSLocationFull=self.runDirectory+'/'+self.readBaseName+'_CSListFull.txt'
        if os.path.isfile(self.CSLocation) and os.path.isfile(self.CSLocationFull):
            self.printToLog('Alignment check complete!')
            self.printToLog('Beginning Alignment of Reads to Host Reference: '+os.path.basename(self.HostRefFa))
            if self.autoMode:
                self.HostAlignment()
            else:
                self.printToLog('Please confirm alignment...(Y/N) button')
                self.app.yesbttn.configure(command=self.HostAlignment)
                self.app.nobttn.configure(command=self.cancel_Alignment)
        else:
            self.printToLog('No Clipped Sequences Found to Align')
            self.GO=False
                     
        
    def HostAlignment(self):
        self.printToLog('Starting Host Alignment...')
        os.chdir(self.HostRefFolder)
        try:
            self.app.update()
            if not os.path.exists(self.runDirectory):
                os.makedirs(self.runDirectory)
            self.printToLog('Aligning...')
            os.chdir(self.HostRefFolder)
            outfile=self.runDirectory+'/'+self.readBaseName+'_hostAlign.sam'
            if os.path.isfile(outfile) and self.loadAlignments:
                self.printToLog('Loaded Prexisting Alignment:\n'+outfile)
            else:
                if self.overlapRequested:
                    if self.platformType == "Windows":
                        temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2.bat','--very-sensitive-local','--threads',str(self.coreOption),'-x',str(self.HostPrefix),'-f','-U',self.CSLocationFull,'-S',outfile],stderr=subprocess.PIPE)
                    else:
                        temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2','--very-sensitive-local','--threads',str(self.coreOption),'-x',str(self.HostPrefix),'-f','-U',self.CSLocationFull,'-S',outfile],stderr=subprocess.PIPE)
                else:
                    if self.platformType == "Windows":
                        temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2.bat','--very-sensitive-local','--threads',str(self.coreOption),'-x',str(self.HostPrefix),'-f','-U',self.CSLocation,'-S',outfile],stderr=subprocess.PIPE)
                    else:
                        temp=subprocess.Popen([self.Bowtie2Folder+'/bowtie2','--very-sensitive-local','--threads',str(self.coreOption),'-x',str(self.HostPrefix),'-f','-U',self.CSLocation,'-S',outfile],stderr=subprocess.PIPE)
                for line in iter(temp.stderr.readline,b''):
                    self.printToLog(line.decode('utf-8'))
                temp.stderr.close()
            self.printToLog('======================================================================\nAlignment Complete')
            self.GO=True
            self.event.set()
        except:
            err=traceback.format_exc()
            self.printToLog('Error Occured\n'+err)  
            self.GO=False
            self.event.set()
        os.chdir(self.workingDirectory)
        self.app.resetYesNoButtons()                


        
    
    def setSingleMode(self):
        if self.folderMode==True:
            self.printToLog('Changed Mode: Single Run')
            self.folderMode=False
    def setDirMode(self):
        if self.folderMode==False:
            self.printToLog('Changed Mode: Directory Run')
            self.folderMode=True
          
        
    def getCS(self):
        self.printToLog('Retrieving Clipped Sequences...')
        try:
            self.app.update()
            self.vSam=self.createSam('v')
            if self.vSam is not None:
                outfile=self.runDirectory+'/'+self.readBaseName+'_CSListFull.txt'
                if os.path.isfile(outfile) and self.loadAlignments:
                    self.printToLog('Loaded Prexisting Clipped Sequences:\n'+outfile)
                else:
                    self.extractClippedSequences()
                self.printToLog('======================================================================\nClipped Sequence Files Generated') 
                self.GO=True
            else:
                self.GO=False
        except:
            err=traceback.format_exc()
            self.printToLog('Error Occured\n'+err)
            self.GO=False
        
           


    def createSam(self,mode):
        switch=True
        if len(self.readLocation)==0:
            switch=False
        if switch:
            self.updateRunDirectory()
            if mode=='v':
                filename=self.runDirectory+'/'+self.readBaseName+'_viralAlign.sam'
                self.printToLog('createSam V, filename = '+ filename)
                if os.path.isfile(filename):
                    vSam=Sam(filename,False,self.readBaseName)
                    if vSam.Quit==False:
                        return vSam
                    else:
                        self.printToLog('Zero Alignment for the Virus:\n')
                        return None
                else:
                    self.printToLog('No Viral Alignment File exists under name:\n'+filename)
                    return None
                
            elif mode=='h':
                filename=self.runDirectory+'/'+self.readBaseName+'_hostAlign.sam'
                self.printToLog('createSam H, filename = '+ filename)
                if os.path.isfile(filename):
                    hSam=Sam(filename,False,self.readBaseName)
                    return hSam
                else:
                    self.printToLog('No Host Alignment File exists under name:\n'+filename)
                    return None
            else:
                self.printToLog('Error in reading sam file: no type switch given (virus/host)')
                return None
        else:
            self.printToLog('No Reads Selected')
            return None

    def sequence_cleaner(self,fastq_file, min_length=0, per_n=100):
        self.printToLog('Cleaning: '+fastq_file)
        sequences=[]
        test=[]
        dup=[]
        count1 = 0
        count2 = 0
        for seq_record in SeqIO.parse(fastq_file, "fastq"):
            count1+=1
            sequence = str(seq_record.seq).upper()
            if (len(sequence) >= min_length and (float(sequence.count("N"))/float(len(sequence)))*100 <= per_n):
                sequences.append(sequence)
                count2+=1
                test.append(seq_record)
                if (sequence.count("N") > 0):
                    self.printToLog("cleaning N="+str(sequence.count("N"))+" at "+str(count1)+" ")
                ### FSS -- the unique check is O(n^2) time. moved to postprocessing for less reads to be checked
                """
                if sequence not in sequences:
                    sequences.append(sequence)
                    count3 = count3 + 1 ### FSS                
                    test.append(seq_record)
                """
            #else: self.printToLog("FSS cleaning remove length = "+str(len(sequence)))
        self.printToLog("cleaning total="+str(count1)+" good="+str(count2)+" ")
        if fastq_file.endswith('1.fq') or fastq_file.endswith('1.fastq'):
            temp='1.fastq'
        elif fastq_file.endswith('2.fq') or fastq_file.endswith('2.fastq'):
            temp='2.fastq'
        else:
            temp='unknown.fastq'
            self.printToLog('Error in fastQ cleanup.  One or more of the file names did not end with 1.fastq/1.fq or 2.fastq/2.fq')
        SeqIO.write(test,self.readDir+'/cleaned_'+self.readBaseName+temp,'fastq')
        temp=fastq_file.split('/')
        temp[-1]='cleaned_'+temp[-1]
        s='/'
        temp2=s.join(temp)
        self.printToLog('Wrote cleaned fastq to: '+temp2)
        self.event.set()


    def extractClippedSequences(self):
        stop=False
        try:
            self.vSam
        except NameError:
            self.printToLog('No alignment file to virus was found')
            stop=True
        try:
            self.CSminLen=int(self.app.clipText.get(1.0,'end-1c'))
        except:
            self.printToLog('Clipped Sequence Min Length not entered correctly in config')
            err=traceback.format_exc()
            self.printToLog(err)
            stop=True
        if not stop:
            self.vSam.extractClippedSeq(self.CSminLen,self.paired)
            
    
    def mapper(self):
        vSam=self.createSam('v')
        hSam=self.createSam('h')

        try:
            if vSam is not None and hSam is not None:
                carryOn=True
            else:
                carryOn=False
        except:
            self.printToLog('One or more .sam files are not loaded')
            carryOn=False
        if carryOn:
            self.Map=Mapper(vSam,hSam,self)
        self.event.set()

    def countUniqueForHighCount(self, countUniqueOnly=False):
        list_data=[]
        for i in range(len(self.Map.data)):
            sequence=self.Map.data[i]['Sequence']
            sequence_length=len(sequence)
            sequence_complement=str(Seq(sequence).reverse_complement())
            sequence_hostLength=int(self.Map.data[i]['Hlength'])
            if (self.Map.data[i]['Similarity'] == True):
                similarity_found=True
            else:
                similarity_found=False
            if (self.Map.data[i]['Microhomology'] == True) or (self.Map.data[i]['Filtered'] == True):
                filtered=True
            else:
                filtered=False
            #................[0].[1]..............[2].......[3]..................[4]..................[5]...[6]...............[7]......
            list_data.append([i, sequence_length, sequence, sequence_complement, sequence_hostLength, True, similarity_found, filtered])

        new_list_data=sorted(list_data, key=lambda item:item[2]) # sorting on sequence 
        new_list_data2=sorted(new_list_data, key=lambda item:item[1]) # sorting on sequence length

        # check unique
        self.printToLog('Checking unique reads...')
        start_time = time.time()
        last_unique_index=0
        unique_count = 1
        dup_count=0
        for i in range(1, len(new_list_data2)):
            if (new_list_data2[i][1] == new_list_data2[last_unique_index][1]): #same length
                if (new_list_data2[i][2] != new_list_data2[last_unique_index][2]): # compare sequence
                    last_unique_index=i
                    unique_count+=1
                else:
                    dup_count +=1
                    new_list_data2[i][5]=False
            else:
                last_unique_index=i
                unique_count+=1
        end_time = time.time()
        self.printToLog(str(unique_count)+' unique reads')

        if (countUniqueOnly == True):
            return
        
        # check complement
        start_time = time.time()
        self.printToLog('Checking complement reads...')
        complement_count=0
        for i in range(0, len(new_list_data2)-1):
            if new_list_data2[i][5]==False:
                continue
            for j in range(i+1, len(new_list_data2)):
                if new_list_data2[j][5]==False:
                    continue
                diff=new_list_data2[j][1]- new_list_data2[i][1]
                if (diff > 0):
                    break
                if (new_list_data2[i][2]==new_list_data2[j][3]):
                    complement_count +=1
                    new_list_data2[j][5]=False
        end_time = time.time()
        self.printToLog(str(unique_count-complement_count)+' unique reads (complement removed)')

        # update original similarity data
        for i in range(0, len(new_list_data2)):
            if (new_list_data2[i][5] == False):
                ndx=new_list_data2[i][0]
                self.Map.data[ndx]['Similarity']=True
            if (new_list_data2[i][6] == True) or (new_list_data2[i][6] == True):
                new_list_data2[i][5] = False
            
        self.readSimilarity=int(self.app.readSimilarityText.get(1.0,'end-1c'))
        if (self.readSimilarity >= 100):
            return
        start_time = time.time()
        self.printToLog("Checking similarity...")
        level=self.readSimilarity
        level_p=(100.0-level)/100.0
        similarity_count=0
        similarity_count_C=0
        #for i in range(len(new_list_data2), len(new_list_data2)-1):
        for i in range(0, len(new_list_data2)-1):
            if (i% 1000)==0:
                self.printToLog("similarity: processing "+str(i)+" of "+str(len(new_list_data2)))
            if new_list_data2[i][5]==False:
                continue
            for j in range(i+1, len(new_list_data2)):
                if new_list_data2[j][5]==False:
                    continue
                up_diff=int(new_list_data2[j][1]*level_p)
                diff=new_list_data2[j][1]- new_list_data2[i][1]
                if (up_diff < 6):
                    up_diff=6
                if (diff) > up_diff: #stop similarity check 
                    break
                if (abs(new_list_data2[j][4]- new_list_data2[i][4])>(up_diff)):
                    break
                similarity=(difflib.SequenceMatcher(None, new_list_data2[i][2], new_list_data2[j][2]).ratio())*100.0
                if similarity >= level:
                    if new_list_data2[j][1]> new_list_data2[i][1]: # keep longer one
                        new_list_data2[i][5]=False
                        similarity_count+=1
                        #update
                        ndx=new_list_data2[i][0]
                        self.Map.data[ndx]['Similarity']=True
                        break
                    else:
                        new_list_data2[j][5]=False
                        #update
                        ndx=new_list_data2[j][0]
                        self.Map.data[ndx]['Similarity']=True
                        similarity_count+=1
                    continue        
                similarity=(difflib.SequenceMatcher(None, new_list_data2[i][2], new_list_data2[j][3]).ratio())*100.0
                if similarity >= level:
                    if new_list_data2[j][1]> new_list_data2[i][1]: # keep longer one
                        new_list_data2[i][5]=False
                        similarity_count+=1
                        similarity_count_C+=1
                        #update
                        ndx=new_list_data2[i][0]
                        self.Map.data[ndx]['Similarity']=True
                        break
                    else:
                        new_list_data2[j][5]=False
                        similarity_count+=1
                        similarity_count_C+=1
                        #update
                        ndx=new_list_data2[j][0]
                        self.Map.data[ndx]['Similarity']=True
                    continue        
        end_time = time.time()
        self.printToLog("similarity check complete")

        return


    
    def countUnique(self):
        if (len(self.Map.data))>1200:
            self.countUniqueForHighCount()
            return
        else:
            self.countUniqueForHighCount(True)
        list1=[]
        list1_length=[0]*len(self.Map.data) # sequence host length
        list1_index=[0]*len(self.Map.data) # sequence index in the map.data
        count=0
        self.readSimilarity=int(self.app.readSimilarityText.get(1.0,'end-1c'))
        #self.printToLog('Similarity set to: '+str(self.readSimilarity))
        """
        self.printToLog('write reads to temp file =============================== ')
        out = out=open('c:/trash/test_reads.txt','w')
        out.write(str(len(self.Map.data))+"\n")
        for i in range(0,len(self.Map.data)):
            seq=self.Map.data[i]['Sequence']
            seq_length=int(self.Map.data[i]['Hlength'])
            out.write(str(seq_length)+"\n")
            out.write(seq+"\n")
        out.close()
        self.printToLog('write reads to temp file complete ====================== ')
        """
        try:
            ndx=-1
            for i in range(0,len(self.Map.data)):
                ndx +=1
                seq=Seq(self.Map.data[i]['Sequence'])
                seq_length=int(self.Map.data[i]['Hlength'])
                temp=str(seq)
                temp2=str(seq.reverse_complement())

                if (self.Map.data[i]['Microhomology'] == True) or (self.Map.data[i]['Filtered'] == True) \
                   or (self.Map.data[i]['Similarity'] == True):
                    continue
                    
                if list1.count(temp)==0:
                    ### FSS check similarity
                    skip=False
                    a = len(temp)
                    index = 0
                    for data in list1:
                        str_length = list1_length[index]
                        str_index = list1_index[index]
                        if self.checkSimilarity(temp, data, int(str_length), int(seq_length))==True:
                            if int(seq_length) > int(str_length):
                                #self.printToLog('similarity list1 size = '+str(len(list1))+' list_index '+str(index)+' length = in list '+str(str_length)+" "+str(str_index)+" :  new seq "+str(seq_length)+" "+str(ndx))
                                list1.remove(data)
                                self.Map.data[str_index]['Similarity'] = True #set old item's similarity to be True
                                list1_index=list1_index[:index]+list1_index[index+1:]
                                list1_length=list1_length[:index]+list1_length[index+1:]
                                skip = False
                            else:    
                                skip = True
                            break                            
                        index +=1
                    if not skip:            
                        list1.append(temp)
                        list1_length[len(list1)-1]=seq_length
                        list1_index[len(list1)-1]=ndx
                        if list1.count(temp2)==0:
                            count+=1
                    else:
                        self.Map.data[i]['Similarity'] = True
                else:
                    self.Map.data[i]['Similarity'] = True

            """
            if self.mode=='paired':
                self.printToLog(str(len(list1))+' reads are unique (counting complements)')
            self.printToLog(str(count)+' reads are unique')
            """
                
        except:
            self.printToLog('Something went wrong')
            err=traceback.format_exc()
            self.printToLog(err)

    def checkSimilarity(self, seq1, seq2, host_length1=0, host_length2=0):
        similarityFound=False
        a = len(seq1)
        b = len(seq2)
        if (abs(a-b))<=6: # only check for length difference <= 6 to save check time
            #if (seq1[0]==seq2[0]) or (seq1[a-1]==seq2[b-1]):
            similarity=(difflib.SequenceMatcher(None, seq1, seq2).ratio())*100.0
            if (similarity >= self.readSimilarity):
                similarityFound = True
            if not similarityFound:
                seq1_r=str(Seq(seq1).reverse_complement())
                #if (seq1_r[0]==seq2[0]) or (seq1_r[a-1]==seq2[b-1]):
                similarity=(difflib.SequenceMatcher(None, seq1_r, seq2).ratio())*100.0
                if (similarity >= self.readSimilarity):
                    similarityFound = True
        return similarityFound

    def finalCheck(self):
        self.printToLog('final similarity check!!!!!!!!!')
        count=0
        count1=0
        for i in range(0, len(self.savedReads)):
            comp=str(Seq(self.savedReads[i]).reverse_complement())
            for j in range(i+1, len(self.savedReads)):
                if (self.checkSimilarity(self.savedReads[i], self.savedReads[j])==True):
                    count+=1
                    self.printToLog('final similarity found '+str(i)+"  "+str(j)+"  count = "+str(count))
                    
    def saveMultipleFinalReads(self):
        # the column titles have been created already
        with open (self.savedReadsFilenameFinal, 'a') as csvfile1:
            writer=csv.writer(csvfile1,delimiter=',',lineterminator = '\n')
            for row in self.multipleFileFinalReads:
                writer.writerow(row)
        self.printToLog('\n\n**********************************************************************')
        self.printToLog('Multiple files final result '+str(len(self.multipleFileFinalReads))+' reads saved to file:')
        self.printToLog(self.savedReadsFilenameFinal)
        self.printToLog('**********************************************************************\n')
        self.Map.data = []
        ndx= 0
        for row in self.savedReadsMapData:
            row['Index']=ndx # reset index for new list to be displayed
            self.Map.data.append(row)
            ndx +=1
        self.processLargeFile = False
        self.app.listRedraw()
        return
    
    def saveData(self):
            
        def target():
            """
            if (len(self.Map.data))==0:
                self.printToLog('No new data!')
                return
            """
            try:
                print(self.Map)
            except AttributeError:
                self.printToLog('Data has not been generated for saving')
                return
            self.printToLog('Saving data...')
            ts=time.time()
            st=datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%Hh%Mm%Ss')
            try:
                out=open(self.runDirectory+'/'+self.readBaseName+'_OUTPUT_'+st+'.csv','w')
                go=True
            except PermissionError:
                self.printToLog('Cannot save data, file currently in use')
                go=False
            if go:
                start_time = time.time()
                sequences=[]
                dataNoC=[]
                data=[]
                dataNoC_L=[]
                
                attributeList= ['ReadName','ReadLength','Sequence','Chromosome','Gene','InsideGene','DistanceToGene','GeneDirection','Focus','GeneObj','FocusObj','HostLocalCords','Hlength','HostRefCords','HostOrientation','ViralAccession','ViralLocalCords','Vlength','ViralRefCords','ViralOrientation','HTM','HTMAdjusted','VTM','VTMAdjusted','Overlap','OverlapTM','OverlapTMAdjusted','Inserted','Microhomology','HostMapFlag','ViralMapFlag','HMapQ','VMapQ','FastqSource','Index']
                for a in attributeList:
                    out.write(a+',')
                out.write('\n')
                #tempList=[]
                cnt=0
                for i in range (0,len(self.Map.data)):
                    if (self.Map.data[i]['Microhomology']==False) and (self.Map.data[i]['Filtered']==False) and (self.Map.data[i]['Similarity']==False):
                        for item in attributeList:
                                out.write('"'+str(self.Map.data[i][item])+'"'+',')
                        out.write('\n')
                        cnt+=1
                out.close()
                records_saved = cnt
                """
                # save all reads for debug purpose
                self.printToLog('save all reads for debug purpose............')
                if not(self.processLargeFile) or (self.processLargeFile):
                    out=open(self.runDirectory+'/'+self.readBaseName+'_DEBUG_'+st+'.csv','w')
                    attributeList_debug= ['ReadName','ReadLength','Sequence','Chromosome','Gene','InsideGene','DistanceToGene','GeneDirection','Focus','GeneObj',\
                                          'FocusObj','HostLocalCords','Hlength','HostRefCords','HostOrientation','ViralAccession','ViralLocalCords','Vlength',\
                                          'ViralRefCords','ViralOrientation','HTM','HTMAdjusted','VTM','VTMAdjusted','Overlap','OverlapTM','OverlapTMAdjusted',\
                                          'Inserted','Microhomology','Filtered','Similarity','HostMapFlag','ViralMapFlag','HMapQ','VMapQ','FastqSource','Index']
                    for a in attributeList_debug:
                        out.write(a+',')
                    out.write('\n')
                    #tempList=[]
                    for i in range (0,len(self.Map.data)):
                        for item in attributeList_debug:
                                out.write('"'+str(self.Map.data[i][item])+'"'+',')
                        out.write('\n')
                    out.close()
                """    
                firstSave=False
                if (self.processLargeFile): # handle large file
                    if len(self.savedReads)==0:
                        firstSave=True
                        self.savedReadsDir = self.processLargeFileDir
                        self.savedReadsFilename= self.savedReadsDir+"/"+self.readBaseName+'_OUTPUT_'+st+'.csv'
                        self.savedReadsFilenameFinal= self.savedReadsDir+"/"+self.readBaseName+'_final_OUTPUT_'+st+'.csv'
                        self.savedReadsLogFilename= self.savedReadsDir+"/"+self.readBaseName+'_log_'+st+'.txt'
                        out_L=open(self.savedReadsFilename,'w')
                        out_L.close()
                        out_LD=open(self.savedReadsFilenameFinal,'w')
                        # create column title for final result file
                        for a in attributeList:
                            out_LD.write(a+',')
                        out_LD.write('\n')
                        out_LD.close()
                        
                with open (self.runDirectory+'/'+self.readBaseName+'_OUTPUT_'+st+'.csv','r') as csvfile:
                    reader = csv.reader(csvfile, delimiter=',')
                    read_index = 0
                    for row in reader:
                        seq=Seq(row[2])
                        try:
                            hostLength=int(row[12])
                        except:
                            hostLength=0
                        if seq not in sequences:
                            #check for similarity
                            skip=False
                            if (records_saved <= 200): #only check similarity again for small number of records, otherwise time to save will be long
                                for sequence in sequences:
                                    if (self.checkSimilarity(str(seq), str(sequence)) == True):
                                        skip=True
                                        break
                            #else:
                            #    self.printToLog('similarity ccheck skip for high records count')                
                            if not skip:
                                data.append(row)
                                if (records_saved > 200) or (seq.reverse_complement() not in sequences):
                                    dataNoC.append(row)
                                    if (self.processLargeFile): # handle large file
                                        temp1 = str(seq)
                                        temp1_r = str(seq.reverse_complement())
                                        if (read_index > 0):
                                            actual_index = int(row[34])
                                        if firstSave:
                                            if (read_index == 0): # keep first row for attribute
                                                dataNoC_L.append(row)
                                            if (read_index > 0): #skip first row (attributes)
                                                if (self.Map.data[actual_index]['Similarity']!=True) \
                                                   and (self.Map.data[actual_index]['Filtered']!=True) \
                                                   and (self.Map.data[actual_index]['Microhomology']!=True):
                                                    dataNoC_L.append(row)
                                                    self.multipleFileFinalReads.append(row)
                                                    self.savedReads.append(str(seq))
                                                    self.savedReadsHlength.append(hostLength)
                                                    self.savedReadsMapData.append(self.Map.data[actual_index])
                                                else: self.printToLog('FSS remove 1st save due to similarity filter homology')
                                        else:
                                            if (read_index > 0):
                                                actual_index = int(row[34])
                                                if (self.Map.data[actual_index]['Similarity']!=True) \
                                                   and (self.Map.data[actual_index]['Filtered']!=True) \
                                                   and (self.Map.data[actual_index]['Microhomology']!=True):
                                                    if (temp1 not in self.savedReads) and (temp1_r not in self.savedReads):
                                                        a = len(temp1)
                                                        skip_L=False
                                                        index=0
                                                        for r in self.savedReads:
                                                            hostLength2=self.savedReadsHlength[index]
                                                            if (self.checkSimilarity(temp1, r)==True):
                                                                if (hostLength > hostLength2):
                                                                    loc=self.savedReads.index(r)
                                                                    self.savedReads=self.savedReads[:loc]+self.savedReads[loc+1:]
                                                                    self.savedReadsHlength=self.savedReadsHlength[:loc]+self.savedReadsHlength[loc+1:]
                                                                    self.savedReadsMapData=self.savedReadsMapData[:loc]+self.savedReadsMapData[loc+1:]
                                                                    for k in range(len(self.multipleFileFinalReads)):
                                                                        read1=self.multipleFileFinalReads[k][2]
                                                                        if (r == read1):
                                                                            self.multipleFileFinalReads = self.multipleFileFinalReads[:k] + self.multipleFileFinalReads[k+1:]
                                                                            break;
                                                                    skip_L = False
                                                                else:
                                                                    skip_L = True                                                            
                                                                break
                                                            index+=1
                                                        if not skip_L:
                                                            dataNoC_L.append(row)
                                                            self.multipleFileFinalReads.append(row)
                                                            self.savedReads.append(temp1)
                                                            self.savedReadsHlength.append(hostLength)
                                                            self.savedReadsMapData.append(self.Map.data[actual_index])
                                                        #else: self.printToLog('FSS remove 1st save due to similarity on previous file')                                                        
                                                    #else: self.printToLog('FSS remove 1st save due to existing on previous file')
                                                else: self.printToLog('FSS remove regular save due to similarity filter homology')
                                #else: self.printToLog('FSS remove due to complement')
                                sequences.append(seq)
                            #else: self.printToLog('FSS remove due to similarity')
                        #else: self.printToLog('FSS remove due to existing duplicate')
                        read_index +=1
                if (self.processLargeFile==False):
                    filename1 = self.runDirectory+'/'+self.readBaseName+'_OUTPUT_'+st+'_no_duplicates.csv'
                    count1 = 0
                    with open (filename1, 'w') as csvfile1:
                        writer=csv.writer(csvfile1,delimiter=',',lineterminator = '\n')
                        for row in data:
                            writer.writerow(row)
                            count1 +=1
                    #self.printToLog(str(count1)+' reads saved to '+filename1+'\n')

                    if len(self.savedReadsMapData)==0:
                        filename2 = self.runDirectory+'/'+self.readBaseName+'_OUTPUT_'+st+'_no_complements.csv'        
                        count2 = 0
                        with open (filename2, 'w') as csvfile2:
                            writer=csv.writer(csvfile2,delimiter=',',lineterminator = '\n')
                            for row in dataNoC:
                                writer.writerow(row)
                                count2 +=1
                        if (count2 > 0):
                            self.printToLog(str(count2-1)+' reads saved to '+filename2+'\n')
                        else:
                            self.printToLog('0 reads saved to '+filename2+'\n')
                    else: # save reads after all multiple files procesed completely
                        filename2 = self.savedReadsDir+"/"+self.readBaseName+'_revised_OUTPUT_'+st+'.csv'        
                        count2 = 0
                        with open (filename2, 'w') as csvfile2:
                            writer=csv.writer(csvfile2,delimiter=',',lineterminator = '\n')
                            for row in dataNoC:
                                writer.writerow(row)
                                count2 +=1
                        if (count2 > 0):
                            self.printToLog(str(count2-1)+' reads saved to '+filename2+'\n')
                        else:
                            self.printToLog('0 reads saved to '+filename2+'\n')
                        
                if (self.processLargeFile):
                    count3 = 0
                    filename3 = self.savedReadsFilename
                    with open (filename3, 'a') as csvfile3:
                        writer=csv.writer(csvfile3,delimiter=',',lineterminator = '\n')
                        for row in dataNoC_L:
                            writer.writerow(row)
                            count3 +=1
                    if (count3 > 0):        
                        self.printToLog(str(count3-1)+' reads appended to '+filename3+'\n')
                    else:
                        self.printToLog('0 reads appended to '+filename3+'\n')
                    accumulated_reads=len(self.savedReads)
                    self.printToLog('accumulated reads = '+str(accumulated_reads)+'\n')
                    self.changeProcessingSaveDataStatus(False)

                    #self.finalCheck()

                end_time = time.time()        
                self.printToLog('Data saved successfully!')
            try: 
                if (self.processLargeFile == False):
                    out2=open(self.runDirectory+'/'+self.readBaseName+'_log.txt','w')
                    go=True
                else:
                    if (self.savedReadsLogFilename != ""):
                        out2=open(self.savedReadsLogFilename, "w")
                        self.printToLog('accumulated time = '+str(end_time - self.splitStartTime)+' seconds')
                        go=True
                    else:
                        go=False
            except:
                self.printToLog('Cannot save log, file currently in use')
                go=False
            if go:
                for line in self.logger:
                    out2.write(line+'\n')
                out2.close()
                self.printToLog('Log saved successfully!')
            self.event.set()
                
        thread=Thread(target=target)
        thread.start()
        self.event.wait()
        self.event.clear()
        self.printToLog('======================================================================\n')
        

    

class Mapper:
    def __init__(self,vSam,hSam,core):
        self.readCount=0
        self.data=[]
        self.vSam=vSam
        self.hSam=hSam
        self.core=core
        self.gatherAligned()
        
    
    def calculateTM(self,seq):
        try:
            self.core.Na=int(self.core.app.saltText.get(1.0,'end-1c'))           
        except:
            return 'Err','Err'
        Na=self.core.Na
        a=seq.count('A')
        c=seq.count('C')
        t=seq.count('T')
        g=seq.count('G')
        if len(seq)<14:
            TMb=(a+t)*2+(g+c)*4
            TMa=(a+t)*2+(g+c)*4-16.6*math.log(0.05,10)+16.6*math.log((Na*10**-3),10)
        else:
            TMb=64.9+41*(g+c-16.4)/(a+t+g+c)
            TMa=100.5+(41*(g+c)/(a+t+g+c))-(820/(a+t+g+c))+16.6*math.log((Na*10**-3),10)
        return TMb,TMa

            
    
    def gatherAligned(self):
        code=1
        stop=False
        if (self.core.app.mh_stringency.get(1.0,'end-1c').strip() == ""): # microhomology stringency % is empty
            microhomology_percent = 95
            microhomology_nt = int(self.core.app.mh_stringency1.get(1.0,'end-1c'))
        else:
            microhomology_percent = int(self.core.app.mh_stringency.get(1.0,'end-1c'))
            microhomology_nt = 0
        self.core.app.listBox1.delete(0,END)
        increment=len(self.hSam.sam)/20.0
        progress=0
        for i in range (0,len(self.hSam.sam)):
            tempProgress=int(i/increment)
            if tempProgress>progress:
                self.printToLog('Progress: '+str(tempProgress*5)+'%')
                progress=tempProgress
            stop=False
            if self.hSam.sam[i][5]=='*':
                stop=True
            if not stop:
                if code==0:
                    name=self.hSam.sam[i][0].split('.')
                    if self.core.paired:
                        name=(name[1]+'.'+name[2]).split('-')
                    else:
                        name=name[1].split('-')
                    name=name[1].split(':')
                    name=':'.join(name[2:len(name)])
                if code==1:
                    name=self.hSam.sam[i][0]
                self.core.app.listBox1.insert(END,name)
                self.readCount+=1
                self.loadData(i, microhomology_percent, microhomology_nt)    
        self.printToLog('Progress: 100%')                          
        self.printToLog(str(self.readCount)+' reads contain host and viral sequences')
    
    
    def loadData(self,number, microhomology_percent, microhomology_nt):
        def getOverlap(hstart,hstop,vstart,vstop):
            HostFirst=False
            if vstart<=hstart:
                if vstop>=hstart and hstop>=vstop:
                    overlap=vstop-hstart+1
                elif vstop>=hstart and hstop<vstop:
                    overlap=hstop-hstart+1
                else:
                    overlap=0
            elif hstart<vstart:
                HostFirst=True
                if hstop>=vstart and vstop>=hstop:
                    overlap=hstop-vstart+1
                elif hstop>=vstart and vstop<hstop:
                    overlap=vstop-vstart+1
                else:
                    overlap=0
            return [overlap,HostFirst]
        def getAlignmentLocal(cigar):
            if cigar!="*":
                pairs=re.findall('\d+[A-Z]+',cigar)
                index=0
                inserted=[]
                count=0
                deleted=0
                stop=None
                start=None
                stopSignal=True
                for j in range (0,len(pairs)):
                    length=int(pairs[j][0:len(pairs[j])-1])
                    if pairs[j].endswith('S'):
                        index+=length
                    elif pairs[j].endswith('M'):
                        stopSignal=True
                        if j<len(pairs):
                            for k in range(j+1,len(pairs)):
                                if pairs[k].endswith('M'):
                                    stopSignal=False
                        if count==0:
                            start=index+1
                        index+=length
                        if stopSignal:
                            stop=index
                        count+=1
                    elif pairs[j].endswith('D'):
                        deleted+=length
                    elif pairs[j].endswith('I'):
                        insertedStart=index
                        index+=length
                        insertedStop=index
                        inserted.append([insertedStart,insertedStop])
                        
                return [start,stop,inserted]
            else:
                print('!!! Error Detected - Cigar String Should be not empty (not *)')
                return [None,None]
        data={}
        code=1#for final code - can change depending on what naming specification is.  default =1
        if code==0:
            name=self.hSam.sam[number][0].split('.')
            if self.core.paired:
                name=(name[1]+'.'+name[2]).split('-')
            else:
                name=name[1].split('-')
            name=name[1].split(':')
            name=':'.join(name[2:len(name)])
        else:
            name=self.hSam.sam[number][0]
        #self.printToLog('FSS data = '+ str(self.hSam.sam[number]))
        data['ReadName']=name
        data['Index']=len(self.data)
        data['Chromosome']=self.hSam.sam[number][2][3:]
        data['HMapQ']=self.hSam.sam[number][4]
        n2=self.hSam.sam[number][0].split('.')
        if len(n2)>2:
            n=int(n2[-1])
            data['FastqSource']=self.core.readLocation[n-1]
        else:
            data['FastqSource']=self.core.readLocation[0]
        n=int(n2[0])
        data['ReadLength']=len(self.vSam.sam[n][9])
        data['VMapQ']=self.vSam.sam[n][4]
        data['HostMapFlag']=self.hSam.sam[number][1]
        data['ViralMapFlag']=self.vSam.sam[n][1]
        
        vswitch=False
        hswitch=False
        if len(self.vSam.flags[n])>=4:
            if self.vSam.flags[n][4]=='1':
                vswitch=True
        if len(self.hSam.flags[number])>=4:
            if self.hSam.flags[number][4]=='1':
                hswitch=True
                
        data['Vswitch']=vswitch
        data['Hswitch']=hswitch
        if not (vswitch):
            data['Sequence']=self.vSam.sam[n][9]
        else:
            fullSeq=Seq(self.vSam.sam[n][9])
            data['Sequence']=str(fullSeq.reverse_complement())
    
        if data['Vswitch']:
            data['ViralOrientation']='Plus/Minus'
        else:
            data['ViralOrientation']='Plus/Plus'
        if data['Hswitch']:
            data['HostOrientation']='Plus/Minus'
        else:
            data['HostOrientation']='Plus/Plus'
            
        cigar=self.vSam.sam[n][5]
        [vstart,vstop,inserted]=getAlignmentLocal(cigar)
        vlen=abs(vstop-vstart)+1
        cigar=self.hSam.sam[number][5]
        [hstart,hstop,inserted]=getAlignmentLocal(cigar)
        hlen=abs(hstop-hstart)+1
        
        data['ViralLocalCords']=[vstart,vstop]
        data['HostLocalCords']=[hstart,hstop]     
        data['Hlength']=hlen
        data['Vlength']=vlen
        
        
        
        if vswitch:
            temp=len(self.hSam.sam[number][9])+1-vstop
            vstop=len(self.hSam.sam[number][9])+1-vstart
            vstart=temp
            data['ViralLocalCords']=[vstart,vstop]
            if hswitch:
                data['HostLocalCords']=[hstart,hstop]
            else:
                temp=len(self.hSam.sam[number][9])+1-hstop
                hstop=len(self.hSam.sam[number][9])+1-hstart
                hstart=temp
                data['HostLocalCords']=[hstart,hstop]
        else:
            data['ViralLocalCords']=[vstart,vstop]
            if hswitch:
                temp=len(self.hSam.sam[number][9])+1-hstop
                hstop=len(self.hSam.sam[number][9])+1-hstart
                hstart=temp
                data['HostLocalCords']=[hstart,hstop]
            else:
                data['HostLocalCords']=[hstart,hstop]
        """
        if hstart>=vstart and hstop<=vstop:
            data['Microhomology']=True
        elif vstart>=hstart and vstop<=hstop:
            data['Microhomology']=True
        else:
            data['Microhomology']=False
        """
        
        
        seq=data['Sequence'][int(hstart)-1:int(hstop)]
        data['HostSeq']=seq
        basicTM,adjustedTM=self.calculateTM(seq)
        data['HTM']=basicTM
        data['HTMAdjusted']=adjustedTM
        data['Inserted']=inserted
        data['ViralAccession']=self.vSam.sam[n][2]
        seq=data['Sequence'][int(vstart)-1:int(vstop)]
        data['ViralSeq']=seq
        basicTM,adjustedTM=self.calculateTM(seq)
        data['VTM']=basicTM
        data['VTMAdjusted']=adjustedTM
        [overlap,HostFirst]=getOverlap(hstart,hstop,vstart,vstop)
        data['Overlap']=overlap
        data['HostFirst']=HostFirst
        if overlap>0: 
            if HostFirst:
                overlapSeq=data['Sequence'][int(hstop)-int(overlap):int(hstop)]
            else:
                overlapSeq=data['Sequence'][int(vstop)-int(overlap):int(vstop)]
            [basicTm,saltadjustedTm]=self.calculateTM(overlapSeq)
            data['OverlapTM']=basicTm
            data['OverlapTMAdjusted']=saltadjustedTm
        else:
            data['OverlapTM']=''
            data['OverlapTMAdjusted']=''
        data['HostRefCords']=[int(self.hSam.sam[number][3]),int(self.hSam.sam[number][3])+hlen]
        data['ViralRefCords']=[int(self.vSam.sam[n][3]),int(self.vSam.sam[n][3])+vlen]

        if (hlen > vlen):
            percent = int(float(overlap)/float(vlen)*100)
        else:
            percent = int(float(overlap)/float(hlen)*100)

        if (percent >= microhomology_percent) and (overlap >= microhomology_nt):
            #self.printToLog('Microhomology: '+str(percent)+"% "+str(overlap)+" <<<"+str(microhomology_percent)+" "+str(microhomology_nt)+">>>")                          
            data['Microhomology']=True
        else:
            data['Microhomology']=False

        # FSS added
        data['Similarity']=0.0
        data['Filtered']=False
        
        if self.core.doGTF:
            hostStart=data['HostRefCords'][0]
            hostStop=data['HostRefCords'][1]
            hostChrom=data['Chromosome']
            #print (hostStart,hostStop,hostChrom)
            [information,geneObj,focusObj]=self.searchGenes(hostChrom,hostStart,hostStop)
            data['Gene']=information['Gene']
            data['InsideGene']=information['Inside']
            data['GeneDirection']=information['Direction']
            data['DistanceToGene']=information['DistanceToGene']
            data['Focus']=information['Focus']
            data['GeneObj']=geneObj
            data['FocusObj']=focusObj
        else:
            data['Gene']='N/A'
            data['InsideGene']='N/A'
            data['GeneDirection']='N/A'
            data['DistanceToGene']='N/A'
            data['Focus']='N/A'
            data['GeneObj']='N/A'
            data['FocusObj']='N/A'
        self.data.append(data)
        
        
        
    def searchGenes(self,chrom,start,stop):
        distanceThreshold=self.core.distanceThreshold
        information={}
        information['Inside']='N/A'
        information['DistanceToGene']='N/A'
        information['Direction']='N/A'
        information['Gene']='N/A'
        information['Focus']='N/A'
        final=None
        finalData=None
        chrom=str(chrom)
        go=True
        focusGene=None
        focus=None
        try:
            index=self.core.ChromList.index(chrom)
        except:
            go=False
        if go:
            temp=[float(x[3]) for x in self.core.FocusData[index]]
            index2=bisect(temp,start)
            temp2=[float(x[3]) for x in self.core.GeneData[index]]
            index_gene=bisect(temp2,start)
            try:
                preData=self.core.FocusData[index][index2-1]
                pre=self.core.GeneData[index][index_gene-1]
            except:
                pre=None
                preData=None
            try:
                postData=self.core.FocusData[index][index2]
                post=self.core.GeneData[index][index_gene]
            except:
                post=None
                postData=None
            if not pre:
                final='post'
            elif not post:
                final='pre'
            if not preData:
                finalData='post'
            elif not postData:
                finalData='pre'
            if final is None:
                if start<=float(pre[4]):
                    gene=pre[-1].split(';')[2].split(' ')[-1].replace("\"","")
                    information['Inside']='True'
                    information['DistanceToGene']=0
                    information['Direction']='N/A'
                    focusGene=pre
                elif stop>=float(post[3]):
                    gene=post[-1].split(';')[2].split(' ')[-1].replace("\"","")
                    information['Inside']='True'
                    information['DistanceToGene']=0
                    information['Direction']='N/A'
                    focusGene=post
                else:
                    distance1=start-float(pre[4])
                    distance2=float(post[3])-stop
                    if distance1<distanceThreshold and distance1<distance2:
                        gene=pre[-1].split(';')[2].split(' ')[-1].replace("\"","")
                        information['DistanceToGene']=distance1
                        information['Direction']='Down stream'
                        information['Inside']='False'
                        focusGene=pre
                        
                    elif  distance2<distanceThreshold and distance2<distance1:
                        gene=post[-1].split(';')[2].split(' ')[-1].replace("\"","")
                        information['DistanceToGene']=distance2
                        information['Direction']='Up stream'
                        information['Inside']='False'
                        focusGene=post
                    else:
                        gene='N/A'
                        focusGene=None
                        information['DistanceToGene']='N/A'
                        information['Direction']='N/A'
                        information['Inside']='N/A'    
                information['Gene']=gene
            else:
                if final=='pre':
                    gene=pre[-1].split(';')[2].split(' ')[-1].replace("\"","")
                    if start<=float(pre[4]):
                        information['Inside']='True'
                        information['DistanceToGene']=0
                        information['Direction']='N/A'
                        distance1=0
                    else:
                        distance1=start-float(pre[4])
                        information['DistanceToGene']=distance1
                        information['Direction']='Down stream'
                        information['Inside']='False'
                    focusGene=pre     
                else:
                    gene=post[-1].split(';')[2].split(' ')[-1].replace("\"","")
                    if stop>=float(post[3]):
                        information['Inside']='True'
                        information['DistanceToGene']=0
                        information['Direction']='N/A'
                        distance1=0
                    else:
                        distance1=float(post[3])-stop
                        information['DistanceToGene']=distance1
                        information['Direction']='Up stream'
                        information['Inside']='False'
                    focusGene=post
                if distance1<distanceThreshold:
                    information['Gene']=gene
                else:
                    information['Gene']='N/A'
            if finalData is None:
                if start<=float(preData[4]):
                    focus=preData
                elif stop>=float(postData[3]):
                    focus=postData
                else:
                    focus=None
            else:
                if finalData=='pre':
                    if start<=float(preData[4]):
                        focus=preData
                    else:
                        focus=None
                elif finalData=='post':
                    if stop>=float(postData[3]):
                        focus=postData
                    else:
                        focus=None
            if focus is not None and focus[2]!='gene':
                information['Focus']=focus[2]
            else:
                information['Focus']='N/A'
        return [information,focusGene,focus]
                    
                    
                    
    def printToLog(self,message):
        self.core.printToLog(message)
        

    
#def getViralReferenceFa():
#    return
#
#def getHostReferenceFa():
#    return
#    
#def getFastQReadsLocation():
#    return
#    
#def installDependencies():
#    
#def checkDependencies():
#    
#def checkIndex():
#

    

#=============================================================================================================================================================
#=============================================================================================================================================================
#=============================================================================================================================================================



class Sam:

    Quit=False
    def __init__(self,filename,CS,readBaseName):
        self.sam=[]
        self.readBaseName=readBaseName
        self.runDirectory=os.path.dirname(filename)
        self.storage=[]
        self.header=[]
        self.flags=[]
        self.minLen=10
        self.loadSam(filename)
        if self.Quit==False:
            if CS:
                self.loadFlagsCSAlignment()
            else:
                self.loadFlags()
                
    def loadSam(self,filename):
        file = open(filename,'r')
        for line in file:
            if line.startswith('@'):
                self.header.append(line.split())
            else:
                self.sam.append(line.split())
        file.close()
        if len(self.sam)==0:
            self.Quit=True
        
        
    def loadFlags(self):
        countPaired=0
        countAligned=0
        countUnaligned=0
        
        for i in range (0,len(self.sam)):
            self.flags.append([])
            flag=bin(int(self.sam[i][1]))[2:]
            for j in range (len(flag),0,-1):
                self.flags[i].append(flag[j-1])
            for k in range (0,len(self.flags[i])):
                if k==0 and self.flags[i][k]=='1':
                    countPaired+=1
                if k==1 and self.flags[i][k]=='1':
                    countAligned+=1
                if k==2 and self.flags[i][k]=='1':
                    countUnaligned+=1
        self.alignmentRate=((len(self.sam)-countUnaligned)/len(self.sam))
                    
    
    def loadFlagsCSAlignment(self):
        countAligned=0
        countUnaligned=0
        for i in range (0,len(self.sam)):
            self.flags.append([])
            flag=bin(int(self.sam[i][1]))[2:]
            for j in range (len(flag),0,-1):
                self.flags[i].append(flag[j-1])
            for k in range (0,len(self.flags[i])):
                if k==1 and self.flags[i][k]=='1':
                    countAligned+=1
                if k==2 and self.flags[i][k]=='1':
                    countUnaligned+=1
        self.alignmentRate=((len(self.sam)-countUnaligned)/len(self.sam))
                
    def extractClippedSeq(self,minLen,paired):
        count=0
        out=open(self.runDirectory+'/'+self.readBaseName+'_CSList.txt','w')
        out2=open(self.runDirectory+'/'+self.readBaseName+'_CSListFull.txt','w')
        for i in range (0,len(self.sam)):
            cigar=self.sam[i][5]
            if cigar!="*" and 'S' in cigar:
                pairs=re.findall('\d+[A-Z]+',cigar)
                index=0
                for j in range (0,len(pairs)):
                    length=int(pairs[j][0:len(pairs[j])-1])
                    if pairs[j].endswith('S') and length>=minLen:
                        seq=Seq(self.sam[i][9][index:index+length])
                        fullseq=Seq(self.sam[i][9])
                        rID=self.sam[i][0]
                        if paired:
                            pair=2
                            if self.flags[i][6]=='1':
                                pair=1
                            if self.flags[i][4]=='0':
                                self.storage.append([rID,i,seq])
                                out.write('>'+str(i)+'.'+str(rID)+'.'+str(pair)+'\n')
                                out.write(str(seq)+'\n')
                                out2.write('>'+str(i)+'.'+str(rID)+'.'+str(pair)+'\n')
                                out2.write(str(fullseq)+'\n')
                            else:
                                self.storage.append([rID,i,seq,i])
                                out.write('>'+str(i)+'.'+str(rID)+'.'+str(pair)+'\n')
                                out.write(str(seq.reverse_complement())+'\n')
                                out2.write('>'+str(i)+'.'+str(rID)+'.'+str(pair)+'\n')
                                out2.write(str(fullseq)+'\n') #.reverse_complement()
                        else:
                            self.storage.append([rID,i,seq])
                            out.write('>'+str(i)+'.'+str(rID)+'\n')
                            out.write(str(seq)+'\n')
                            out2.write('>'+str(i)+'.'+str(rID)+'\n')
                            out2.write(str(fullseq)+'\n')
                        count+=1
                        break
                    index+=length
        
        out.close()
        out2.close()
        
                
    def returnAlignedCS(self,s,mode):
        for i in range (0,len(self.sam)):
            if int(self.sam[i][1]) != 4:
                sid=int(self.sam[i][0].split('.')[0])
                
                

#=============================================================================================================================================================
#=============================================================================================================================================================
#=============================================================================================================================================================

def hello():
    pass





class Interface(Frame):
    pairedReads=True
    filtered=False
    useBasic=False
    filteredReroute=[]
    releaseVersion = "1.13"
    currentConfigEntry=None
    selectSplitFiles_Window=0
    #selectDirVar = tki.IntVar()
    platformType=platform.system()
    def __init__(self, master,c):
        Frame.__init__(self, master)
        self.c=c
        self.defaultbg=master.cget('bg')
        self.create_widgets()
        self.master=master
        
    def create_widgets(self):
        self.master.title("JBS ChimericSeq "+ self.releaseVersion)

        self.label1=Text(self.master,height=1,width=20,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.label1,'Fastq File(s):')
        self.label1.grid(row=0,column=0,columnspan=1,sticky=W)
        self.selectDirVar = IntVar()
        self.selectDirButton = Checkbutton(self.master, text = "Select Directory", variable=self.selectDirVar, onvalue=1, offvalue=0, command=self.selectDirToggle)
        self.selectDirButton.grid(row=0,column=1,columnspan=1,sticky=W)
        
        self.frame1=Frame(self.master)
        self.frame1.grid(row=1,column=0,sticky='nw',columnspan=2)         
        self.bttn1=Button(self.frame1,text = '. . .', command=self.getReads)
        self.bttn1.grid(row=0,column=0,sticky='',ipadx=3,ipady=5)
        self.fastqBox=Text(self.frame1,height=2,width=70,state=DISABLED,background=self.defaultbg)
        self.fastqBox.bind('<Control-F1>',self.selectSplitFiles)
        self.fastqBox.grid(row=0,column=1,sticky=W,columnspan=1)
        self.scrollbarFile=ttk.Scrollbar(self.frame1,command=self.fastqBox.yview)
        self.fastqBox.config(yscrollcommand=self.scrollbarFile.set)
        self.scrollbarFile.grid(row=0,column=2,sticky='ns')
      
        self.label2=Text(self.master,height=2,width=65,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.label2,'Reads Containing\nChimeric Sequences:')
        self.label2.grid(row=2,column=0,columnspan=2,sticky=W)

        self.listBoxFrame=Frame(self.master)
        self.listBoxFrame.grid(row=3,column=0,columnspan=1,rowspan=1,sticky='wens')
        self.listBox1=Listbox(self.listBoxFrame,selectmode=SINGLE,height=20,width=25,activestyle=NONE)
        self.listBox1.grid(row=0,column=0,sticky='wens')
#        self.listBox1.insert(END, "a list entry dfsgdfg dfg dfg dfg dfg dfg dfg dfg dfg df ")        
#        for item in ["one", "two", "three", "four","one", "two", "three", "four","one", "two", "three", "four","one", "two", "three", "four"]:
#            self.listBox1.insert(END, item)
        self.scrollbar1=ttk.Scrollbar(self.listBoxFrame,command=self.listBox1.yview)
        self.scrollbar2=ttk.Scrollbar(self.listBoxFrame,command=self.listBox1.xview,orient=HORIZONTAL)
        self.listBox1.config(xscrollcommand=self.scrollbar2.set)
        self.listBox1.config(yscrollcommand=self.scrollbar1.set)
        self.list1Selection=0
        self.listBox1.select_set(self.list1Selection)
        self.scrollbar1.grid(row=0,column=1,sticky='wns')
        self.scrollbar2.grid(row=1,column=0,sticky='ew')
        self.listBox1.bind('<Button-1>',self.listChoiceClick)
        self.listBox1.bind('<Up>',self.listChoice)
        self.listBox1.bind('<Down>',self.listChoice)
        
        self.statsFrame=Frame(self.master)
        self.statsFrame.grid(row=3,column=1,sticky='nswe',rowspan=1)  
        
        self.statsTreeCols=['Attribute','Value']
        self.statsTree=ttk.Treeview(self.statsFrame,columns=self.statsTreeCols,show='headings',height=15)
        self.statsTree.bind("<Button-3>",self.popTreeMenu)
        self.ysb=ttk.Scrollbar(self.statsFrame,orient=VERTICAL, command= self.statsTree.yview)
        self.xsb=ttk.Scrollbar(self.statsFrame,orient=HORIZONTAL, command= self.statsTree.xview)
        self.statsTree.column(0,width=170)
        self.statsTree.column(1,width=250)
        self.statsTree['yscroll'] = self.ysb.set
        self.statsTree['xscroll'] = self.xsb.set
        self.statsTree.grid(row=0,column=0,sticky='ns')
        self.ysb.grid(row=0,column=1,sticky='ns')
        self.xsb.grid(row=1,column=0,sticky='ew')
        for c in self.statsTreeCols:
            self.statsTree.heading(c,text=c.title())#command=lambda c=c: self._column_sort(c,True)
        self.attributeList=['{Read Name}','{Read Length}','{Fastq Source}','{Index}','{Viral Component Length}','{Local Viral Coordinates}',
        '{Viral Reference Coordinates}','{Viral Accession}','{Viral Orientation}',
        '{Viral TM (Basic)}','{Viral TM (Salt Adjusted)}',
        '{Host Component Length}','{Local Host Coordinates}','{Host Reference Coordinates}','{Chromosome}','{Gene}',
        '{Inside Gene}','{Distance To Gene}','{Direction}','{Focus Region}','{Host Orientation}',
        '{Host TM (Basic)}','{Host TM (Salt Adjusted)}','{Inserted Regions}','{Overlap}','{Overlap TM (Basic)}',
        '{Overlap TM (Salt Adjusted)}','{Host Map Quality}','{Viral Map Quality}','{Fastq score}']
        self.genericAttributeList=['{Max Mapping Score}','{Min Mapping Score}','{Aligned to Chromosome 1}','{Aligned to Chromosome 2}','{Aligned to Chromosome 3}',
        '{Aligned to Chromosome 4}','{Aligned to Chromosome 5}','{Aligned to Chromosome 6}','{Aligned to Chromosome 7}','{Aligned to Chromosome 8}','{Aligned to Chromosome 9}',
        '{Aligned to Chromosome 10}','{Aligned to Chromosome 11}','{Aligned to Chromosome 12}','{Aligned to Chromosome 13}','{Aligned to Chromosome 14}','{Aligned to Chromosome 15}',
        '{Aligned to Chromosome 16}','{Aligned to Chromosome 17}','{Aligned to Chromosome 18}','{Aligned to Chromosome 19}','{Aligned to Chromosome 20}','{Aligned to Chromosome 21}',
        '{Aligned to Chromosome 22}','{Aligned to Chromosome X}','{Aligned to Chromosome Y}',]
        for item in self.attributeList:
            self.statsTree.insert('','end',item,values=item)

                
        self.seqFrame=Frame(self.master)
        self.seqFrame.grid(row=0,column=2,ipadx=0,rowspan=3,columnspan=5,sticky='nw')
        self.seqLabel=Text(self.seqFrame,height=1,width=56,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.seqLabel.tag_configure("Virus", foreground="red")
        self.seqLabel.tag_configure('Host',background='yellow')
        self.seqLabel.tag_configure('Overlap',underline=True)
        self.changeLabel(self.seqLabel,'Sequence: Red=Viral Highlighted=Host Underlined=Overlap')
        self.seqLabel.tag_add('Virus','1.10','1.13')
        self.seqLabel.tag_add('Host','1.20','1.31')
        self.seqLabel.tag_add('Virus','1.37','1.47')
        self.seqLabel.tag_add('Host','1.37','1.47')
        self.seqLabel.tag_add('Overlap','1.37','1.47')
        self.seqLabel.grid(row=0,column=0,sticky='nw')
        self.seqBox=Text(self.seqFrame,height=4,width=70,state=DISABLED,relief=SUNKEN)
        self.seqBox.grid(row=1,column=0,sticky='nw')
        self.seqBox.tag_configure("Virus", foreground="red")
        self.seqBox.tag_configure('Host',background='yellow')
        self.seqBox.tag_configure('Overlap',underline=True)
        self.scrollbar3=ttk.Scrollbar(self.seqFrame,command=self.seqBox.yview)
        self.seqBox.config(yscrollcommand=self.scrollbar3.set)
        self.scrollbar3.grid(row=1,column=1,sticky='ns')
              
        self.logFrame=Frame(self.master)
        self.logFrame.grid(row=3,column=2,sticky='nwse',ipadx=0,ipady=0,rowspan=1,columnspan=5)
        self.log=Text(self.logFrame,height=20,width=70,state=DISABLED)
        self.log.tag_configure("Error", foreground="red")
        self.log.tag_configure("Important", background='yellow')
        self.log.grid(row=0,column=0,sticky='e')      
        self.scrollbar=ttk.Scrollbar(self.logFrame,command=self.log.yview)
        self.log.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.grid(row=0,column=1,sticky='nse')
        
        self.logoFrame=Frame(self.master)
        self.logoL=Text(self.logoFrame,height=1,width=45,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.logoL,'')
        self.logoL.grid(row=0,column=0,columnspan=1,sticky=W)
        self.logoCanvas=Canvas(self.logoFrame)
        self.logo=PhotoImage(file=self.c.workingDirectory+'/JBS.gif')
        self.logoCanvas.create_image(0,0,image=self.logo,anchor=NW)
        self.logoCanvas.grid(row=0,column=1)
        self.logoFrame.grid(row=4,column=3)  
        
        
        self.buttonFrame=Frame(self.master)
        self.buttonFrame.grid(row=4,column=0,sticky='nw',columnspan=2)
        self.buttonL=Text(self.buttonFrame,height=1,width=15,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.buttonL,'Controls:')
        self.buttonL.grid(row=0,column=0,columnspan=2,sticky=W)
        self.bttnStart = Button(self.buttonFrame, text = "Start Run", command=self.c.run)
        self.bttnStart.bind("<Button-3>",self.c.run1)
        self.bttnStart.grid(row=1,column=0,columnspan=2,sticky='nesw')
        self.yesbttn=Button(self.buttonFrame,text='Yes',command=hello)
        self.yesbttn.grid(row=2,column=0,sticky='nwes')
        self.nobttn=Button(self.buttonFrame,text='No',command=hello)
        self.nobttn.grid(row=2,column=1,sticky='nwes')
        
        
        self.modeIntVar=IntVar()
        self.menubar = Menu(self.master)
        self.filemenu = Menu(self.menubar, tearoff=0)
        self.modeMenu=Menu(self.filemenu,tearoff=0)   
        self.modeMenu.add_radiobutton(label='Integration',command=self.changeMode,variable=self.modeIntVar,value=0)
        self.modeMenu.add_radiobutton(label='Translocation',command=self.changeMode,variable=self.modeIntVar,value=1)
        self.modeIntVar.set(0)
        self.filemenu.add_cascade(label="Mode", menu=self.modeMenu)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Load Settings", command=self.c.loadSettings)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Save", command=self.c.saveData)
        self.filemenu.add_command(label='Save Settings', command=self.c.saveSettings)
        self.filemenu.add_separator()
        self.filemenu.add_command(label='Split Large File', command=self.c.splitLargeFile)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=self.EXIT)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.editmenu = Menu(self.menubar, tearoff=0)
        self.cpumenu = Menu(self.editmenu,tearoff=0)
        self.cpuIntVar=IntVar()
        self.autoModeVar=IntVar()
        self.autoModeVar.set(1)  ##set to 1 for final code
        self.loadAlignmentBoolVar=IntVar()
        self.loadAlignmentBoolVar.set(0)  ##set to 0 for final code
        self.editmenu.add_command(label="Set Locations", command=self.newWindow)
        self.editmenu.add_command(label='Configurations',command=self.optionsWindow)
        self.otempFilter=IntVar()
        self.otempFilter.set(0)
        self.olenFilter=IntVar()
        self.olenFilter.set(0)
        self.htempFilter=IntVar()
        self.htempFilter.set(0)
        self.hlenFilter=IntVar()
        self.hlenFilter.set(0)
        self.vtempFilter=IntVar()
        self.vtempFilter.set(0)
        self.vlenFilter=IntVar()
        self.vlenFilter.set(0)
        self.nomicrohomologyFilter=IntVar()### FSS
        self.nomicrohomologyFilter.set(0) ### FSS
        self.doGTFBoolVar=IntVar()
        self.doGTFBoolVar.set(1) # set to 1 for final code
        self.doReadCleanBoolVar=IntVar()
        self.doReadCleanBoolVar.set(0) # set to 1 for final code 
        self.filtermenu = Menu(self.editmenu,tearoff=0)
        self.editmenu.add_cascade(label='Filtering',menu=self.filtermenu)
        self.filtermenu.add_checkbutton(label='Overlap TM',command=lambda: self.changeFilter(1),onvalue=1,offvalue=0,variable=self.otempFilter)
        self.filtermenu.add_checkbutton(label='Overlap Length',command=lambda: self.changeFilter(2),onvalue=1,offvalue=0,variable=self.olenFilter)
        self.filtermenu.add_checkbutton(label='Host TM',command=lambda: self.changeFilter(3),onvalue=1,offvalue=0,variable=self.htempFilter)
        self.filtermenu.add_checkbutton(label='Host Length',command=lambda: self.changeFilter(4),onvalue=1,offvalue=0,variable=self.hlenFilter)
        self.filtermenu.add_checkbutton(label='Viral TM',command=lambda: self.changeFilter(5),onvalue=1,offvalue=0,variable=self.vtempFilter)
        self.filtermenu.add_checkbutton(label='Viral Length',command=lambda: self.changeFilter(6),onvalue=1,offvalue=0,variable=self.vlenFilter)
        ### FSS
        # self.filtermenu.add_checkbutton(label='No Microhomology',command=lambda: self.changeFilter(7),onvalue=1,offvalue=0,variable=self.nomicrohomologyFilter)
        
        self.filtermenu.add_separator()
        self.filtermenu.add_command(label='Reset Filters',command=self.resetFilters)
        self.editmenu.add_cascade(label='Thread Count',menu=self.cpumenu)
        self.editmenu.add_checkbutton(label='Prompting',command=self.changePrompting,onvalue=1,offvalue=0,variable=self.autoModeVar)
        #self.editmenu.add_checkbutton(label='Read cleaning',command=self.doReadClean,onvalue=1,offvalue=0,variable=self.doReadCleanBoolVar)
        self.editmenu.add_checkbutton(label='Load Alignments When Possible',command=self.loadAlignment,onvalue=1,offvalue=0,variable=self.loadAlignmentBoolVar)
        self.editmenu.add_checkbutton(label='Get Gene info',command=self.doGTF,onvalue=1,offvalue=0,variable=self.doGTFBoolVar)
        self.menubar.add_cascade(label="Options", menu=self.editmenu)
        self.helpmenu = Menu(self.menubar, tearoff=0)
        self.helpmenu.add_command(label="About", command=self.openWeb)
        self.menubar.add_cascade(label="Help", menu=self.helpmenu)
        self.master.config(menu=self.menubar)
        for i in self.c.Cores:
            self.cpumenu.add_radiobutton(label=i,command=self.changeCores,variable=self.cpuIntVar,value=i)
            temp=i
        if (temp > 1):
            temp = 2 # default 2 cores
        self.cpuIntVar.set(temp)
        self.c.coreOption=int(self.cpuIntVar.get())
        
        self.treeMenu=Menu(self.statsFrame,tearoff=0)
        self.treeMenu.add_command(label='Copy',command=self.copySelection)
                

        self.OptionsWin=tki.Toplevel(self.master)
        self.OptionsWin.title('Configurations')
        self.OptionsWin.geometry('390x360')   
        self.OptionsWin.wm_attributes("-topmost", 1)   
        self.OptionsWin.protocol('WM_DELETE_WINDOW',self.OptionsWin.withdraw)
        self.OptionsWin.withdraw()

        """
        self.readminlength=Text(self.OptionsWin,height=1,width=40,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.readminlength.grid(row=0,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.readminlength,'Read Min Length Preprocessing (0 - 100):')
        self.readminlength=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.readminlength.bind('<Return>',self.noReturn)
        self.readminlength.bind('<space>',self.noReturn)
        self.readminlength.bind('<FocusIn>',self.entry_Enter)
        self.readminlength.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,100, 10))
        self.readminlength.insert(END,self.c.preprocessReadLength)
        self.readminlength.grid(row=0,column=1,columnspan=1,sticky=W)
        
        self.maxunknown=Text(self.OptionsWin,height=1,width=40,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.maxunknown.grid(row=1,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.maxunknown,'Max % Unknown Nucleotides (0 - 100):')
        self.maxunknown=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.maxunknown.bind('<Return>',self.noReturn)
        self.maxunknown.bind('<space>',self.noReturn)
        self.maxunknown.bind('<FocusIn>',self.entry_Enter)
        self.maxunknown.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,100, 100))
        self.maxunknown.insert(END,self.c.preprocessUnknownPercent)
        self.maxunknown.grid(row=1,column=1,columnspan=1,sticky=W)
        """        
        
        self.clipLabel=Text(self.OptionsWin,height=1,width=40,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.clipLabel.grid(row=0,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.clipLabel,'Clipped Sequence Min Length (0 - 50):')
        self.clipText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.clipText.bind('<Return>',self.noReturn)
        self.clipText.bind('<space>',self.noReturn)
        self.clipText.bind('<FocusIn>',self.entry_Enter)
        self.clipText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,50, 10))
        self.clipText.insert(END,self.c.CSminLen)
        self.clipText.grid(row=0,column=1,columnspan=1,sticky=W)
        
        self.saltLabel=Text(self.OptionsWin,height=1,width=40,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.saltLabel.grid(row=1,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.saltLabel,'Salt Concentration (mM) (0 - 500):')
        self.saltText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.saltText.bind('<Return>',self.noReturn)
        self.saltText.bind('<space>',self.noReturn)
        self.saltText.bind('<FocusIn>',self.entry_Enter)
        self.saltText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,500, 115))
        self.saltText.insert(END,self.c.Na)
        self.saltText.grid(row=1,column=1,columnspan=1,sticky=W)
        
        self.TrimmingLabel=Text(self.OptionsWin,height=1,width=40,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.TrimmingLabel.grid(row=2,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.TrimmingLabel,'Trimming 5\' (bp) (0 - 70):')
        self.TrimmingText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.TrimmingText.bind('<Return>',self.noReturn)
        self.TrimmingText.bind('<space>',self.noReturn)
        self.TrimmingText.bind('<FocusIn>',self.entry_Enter)
        self.TrimmingText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,70, 0))
        self.TrimmingText.insert(END,self.c.Trim)
        self.TrimmingText.grid(row=2,column=1,columnspan=1,sticky=W)
        
        self.TrimmingLabel2=Text(self.OptionsWin,height=1,width=40,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.TrimmingLabel2.grid(row=3,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.TrimmingLabel2,'Trimming 3\' (bp) (0 - 70):')
        self.TrimmingText2=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.TrimmingText2.bind('<Return>',self.noReturn)
        self.TrimmingText2.bind('<space>',self.noReturn)
        self.TrimmingText2.bind('<FocusIn>',self.entry_Enter)
        self.TrimmingText2.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,70, 0))
        self.TrimmingText2.insert(END,self.c.Trim2)
        self.TrimmingText2.grid(row=3,column=1,columnspan=1,sticky=W)

        self.DistanceLabel=Text(self.OptionsWin,height=1,width=40,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.DistanceLabel.grid(row=4,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.DistanceLabel,'Gene Distance Threshold (0 - 1000000):')
        self.DistanceText=Text(self.OptionsWin,height=1,width=7,state=NORMAL)
        self.DistanceText.bind('<Return>',self.noReturn)
        self.DistanceText.bind('<space>',self.noReturn)
        self.DistanceText.bind('<FocusIn>',self.entry_Enter)
        self.DistanceText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,1000000, 10000))
        self.DistanceText.insert(END,self.c.distanceThreshold)
        self.DistanceText.grid(row=4,column=1,columnspan=1,sticky=E+W)

        self.SimilarityLabel=Text(self.OptionsWin,height=1,width=40,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.SimilarityLabel.grid(row=5,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.SimilarityLabel,'Similarity Max% (80 - 100):')
        self.readSimilarityText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.readSimilarityText.bind('<Return>',self.noReturn)
        self.readSimilarityText.bind('<space>',self.noReturn)
        self.readSimilarityText.bind('<FocusIn>',self.entry_Enter)
        self.readSimilarityText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',80,100,95))
        self.readSimilarityText.insert(END,self.c.readSimilarity)
        self.readSimilarityText.grid(row=5,column=1,columnspan=1,sticky=W)
        
        self.StretchLabel=Text(self.OptionsWin,height=1,width=40,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.StretchLabel.grid(row=6,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.StretchLabel,'nt Stretch Count Max (5 - 30):')
        self.ntStretchCountText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.ntStretchCountText.bind('<Return>',self.noReturn)
        self.ntStretchCountText.bind('<space>',self.noReturn)
        self.ntStretchCountText.bind('<FocusIn>',self.entry_Enter)
        self.ntStretchCountText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',5,30,8))
        self.ntStretchCountText.insert(END,self.c.ntStretchCount)
        self.ntStretchCountText.grid(row=6,column=1,columnspan=1,sticky=W)

        self.UnidentifyLabel=Text(self.OptionsWin,height=1,width=40,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.UnidentifyLabel.grid(row=7,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.UnidentifyLabel,'nt Unidentified Count Max (5 - 25):')
        self.ntUnidentifiedCountText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.ntUnidentifiedCountText.bind('<Return>',self.noReturn)
        self.ntUnidentifiedCountText.bind('<space>',self.noReturn)
        self.ntUnidentifiedCountText.bind('<FocusIn>',self.entry_Enter)
        self.ntUnidentifiedCountText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',5,25,10))
        self.ntUnidentifiedCountText.insert(END,self.c.ntUnidentifiedCount)
        self.ntUnidentifiedCountText.grid(row=7,column=1,columnspan=1,sticky=W)
        
        
        self.mh_stringency=Text(self.OptionsWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.mh_stringency.grid(row=8,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.mh_stringency,'Microhomology Max% (80 - 100):')
        self.mh_stringency=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.mh_stringency.bind('<Return>',self.noReturn)
        self.mh_stringency.bind('<space>',self.noReturn)
        self.mh_stringency.bind('<FocusIn>',self.mh_stringency_Enter)
        self.mh_stringency.bind('<FocusOut>',self.mh_stringency_Leave)
        #self.mh_stringency.insert(END,self.c.microhomologystringency)
        self.mh_stringency.grid(row=8,column=1,columnspan=1,sticky=W)
        
        self.mh_stringency1=Text(self.OptionsWin,height=1,width=31,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.mh_stringency1.grid(row=9,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.mh_stringency1,'Microhomology nt (15 - 35):')
        self.mh_stringency1=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.mh_stringency1.bind('<Return>',self.noReturn)
        self.mh_stringency1.bind('<space>',self.noReturn)
        self.mh_stringency1.bind('<FocusIn>',self.mh_stringency1_Enter)
        self.mh_stringency1.bind('<FocusOut>',self.mh_stringency1_Leave)
        #self.mh_stringency1.insert(END,self.c.microhomologystringency1)
        self.mh_stringency1.grid(row=9,column=1,columnspan=1,sticky=W)
        
        self.printToLog('===== set microhomology')
        if (self.c.microhomologystringency_selection == 0):
            self.mh_stringency.insert(END,self.c.microhomologystringency)
            self.mh_stringency1.insert(END,"")
        else:
            self.mh_stringency1.insert(END,self.c.microhomologystringency1)
            self.mh_stringency.insert(END,"")
        
        ### FSS filtering section
        self.tempFilterLabel=Text(self.OptionsWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.tempFilterLabel.grid(row=10,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.tempFilterLabel,'Overlap TM Max (C) (0 - 95):')
        self.tempFilterText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.tempFilterText.bind('<Return>',self.noReturn)
        self.tempFilterText.bind('<space>',self.noReturn)
        self.tempFilterText.bind('<FocusIn>',self.entry_Enter)
        self.tempFilterText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,95, 70))
        self.tempFilterText.insert(END,self.c.otempFilter)
        self.tempFilterText.grid(row=10,column=1,columnspan=1,sticky=W)
        
        self.lenFilterLabel=Text(self.OptionsWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.lenFilterLabel.grid(row=11,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.lenFilterLabel,'Overlap Length Max (0 - 100):')
        self.lenFilterText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.lenFilterText.bind('<Return>',self.noReturn)
        self.lenFilterText.bind('<space>',self.noReturn)
        self.lenFilterText.bind('<FocusIn>',self.entry_Enter)
        self.lenFilterText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,100, 24))
        self.lenFilterText.insert(END,self.c.olenFilter)
        self.lenFilterText.grid(row=11,column=1,columnspan=1,sticky=W)
        
        self.htempFilterLabel=Text(self.OptionsWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.htempFilterLabel.grid(row=12,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.htempFilterLabel,'Host TM Min (C) (0 - 80):')
        self.htempFilterText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.htempFilterText.bind('<Return>',self.noReturn)
        self.htempFilterText.bind('<space>',self.noReturn)
        self.htempFilterText.bind('<FocusIn>',self.entry_Enter)
        self.htempFilterText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,80, 25))
        self.htempFilterText.insert(END,self.c.htempFilter)
        self.htempFilterText.grid(row=12,column=1,columnspan=1,sticky=W)
        
        self.hlenFilterLabel=Text(self.OptionsWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.hlenFilterLabel.grid(row=13,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.hlenFilterLabel,'Host Length Min (0 - 100):')
        self.hlenFilterText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.hlenFilterText.bind('<Return>',self.noReturn)
        self.hlenFilterText.bind('<space>',self.noReturn)
        self.hlenFilterText.bind('<FocusIn>',self.entry_Enter)
        self.hlenFilterText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,100, 10))
        self.hlenFilterText.insert(END,self.c.hlenFilter)
        self.hlenFilterText.grid(row=13,column=1,columnspan=1,sticky=W)
        
        self.vtempFilterLabel=Text(self.OptionsWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.vtempFilterLabel.grid(row=14,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.vtempFilterLabel,'Viral TM Min (C) (0 - 80):')
        self.vtempFilterText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.vtempFilterText.bind('<Return>',self.noReturn)
        self.vtempFilterText.bind('<FocusIn>',self.entry_Enter)
        self.vtempFilterText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,80, 25))
        self.vtempFilterText.bind('<space>',self.noReturn)
        self.vtempFilterText.insert(END,self.c.vtempFilter)
        self.vtempFilterText.grid(row=14,column=1,columnspan=1,sticky=W)
        
        self.vlenFilterLabel=Text(self.OptionsWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.vlenFilterLabel.grid(row=15,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.vlenFilterLabel,'Viral Length Min (0 - 100):')
        self.vlenFilterText=Text(self.OptionsWin,height=1,width=5,state=NORMAL)
        self.vlenFilterText.bind('<Return>',self.noReturn)
        self.vlenFilterText.bind('<space>',self.noReturn)
        self.vlenFilterText.bind('<FocusIn>',self.entry_Enter)
        self.vlenFilterText.bind('<FocusOut>',lambda event:self.entry_Leave('<FocusOut>',0,100, 10))
        self.vlenFilterText.insert(END,self.c.vlenFilter)
        self.vlenFilterText.grid(row=15,column=1,columnspan=1,sticky=W)
        
        self.useBasicLabel=Text(self.OptionsWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.useBasicLabel.grid(row=16,column=0,columnspan=1,sticky=W)
        self.changeLabel(self.useBasicLabel,'Filter with Basic Temperature:')
        self.useBasicBttn=Button(self.OptionsWin,text = "False", command=self.switchBasic)
        self.useBasicBttn.grid(row=16,column=1,columnspan=1,sticky=W+E)
        
        
        self.NewWin=tki.Toplevel(self.master)
        self.NewWin.title('Location Options')
        self.NewWin.geometry('620x490')
        self.NewWin.protocol('WM_DELETE_WINDOW',self.NewWin.withdraw)
        self.NewWin.withdraw()
        self.bt2locLabel=Text(self.NewWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.bt2locLabel,'Bowtie2 Directory:')
        self.bt2locLabel.grid(row=0,column=0,columnspan=2,sticky=W)
        self.bttn2=Button(self.NewWin,text = '. . .', command=self.setBT2Folder)
        self.bttn2.grid(row=1,column=0,sticky=W,ipadx=3)
        self.bt2Text=Text(self.NewWin,height=1,width=70,state=DISABLED,background=self.defaultbg)
        self.changeLabel(self.bt2Text,self.c.Bowtie2Folder)
        self.bt2Text.grid(row=1,column=1,sticky=W,columnspan=4)
        self.vrfLabel=Text(self.NewWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.vrfLabel,'Viral Reference *.fa:')
        self.vrfLabel.grid(row=2,column=0,columnspan=2,sticky=W)
        self.bttn3=Button(self.NewWin,text = '. . .', command=self.setViralRefFa)
        self.bttn3.grid(row=3,column=0,sticky=W,ipadx=3)
        self.vrf=Text(self.NewWin,height=1,width=70,state=DISABLED,background=self.defaultbg)
        self.changeLabel(self.vrf,self.c.ViralRefFa)
        self.vrf.grid(row=3,column=1,sticky=W,columnspan=4)
        self.vddLabel=Text(self.NewWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.vddLabel,'Viral Index Directory:')
        self.vddLabel.grid(row=4,column=0,columnspan=2,sticky=W)
        self.bttn4=Button(self.NewWin,text = '. . .', command=self.setViralRefFolder)
        self.bttn4.grid(row=5,column=0,sticky=W,ipadx=3)
        self.vdd=Text(self.NewWin,height=1,width=70,state=DISABLED,background=self.defaultbg)
        self.changeLabel(self.vdd,self.c.ViralRefFolder)
        self.vdd.grid(row=5,column=1,sticky=W,columnspan=4)
        self.vpLabel=Text(self.NewWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.vpLabel,'Viral Bowtie2 Prefix:')
        self.vpLabel.grid(row=6,column=0,columnspan=2,sticky=W)
        self.vp=Text(self.NewWin,height=1,width=70,state=NORMAL)
        self.vp.insert(INSERT,self.c.viralPrefix)
        self.vp.bind('<Return>',self.noReturn)
        self.vp.bind('<space>',self.noReturn)
        self.vp.grid(row=7,column=1,sticky=W,columnspan=4)
        self.vpbttn=Button(self.NewWin,text = 'Build', command=lambda: self.c.askbuildBowtie2ViralIndex(self.vp.get(1.0,'end-1c')))
        self.vpbttn.grid(row=7,column=0,sticky=W)
        self.hrfLabel=Text(self.NewWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.hrfLabel,'Host Reference *.fa:')
        self.hrfLabel.grid(row=8,column=0,columnspan=2,sticky=W)
        self.bttn5=Button(self.NewWin,text = '. . .', command=self.setHostRefFa)
        self.bttn5.grid(row=9,column=0,sticky=W,ipadx=3)
        self.hrf=Text(self.NewWin,height=1,width=70,state=DISABLED,background=self.defaultbg)
        self.changeLabel(self.hrf,self.c.HostRefFa)
        self.hrf.grid(row=9,column=1,sticky=W,columnspan=4)
        self.hddLabel=Text(self.NewWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.hddLabel,'Host Index Directory:')
        self.hddLabel.grid(row=10,column=0,columnspan=2,sticky=W)
        self.bttn6=Button(self.NewWin,text = '. . .', command=self.setHostRefFolder)
        self.bttn6.grid(row=11,column=0,sticky=W,ipadx=3)
        self.hdd=Text(self.NewWin,height=1,width=70,state=DISABLED,background=self.defaultbg)
        self.changeLabel(self.hdd,self.c.HostRefFolder)
        self.hdd.grid(row=11,column=1,sticky=W,columnspan=4)
        self.hpLabel=Text(self.NewWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.hpLabel,'Host Bowtie2 Prefix:')
        self.hpLabel.grid(row=12,column=0,columnspan=2,sticky=W)
        self.hp=Text(self.NewWin,height=1,width=70,state=NORMAL)
        self.hp.insert(INSERT,self.c.HostPrefix)
        self.hp.grid(row=13,column=1,sticky=W,columnspan=4)
        self.hp.bind('<Return>',self.noReturn)
        self.hp.bind('<space>',self.noReturn)
        self.hpbttn=Button(self.NewWin,text = 'Build', command=lambda: self.c.askbuildBowtie2HostIndex(self.hp.get(1.0,'end-1c')))
        self.hpbttn.grid(row=13,column=0,sticky=W)

        self.hgtfLabel=Text(self.NewWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)
        self.changeLabel(self.hgtfLabel,'Host GTF File:')
        self.hgtfLabel.grid(row=14,column=0,columnspan=2,sticky=W)
        self.hgtf=Text(self.NewWin,height=1,width=70,state=DISABLED,background=self.defaultbg)
        self.changeLabel(self.hgtf,'None')
        self.hgtf.grid(row=15,column=1,sticky=W,columnspan=4)
        self.hgtfbttn=Button(self.NewWin,text = '. . .', command=self.setHostGTF)
        self.hgtfbttn.grid(row=15,column=0,sticky=W,ipadx=3)
        
        self.outLabel=Text(self.NewWin,height=1,width=30,state=DISABLED,relief=FLAT,background=self.defaultbg)        
        self.changeLabel(self.outLabel,'Output Directory:')
        self.outLabel.grid(row=16,column=0,columnspan=2,sticky=W)
        self.bttn7=Button(self.NewWin,text = '. . .', command=self.setOutputDir)
        self.bttn7.grid(row=17,column=0,sticky=W,ipadx=3)
        self.out=Text(self.NewWin,height=1,width=70,state=DISABLED,background=self.defaultbg)
        self.changeLabel(self.out,self.c.outDir)
        self.out.grid(row=17,column=1,sticky=W,columnspan=4)

    def changeMode(self):
        self.c.runmode=int(self.modeIntVar.get())
        if self.c.runmode==0:
            self.printToLog('Switching to integration discovery mode')
        else:
            self.printToLog('Switching to translocation discovery mode')
    
    
    def popTreeMenu(self,event):
        iid=self.statsTree.identify_row(event.y)
        if iid:
            iid='\"'+iid+'\"'
            self.statsTree.selection_set(iid)
            self.treeMenu.post(event.x_root,event.y_root)
        
    def selectSplitFiles(self, event):
        if self.selectDirVar.get() == 1:
            self.printToLog('FSS == <Ctrl-F1> pressed!!!')
        else:
            return
        
        if (self.selectSplitFiles_Window != 0):
            self.selectSplitFiles_Window.destroy()
            self.selectSplitFiles_Window = 0

        self.selectSplitFiles_Window = Toplevel(self.master)
        self.selectSplitFiles_Window.wm_attributes("-topmost", 1)
        self.selectSplitFiles_Window.title("Select split files")
        self.selectSplitFiles_lb = Label(self.selectSplitFiles_Window)
        self.selectSplitFiles_lb["text"] = " Specify file number: n1, n2, ... "
        self.selectSplitFiles_lb.grid(row = 0, column = 0)
        # file numbers
        self.selectSplitFiles_Numbers = Entry(self.selectSplitFiles_Window, width=120)
        self.selectSplitFiles_Numbers.grid(row = 1, column = 0)
        
        self.selectSplitFiles_btn = Button(self.selectSplitFiles_Window, text=" OK ", \
                                              command=self.close_selectSplitFiles_Window)
        self.selectSplitFiles_btn.grid(row = 2, column = 0, sticky = "W"+"E"+"N"+"S")
        self.selectSplitFiles_Window.protocol('WM_DELETE_WINDOW', self.close_selectSplitFiles_Window)
        self.selectSplitFiles_Window.grab_set()
        return True
            
    def close_selectSplitFiles_Window(self):
        numbers=self.selectSplitFiles_Numbers.get().replace(" ", "")
        self.printToLog('FSS == numbers length = '+str(len(numbers)))
        if numbers == "":
            self.c.splitSelectedFilesUI=[]
        else:
            self.c.splitSelectedFilesUI=numbers.split(",")
        self.selectSplitFiles_Window.grab_release()
        self.selectSplitFiles_Window.destroy()
        self.selectSplitFiles_Window = 0
        return True
        
        return
    
    def copySelection(self):
        item = self.statsTree.focus()
        copy=(self.statsTree.item(item)['values'])
        self.c.topLevelTK.clipboard_clear()
        self.c.topLevelTK.clipboard_append(str(copy))
    
    def resetYesNoButtons(self):
        self.yesbttn.configure(command=hello)
        self.nobttn.configure(command=hello)
    
    def listChoiceClick(self,event):
        self.listBox1.select_clear(self.list1Selection)
        self.list1Selection=self.listBox1.nearest(event.y)
        self.listBox1.select_set(self.list1Selection)
        self.loadSelection()
    
    def doGTF(self):
        if self.doGTFBoolVar.get()==1:
            self.c.doGTF=True
            self.printToLog('Including gene characterization')
        else:
            self.c.doGTF=False
            self.printToLog('Excluding gene characterization')
    
    def doReadClean(self):
        if self.doReadCleanBoolVar.get()==1:
            self.c.doReadClean=True
            self.printToLog('Read Preprocessing Enabled')
        else:
            self.c.doReadClean=False
            self.printToLog('Read Preprocessing Disabled')
            
            
    def EXIT(self):
        self.master.destroy()
        sys.exit()
    
    def openWeb(self):
        webbrowser.open('http://www.jbs-science.com/ChimericSeq.php')
    
    def changeFilter(self,num):
        if num==1:
            if self.otempFilter.get()==1:
                self.printToLog('Overlap TM Filter: Enabled')
            else:
                self.printToLog('Overlap TM Filter: Disabled')
        elif num==2:
            if self.olenFilter.get()==1:
                self.printToLog('Overlap Length Filter: Enabled')
            else:
                self.printToLog('Overlap Length Filter: Disabled')
        elif num==3:
            if self.htempFilter.get()==1:
                self.printToLog('Host TM Filter: Enabled')
            else:
                self.printToLog('Host TM Filter: Disabled')
        elif num==4:
            if self.hlenFilter.get()==1:
                self.printToLog('Host Length Filter: Enabled')
            else:
                self.printToLog('Host Length Filter: Disabled')
        elif num==5:
            if self.vtempFilter.get()==1:
                self.printToLog('Viral TM Filter: Enabled')
            else:
                self.printToLog('Viral TM Filter: Disabled')
        elif num==6:
            if self.vlenFilter.get()==1:
                self.printToLog('Viral Length Filter: Enabled')
            else:
                self.printToLog('Viral Length Filter: Disabled')
        elif num==7:
            if self.nomicrohomologyFilter.get()==1:
                self.printToLog('No Microhomology Filter: Enabled')
            else:
                self.printToLog('No Microhomology Filter: Disabled')
        self.updateFilterValues()
        self.listRedraw()
    
    def updateFilterValues(self):
        self.c.olenFilter=float(self.lenFilterText.get(1.0,'end-1c'))
        self.c.otempFilter=float(self.tempFilterText.get(1.0,'end-1c'))
        self.c.hlenFilter=float(self.hlenFilterText.get(1.0,'end-1c'))
        self.c.htempFilter=float(self.htempFilterText.get(1.0,'end-1c'))
        self.c.vlenFilter=float(self.vlenFilterText.get(1.0,'end-1c'))
        self.c.vtempFilter=float(self.vtempFilterText.get(1.0,'end-1c'))
   
    def resetFilters(self):
        self.otempFilter.set(0)
        self.olenFilter.set(0)
        self.htempFilter.set(0)
        self.hlenFilter.set(0)
        self.vtempFilter.set(0)
        self.vlenFilter.set(0)
        self.nomicrohomologyFilter.set(0) ### FSS
        self.printToLog('Filters Reset')
        try:
            len(self.c.Map.data)
            for i in range(0, len(self.c.Map.data)):
                self.c.Map.data[i]['Filtered']=False
            self.listRedraw()
        except:
            pass
        
    def switchBasic(self):
        self.useBasic=not self.useBasic
        self.c.useBasic=self.useBasic
        if self.useBasic:
            self.changeBttnText(self.useBasicBttn,'True')
            self.printToLog('Using Basic Melt Temperatures for Filtering')
        else:
            self.changeBttnText(self.useBasicBttn,'False')
            self.printToLog('Using Salt Adjusted Melt Temperatures for Filtering')
        self.listRedraw()
            
    def selectDirToggle(self):
        self.changeLabel(self.fastqBox,"")
        self.c.resetLargeFileProcess()
        return
            
    def listChoice(self,event):
        #self.listBox1.yview_scroll(1,'units')
        if self.list1Selection<self.listBox1.size()-1 and event.keycode==40:
            self.listBox1.select_clear(self.list1Selection)
            self.list1Selection+=1
            self.listBox1.select_set(self.list1Selection)
        elif self.list1Selection>0 and event.keycode==38:
            self.listBox1.select_clear(self.list1Selection)
            self.list1Selection-=1
            self.listBox1.select_set(self.list1Selection)
        self.loadSelection()
    
    def loadSelection(self):  ##update to use property index instead of position in list
        try:
            self.c.Map
            carryon=True
        except:
            carryon=False
        if carryon:
            try:
                if self.filtered:
                    #index=self.filteredReroute[self.listBox1.curselection()[0]]##self.filteredReroute[self.c.Map.data[self.listBox1.curselection()[0]]['Index']]
                    index=self.filteredReroute[self.c.Map.data[self.listBox1.curselection()[0]]['Index']]
                else:
                    #index=self.listBox1.curselection()[0]##self.c.Map.data[self.listBox1.curselection()[0]]['Index']
                    index=self.filteredReroute[self.c.Map.data[self.listBox1.curselection()[0]]['Index']]
                #print self.c.Map.data.keys()[self.c.Map.data.values().index()
                self.statsTree.set('{Read Name}',column='Value',value=self.c.Map.data[index]['ReadName'])
                self.statsTree.set('{Read Length}',column='Value',value=self.c.Map.data[index]['ReadLength'])
                self.statsTree.set('{Fastq Source}',column='Value',value=self.c.Map.data[index]['FastqSource'])
                self.statsTree.set('{Index}',column='Value',value=self.c.Map.data[index]['Index'])
                self.statsTree.set('{Viral Component Length}',column='Value',value=self.c.Map.data[index]['Vlength'])
                self.statsTree.set('{Host Component Length}',column='Value',value=self.c.Map.data[index]['Hlength'])
                self.statsTree.set('{Local Viral Coordinates}',column='Value',value=self.c.Map.data[index]['ViralLocalCords'])
                self.statsTree.set('{Local Host Coordinates}',column='Value',value=self.c.Map.data[index]['HostLocalCords'])
                self.statsTree.set('{Chromosome}',column='Value',value=self.c.Map.data[index]['Chromosome'])
                
                self.statsTree.set('{Overlap}',column='Value',value=self.c.Map.data[index]['Overlap'])
                self.changeLabel(self.seqBox,self.c.Map.data[index]['Sequence'])
                self.statsTree.set('{Viral Reference Coordinates}',column='Value',value=self.c.Map.data[index]['ViralRefCords'])
                self.statsTree.set('{Host Reference Coordinates}',column='Value',value=self.c.Map.data[index]['HostRefCords'])
                self.statsTree.set('{Inserted Regions}',column='Value',value=self.c.Map.data[index]['Inserted'])
                self.statsTree.set('{Overlap TM (Basic)}',column='Value',value=str(self.c.Map.data[index]['OverlapTM'])+'°C')
                self.statsTree.set('{Overlap TM (Salt Adjusted)}',column='Value',value=str(self.c.Map.data[index]['OverlapTMAdjusted'])+'°C')
                self.statsTree.set('{Viral Accession}',column='Value',value=self.c.Map.data[index]['ViralAccession'])
                self.statsTree.set('{Viral TM (Basic)}',column='Value',value=str(self.c.Map.data[index]['VTM'])+'°C')
                self.statsTree.set('{Host TM (Basic)}',column='Value',value=str(self.c.Map.data[index]['HTM'])+'°C')
                self.statsTree.set('{Viral TM (Salt Adjusted)}',column='Value',value=str(self.c.Map.data[index]['VTMAdjusted'])+'°C')
                self.statsTree.set('{Host TM (Salt Adjusted)}',column='Value',value=str(self.c.Map.data[index]['HTMAdjusted'])+'°C')
                self.statsTree.set('{Viral Orientation}',column='Value',value=self.c.Map.data[index]['ViralOrientation'])
                self.statsTree.set('{Host Orientation}',column='Value',value=self.c.Map.data[index]['HostOrientation'])
                self.statsTree.set('{Host Map Quality}',column='Value',value=self.c.Map.data[index]['HMapQ'])
                self.statsTree.set('{Viral Map Quality}',column='Value',value=self.c.Map.data[index]['VMapQ'])
                self.statsTree.set('{Gene}',column='Value',value=self.c.Map.data[index]['Gene'])
                self.statsTree.set('{Inside Gene}',column='Value',value=self.c.Map.data[index]['InsideGene'])
                self.statsTree.set('{Direction}',column='Value',value=self.c.Map.data[index]['GeneDirection'])
                self.statsTree.set('{Distance To Gene}',column='Value',value=self.c.Map.data[index]['DistanceToGene'])
                self.statsTree.set('{Focus Region}',column='Value',value=self.c.Map.data[index]['Focus'])   
                viralIndexStart='1.'+str(self.c.Map.data[index]['ViralLocalCords'][0]-1)
                viralIndexStop='1.'+str(self.c.Map.data[index]['ViralLocalCords'][1])
                self.seqBox.tag_add('Virus',viralIndexStart, viralIndexStop)
                HostIndexStart='1.'+str(self.c.Map.data[index]['HostLocalCords'][0]-1)
                HostIndexStop='1.'+str(self.c.Map.data[index]['HostLocalCords'][1])
                self.seqBox.tag_add('Host',HostIndexStart, HostIndexStop)
                overlap = self.c.Map.data[index]['Overlap']
                if self.c.Map.data[index]['Overlap']>0:
                    if self.c.Map.data[index]['HostFirst']:
                        if  self.c.Map.data[index]['ViralLocalCords'][0]+overlap>= self.c.Map.data[index]['HostLocalCords'][1]:
                            overlapIndexStop=HostIndexStop
                            overlapIndexStart='1.'+str(self.c.Map.data[index]['HostLocalCords'][1]-self.c.Map.data[index]['Overlap'])
                        else:
                            overlapIndexStop=viralIndexStop
                            overlapIndexStart='1.'+str(self.c.Map.data[index]['ViralLocalCords'][1]-self.c.Map.data[index]['Overlap'])
                    else:
                        if self.c.Map.data[index]['HostLocalCords'][0]+overlap >= self.c.Map.data[index]['ViralLocalCords'][1]:
                            overlapIndexStop=viralIndexStop
                            overlapIndexStart='1.'+str(self.c.Map.data[index]['ViralLocalCords'][1]-self.c.Map.data[index]['Overlap'])
                        else:
                            overlapIndexStop=HostIndexStop
                            overlapIndexStart='1.'+str(self.c.Map.data[index]['HostLocalCords'][1]-self.c.Map.data[index]['Overlap'])
                    self.seqBox.tag_add('Overlap',overlapIndexStart,overlapIndexStop)
            except IndexError:
                self.printToLog('Index Out of Range: Please Refresh List')
                
    def listRedraw(self):
        total_loaded = 0
        try:
            len(self.c.Map.data)
            GO=True
        except AttributeError:
            self.printToLog('Please do not apply filters until after data has been processed')
            self.resetFilters()
            GO=False
        stretchCount = int(self.ntStretchCountText.get(1.0,'end-1c'))
        stretchA=""
        stretchC=""
        stretchG=""
        stretchT=""
        for i in range(0, stretchCount):
            stretchA=stretchA+"A"
            stretchC=stretchC+"C"
            stretchG=stretchG+"G"
            stretchT=stretchT+"T"
        unidentifiedCount = int(self.ntUnidentifiedCountText.get(1.0,'end-1c'))
        if GO:
            self.listBox1.delete(0,END)
            self.filteredReroute=[]
            if self.otempFilter.get()==1 or self.olenFilter.get()==1 or self.htempFilter.get()==1 or self.hlenFilter.get()==1 \
               or self.vtempFilter.get()==1 or self.vlenFilter.get()==1:
                self.filtered=True
            else:
                self.filtered=False
            if self.filtered:
                for i in range(0,len(self.c.Map.data)):
                    stop=False
                    
                    # remove similarity
                    if self.c.Map.data[i]['Similarity']==True:
                        stop=True
                        
                    # remove microhomology
                    try:
                        if self.c.Map.data[i]['Microhomology']==True:
                            stop=True
                            # self.printToLog('remove Microhomology')
                    except ValueError:
                        stop=True

                    if self.otempFilter.get()==1 and (not stop):
                        try:
                            if self.useBasic:
                                temp=float(self.c.Map.data[i]['OverlapTM'])
                            else:
                                temp=float(self.c.Map.data[i]['OverlapTMAdjusted'])
                            if temp>float(self.c.otempFilter): ##switch < to a > for final code
                                stop=True
                        except ValueError:
                            stop=True
                    if self.olenFilter.get()==1 and (not stop):
                        try:
                            temp=float(self.c.Map.data[i]['Overlap'])
                            if temp>float(self.c.olenFilter):   ##switch < to a > for final code
                                stop=True
                        except ValueError:
                            stop=True
                    if self.htempFilter.get()==1 and (not stop):
                        try:
                            if self.useBasic:
                                temp=float(self.c.Map.data[i]['HTM'])
                            else:
                                temp=float(self.c.Map.data[i]['HTMAdjusted'])
                            if temp<float(self.c.htempFilter):
                                stop=True
                        except ValueError:
                            stop=True
                    if self.hlenFilter.get()==1 and (not stop):
                        try:
                            temp=float(self.c.Map.data[i]['Hlength'])
                            if temp<float(self.c.hlenFilter):
                                stop=True
                        except ValueError:
                            stop=True
                    if self.vtempFilter.get()==1 and (not stop):
                        try:
                            if self.useBasic:
                                temp=float(self.c.Map.data[i]['VTM'])
                            else:
                                temp=float(self.c.Map.data[i]['VTMAdjusted'])
                            if temp<float(self.c.vtempFilter):
                                stop=True
                        except ValueError:
                            stop=True
                    if self.vlenFilter.get()==1 and (not stop):
                        try:
                            temp=float(self.c.Map.data[i]['Vlength'])
                            if temp<float(self.c.vlenFilter):
                                stop=True
                        except ValueError:
                            stop=True
                    """        
                    if self.c.Map.data[i]['Similarity']==True:
                        stop=True
                    """
                    if (not stop) and (stretchA in self.c.Map.data[i]['Sequence']):
                        stop=True
                    if (not stop) and (stretchC in self.c.Map.data[i]['Sequence']):
                        stop=True
                    if (not stop) and (stretchG in self.c.Map.data[i]['Sequence']):
                        stop=True
                    if (not stop) and (stretchT in self.c.Map.data[i]['Sequence']):
                        stop=True

                    # check unidentifed nt
                    if (not stop):
                        [vstart,vstop] = self.c.Map.data[i]['ViralLocalCords']
                        [hstart,hstop] = self.c.Map.data[i]['HostLocalCords']
                        readLength = self.c.Map.data[i]['ReadLength']
                        if (vstop > hstop):
                            if ((readLength - vstop) > unidentifiedCount):
                                stop=True
                        else:
                            if ((readLength - hstop) > unidentifiedCount):
                                stop=True
                        if (not stop):
                            if (vstart > hstart):
                                if (hstart > unidentifiedCount):
                                    stop=True
                            else:
                                if (vstart > unidentifiedCount):
                                    stop=True
                    """
                    # remove microhomology
                    try:
                        if self.c.Map.data[i]['Microhomology']==True:
                            stop=True
                            # self.printToLog('remove Microhomology')
                    except ValueError:
                        stop=True
                    """
                    if stop:
                        self.c.Map.data[i]['Filtered']=True
                    else:
                        self.c.Map.data[i]['Filtered']=False
                    if not stop:
                        #print ('adding ',i,self.c.Map.data[i]['Index'],self.c.Map.data[i]['Overlap'])
                        self.filteredReroute.append(i)
                        name=self.c.Map.data[i]['ReadName']
                        self.listBox1.insert(END,name)
                        total_loaded +=1
            else:
                for i in range(0,len(self.c.Map.data)):
                    seq = self.c.Map.data[i]['Sequence']
                    #self.c.Map.data[i]['Filtered']=False
                    if (stretchA in seq) and (self.c.Map.data[i]['Filtered']==False):
                        self.c.Map.data[i]['Filtered']=True
                    if (stretchC in seq) and (self.c.Map.data[i]['Filtered']==False):
                        self.c.Map.data[i]['Filtered']=True
                    if (stretchG in seq) and (self.c.Map.data[i]['Filtered']==False):
                        self.c.Map.data[i]['Filtered']=True
                    if (stretchT in seq) and (self.c.Map.data[i]['Filtered']==False):
                        self.c.Map.data[i]['Filtered']=True

                    # check unidentifed nt
                    if (self.c.Map.data[i]['Filtered']==False):
                        [vstart,vstop] = self.c.Map.data[i]['ViralLocalCords']
                        [hstart,hstop] = self.c.Map.data[i]['HostLocalCords']
                        readLength = self.c.Map.data[i]['ReadLength']
                        if (vstop > hstop):
                            if ((readLength - vstop) > unidentifiedCount):
                                self.c.Map.data[i]['Filtered']=True
                        else:
                            if ((readLength - hstop) > unidentifiedCount):
                                self.c.Map.data[i]['Filtered']=True
                        if (self.c.Map.data[i]['Filtered']==False):
                            if (vstart > hstart):
                                if (hstart > unidentifiedCount):
                                    self.c.Map.data[i]['Filtered']=True
                            else:
                                if (vstart > unidentifiedCount):
                                    self.c.Map.data[i]['Filtered']=True
                        
                    name=self.c.Map.data[i]['ReadName']
                    if (self.c.Map.data[i]['Microhomology']!=True) and (self.c.Map.data[i]['Similarity']!=True) and (self.c.Map.data[i]['Filtered']!=True):
                        self.listBox1.insert(END,name)
                        self.filteredReroute.append(i)
                        total_loaded +=1
                    # list only unique & not microhomology
            self.printToLog(str(total_loaded)+" reads displayed")
    
    
    def loadAlignment(self):
        if self.loadAlignmentBoolVar.get()==1:
            self.printToLog('Loading prexisting alignments, if any')
            self.c.loadAlignments=True
        else:
            self.printToLog('Not loading prexisting alignments')
            self.c.loadAlignments=False
    def changePrompting(self):
        if self.autoModeVar.get()==1:
            self.printToLog('Prompting Enabled')
            self.c.autoMode=False
        else:
            self.printToLog('Prompting Disabled')
            self.c.autoMode=True
            
    
    def getReads(self):
        readNames=self.c.getReads()
        if readNames=='break':
            self.STOP=True
            return
        temp=''
        if self.c.folderMode==False and self.c.processLargeFile==False:
            for i in range (0,len(readNames)):
                if i<len(readNames)-1:
                    temp=temp+readNames[i]+'\n'
                else:
                    temp=temp+readNames[i]
            self.changeLabel(self.fastqBox,temp)
        else:
            self.changeLabel(self.fastqBox,self.c.processLargeFileDir)
        self.STOP=False
        
    def switchRunMode(self):
        if self.c.folderMode:
            self.editmenu.entryconfig(2,label='Directory Mode')
            self.c.folderMode=False
            self.changeLabel(self.label1,'Fastq File(s):')
            self.changeLabel(self.fastqBox,'')
            self.printToLog('Mode: Single Run')
        else:
            self.editmenu.entryconfig(2,label='Single Run Mode')
            self.c.folderMode=True
            self.changeLabel(self.label1,'Fastq Directory:')
            self.changeLabel(self.fastqBox,'')
            self.printToLog('Mode: Directory Run')
            
#    def switchAutoMode(self):            
#        if self.c.autoMode:
#            self.autoBttn.config(text='Turn prompting off')
#            self.c.autoMode=False
#            self.printToLog('Prompting OFF')
#        else:
#            self.autoBttn.config(text='Turn prompting on')
#            self.c.autoMode=True
#            self.printToLog('Prompting ON')
            
    def changeCores(self):
        self.c.coreOption=int(self.cpuIntVar.get())
        self.printToLog('Using '+str(self.cpuIntVar.get())+' CPU cores')
        
    def setBT2Folder(self):
        temp=filedialog.askdirectory(title='Select Bowtie2 Root Directory',initialdir=self.c.workingDirectory)
        if temp:
            self.c.Bowtie2Folder=temp
            self.printToLog('Bowtie2 Folder changed to:\n'+self.c.Bowtie2Folder)
            self.changeLabel(self.bt2Text,self.c.Bowtie2Folder)
    
    def setOutputDir(self):
        temp=filedialog.askdirectory(title='Select Output Directory',initialdir=self.c.workingDirectory)
        if temp:
            self.c.outDir=temp
            self.printToLog('Output Directory Selected:\n'+self.c.outDir)
            self.changeLabel(self.out,self.c.outDir)
            
    def setViralRefFa(self):
        if self.platformType == "Windows":
            temp=filedialog.askopenfilename(title='Select Viral Reference Fa',initialdir=self.c.ViralRefFolder,filetypes=(("Fasta Files", "*.fa;*.txt"),("All files", "*.*")))
        else:
            temp=filedialog.askopenfilename(title='Select Viral Reference Fa',initialdir=self.c.ViralRefFolder)            
        if temp:
            self.c.ViralRefFa=temp
            self.c.ViralRefFolder=os.path.dirname(self.c.ViralRefFa)
            self.changeLabel(self.vrf,self.c.ViralRefFa)
            self.changeLabel(self.vdd,self.c.ViralRefFolder)
            self.printToLog('Viral Reference .fa set to:\n'+self.c.ViralRefFa)
    
    def setHostGTF(self):
        if self.platformType == "Windows":        
            temp=filedialog.askopenfilename(title='Select Host GTF',initialdir=self.c.HostRefFolder,filetypes=(("GTF Files", "*.gtf;*.gff;*.txt"),("All files", "*.*")))
        else:
            temp=filedialog.askopenfilename(title='Select Host GTF',initialdir=self.c.HostRefFolder)
        if temp:
            self.c.HostGTF=temp
            self.changeLabel(self.hgtf,self.c.HostGTF)
            self.printToLog('Host GTF set to:\n'+self.c.HostGTF)
    
    def setHostRefFa(self):
        if self.platformType == "Windows":        
            temp=filedialog.askopenfilename(title='Select Host Reference Fa',initialdir=self.c.HostRefFolder,filetypes=(("Fasta Files", "*.fa;*.txt"),("All files", "*.*")))
        else:
            temp=filedialog.askopenfilename(title='Select Host Reference Fa',initialdir=self.c.HostRefFolder)
        if temp:
            self.c.HostRefFa=temp
            self.c.HostRefFolder=os.path.dirname(self.c.HostRefFa)    
            self.changeLabel(self.hrf,self.c.HostRefFa)
            self.changeLabel(self.hdd,self.c.HostRefFolder)
            self.printToLog('Host Reference .fa set to:\n'+self.c.HostRefFa) 
    
    def setHostRefFolder(self):
        temp=filedialog.askdirectory(title='Select Host Reference Folder',initialdir=self.c.workingDirectory)
        if temp:
            self.c.HostRefFolder=temp
            self.printToLog('Host Reference Folder set to:\n'+self.c.HostRefFolder)
            self.changeLabel(self.hdd,self.c.HostRefFolder)
    
    def setViralRefFolder(self):
        temp=filedialog.askdirectory(title='Select Viral Reference Folder',initialdir=self.c.workingDirectory)
        if temp:
            self.c.ViralRefFolder=temp
            self.printToLog('Viral Reference Folder set to:\n'+self.c.ViralRefFolder)
            self.changeLabel(self.vdd,self.c.ViralRefFolder)
        
        
    def changeLabel(self,label,message):
        label.configure(state=NORMAL)
        label.delete(1.0,END)
        label.insert(INSERT,message)
        label.configure(state=DISABLED)
        
    def changeBttnText(self,bttn,mytext):
        bttn.config(text=mytext)
        
        
    def changeLabelFree(self,label,message):
        label.configure(state=NORMAL)
        label.delete(1.0,END)
        label.insert(INSERT,message)
    def noReturn(self,x):
        return 'break'
    def entry_Enter(self, key):
        self.currentConfigEntry=self.OptionsWin.focus_get()
        return        
    def entry_Leave(self, key, minimum, maximum, average):
        if self.currentConfigEntry == None:
            return
        data=self.currentConfigEntry.get(1.0,'end-1c').strip()
        try:
            value=int(data)
            if value < minimum:
                self.currentConfigEntry.delete(1.0, END)
                self.currentConfigEntry.insert(1.0, str(minimum))            
            elif value > maximum:
                self.currentConfigEntry.delete(1.0, END)
                self.currentConfigEntry.insert(1.0, str(maximum))            
        except ValueError:
            self.currentConfigEntry.delete(1.0, END)
            self.currentConfigEntry.insert(1.0, str(average))            
        data=self.currentConfigEntry.get(1.0,'end-1c').strip()
        #self.printToLog('entry_Leave '+data+" ")
        return        
    def mh_stringency_Enter(self, x):
        self.currentConfigEntry=self.OptionsWin.focus_get()
        self.mh_stringency1.delete(1.0, END)
        self.mh_stringency1.insert(1.0, "")
        if self.mh_stringency.get(1.0,'end-1c').strip() == "":
            if self.c.microhomologystringency == "":
                self.mh_stringency.insert(1.0, self.c.microhomologystringency_avg)
            else:    
                self.mh_stringency.insert(1.0, self.c.microhomologystringency)
        else:
            self.mh_stringency.delete(1.0, END)
            self.mh_stringency.insert(1.0, self.c.microhomologystringency)
        return
    def mh_stringency_Leave(self, x):
        self.c.microhomologystringency_selection = 0
        try:
            value = int(self.mh_stringency.get(1.0,'end-1c'))
        except ValueError:
            value = self.c.microhomologystringency_avg
        self.c.microhomologystringency = value
        if (value < self.c.microhomologystringency_min):
            self.mh_stringency.delete(1.0, END)
            self.mh_stringency.insert(1.0, self.c.microhomologystringency_min)
            self.c.microhomologystringency = self.c.microhomologystringency_min
        if (value > self.c.microhomologystringency_max):
            self.mh_stringency.delete(1.0, END)
            self.mh_stringency.insert(1.0, self.c.microhomologystringency_max)            
            self.c.microhomologystringency = self.c.microhomologystringency_max
        return
    def mh_stringency1_Enter(self, x):
        self.currentConfigEntry=self.OptionsWin.focus_get()
        self.mh_stringency.delete(1.0, END)
        self.mh_stringency.insert(1.0, "")
        if self.mh_stringency1.get(1.0,'end-1c').strip() == "":
            if self.c.microhomologystringency1 == "":
                self.mh_stringency1.insert(1.0, self.c.microhomologystringency1_avg)
            else:    
                self.mh_stringency1.insert(1.0, self.c.microhomologystringency1)
        else:
            self.mh_stringency1.delete(1.0, END)
            self.mh_stringency1.insert(1.0, self.c.microhomologystringency1)
        return
    def mh_stringency1_Leave(self, x):
        self.c.microhomologystringency_selection = 1
        try:
            value = int(self.mh_stringency1.get(1.0,'end-1c'))
        except ValueError:
            value = self.c.microhomologystringency1_avg
        self.c.microhomologystringency1 = value
        if (value < self.c.microhomologystringency1_min):
            self.mh_stringency1.delete(1.0, END)
            self.mh_stringency1.insert(1.0, self.c.microhomologystringency1_min)
            self.c.microhomologystringency1 = self.c.microhomologystringency1_min
        if (value > self.c.microhomologystringency1_max):
            self.mh_stringency1.delete(1.0, END)
            self.mh_stringency1.insert(1.0, self.c.microhomologystringency1_max)
            self.c.microhomologystringency1 = self.c.microhomologystringency1_max
        return
    
    def newWindow(self):
#        NewWin=tki.Toplevel(self.master)
#        NewWin.title('Location Options')
#        NewWin.geometry('400x600')
#        NewWin.config(state='disable')
        self.NewWin.deiconify()
#        def quit_win():
#            NewWin.destroy()
#            NewWin.config(state='normal')       
#        NewWin.protocol('WM_DELETE_WINDOW',quit_win)
    
    def optionsWindow(self):
        self.OptionsWin.deiconify()
    

#    def toggle_text(self):
#    
#        if self.bttn3['text'] == "Hi":
#        # switch to Goodbye
#            self.bttn3["text"] = "Goodbye"
#        else:
#        # reset to Hi
#            self.bttn3["text"] = "Hi"   
#        
#        self.printToLog('hi')
        
    def printToLog(self,message):
        self.log.configure(state=NORMAL)
        self.log.insert(END,message+'\n')
        self.log.configure(state=DISABLED)
        self.log.see(END)
        self.c.logger.append(message)





c=Core()
