#!/usr/bin/python


import sys, time, os, argparse, subprocess


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample", help="Sample Name", action = "store", required = True)
parser.add_argument("-l", "--library", help="Library", action = "store", required = True)
parser.add_argument("-d", "--date", help="Sequencing Date (YYMMDD format)", action = "store", required = True)
parser.add_argument("-f", "--fastq", help="fastq file (one arg per file)", action='append', required = True)
parser.add_argument("-r", "--reference", help="reference name and file", action='store', nargs=2, required = True, metavar=("name","file"))

parser.add_argument("--threads", help="Number of threads for bwa", type=int, action='store', default=1)
parser.add_argument("--seed", help="bwa mem seed", type=int, action='store', default = 19)
parser.add_argument("--matchScore", help="bwa mem match score", type=int, action='store', default = 1)
parser.add_argument("--misMatchPen", help="bwa mem mis-match penalty", type=int, action='store', default = 4)
parser.add_argument("--gapOpenPen", help="bwa mem gap open penalty", type=int, action='store', default = 6)
parser.add_argument("--gapExtPen", help="bwa mem gap extension penalty", type=int, action='store', default = 1)
parser.add_argument("--endClipPen", help="bwa mem end clipping penalty", type=int, action='store', default = 5)
parser.add_argument("--unpairedPen", help="bwa mem unpaired read penalty", type=int, action='store', default = 17)

parser.add_argument("--illumina", help="Quality format is Illumina 1.3+ (ASCII 64)", action='store_true')

parser.add_argument("--outDir", help="Output directory", action='store', default = "./")
parser.add_argument("--tmpDir", help="Temporary directory for Java memory overflow", action='store', default = "./")

parser.add_argument("--javaXmx", help="Heap sixe for Picard", action='store', default = "2g")

parser.add_argument("--subSample", help="proportion to subsample bam", action="store", required = False, type=float, default=1.0)

parser.add_argument("--flagRetain", help="keep reads with this bitwise flag", action="store", required = False, type=int)
parser.add_argument("--flagRemove", help="discard reads with this bitwise flag", action="store", required = False, type=int)

parser.add_argument("--minMAPQ", help="discard reads with mapping quality below INT", action="store", required = False, type=int)

parser.add_argument("--skipMapping", help="Skip mapping step", action="store_true", required = False)
parser.add_argument("--skipSort", help="Skip bam sort step", action="store_true", required = False)
parser.add_argument("--skipRmdup", help="Skip duplicate removal step", action="store_true", required = False)
parser.add_argument("--skipFlagstat", help="Skip samtools flagstat", action="store_true", required = False)


parser.add_argument("--cleanUpBams", help="Remove intermediate bam files", action="store_true", required = False)
parser.add_argument("--cleanUpLogs", help="Remove log files", action="store_true", required = False)
parser.add_argument("--test", help="Just print commands", action="store_true", required = False)

parser.add_argument("--runName", help="Run name", action="store", default = "")



args = parser.parse_args()

SM = args.sample
LB = args.library
DT = args.date

reads = args.fastq

refName, refFile = args.reference

if args.illumina:
    I = "-I"
else:
    I = ""

outDir = args.outDir + "/"
tmpDir = args.tmpDir + "/"


SAMPLE=".".join([SM,LB,DT])
prefix=".".join([SM,LB,DT,refName])
fileName = outDir + ".".join([prefix,"bwa",args.runName])

flagRetain = "" if not args.flagRetain else " -f " + str(args.flagRetain)
flagRemove = "" if not args.flagRemove else " -F " + str(args.flagRemove)
mapqArg = "" if not args.minMAPQ else " -q " + str(args.minMAPQ)

#############################################################################################################################

#mapping

if not args.skipMapping:
    bwaCommand = ["bwa mem", "-M", I, "-t", str(args.threads), "-R", "'@RG\\tID:" + SAMPLE + "\\tSM:" + SM + "'", "-k", str(args.seed), "-A", str(args.matchScore), "-B", str(args.misMatchPen), "-O", str(args.gapOpenPen), "-E", str(args.gapExtPen), "-L", str(args.endClipPen), "-U", str(args.unpairedPen), refFile, " ".join(reads), "2>", fileName + ".log", "| samtools view -bS -s " + str(args.subSample) + flagRetain + flagRemove + mapqArg + " - >", fileName + ".bam"]
    
    bwaCommand = " ".join(bwaCommand)
    
    sys.stderr.write("\nbwa command:\n" + bwaCommand + "\n")
    
    if not args.test:
        os.system(bwaCommand)

        #check for success
        if subprocess.check_output(["tail", "-1", fileName + ".log"]).split()[1] != b"Real":
            sys.stderr.write("\nbwa failed - check log file.\n")
            sys.exit()

        sys.stderr.write("\nbwa complete at {} {}\n".format(time.strftime("%d/%m/%Y"), time.strftime("%H:%M:%S")))



        if args.cleanUpLogs:
            os.remove(fileName + ".log")


#sort

if not args.skipSort:
    sortCommand = ["picard", "SortSam", "-Xmx"+args.javaXmx, "TMP_DIR=" + tmpDir, "MAX_RECORDS_IN_RAM=1000000", "SORT_ORDER=coordinate", "INPUT=" + fileName + ".bam", "OUTPUT=" + fileName + ".sort.bam", "2>", fileName + ".sort.log"]
    
    sortCommand = " ".join(sortCommand)
    
    sys.stderr.write("\nSort command:\n" + sortCommand + "\n")
    
    if not args.test:
        os.system(sortCommand)
        
        #check for success
        if subprocess.check_output(["tail", "-1", fileName + ".sort.log"]).split(b".")[0] != b"Runtime":
            sys.stderr.write("\nSort failed - check log file.\n")
            sys.exit()
        
        sys.stderr.write("\nSort complete at {} {}\n".format(time.strftime("%d/%m/%Y"), time.strftime("%H:%M:%S")))
        
        if args.cleanUpBams:
            os.remove(fileName + ".bam")
        
        if args.cleanUpLogs:
            os.remove(fileName + ".sort.log")


#flagstat

if not args.skipFlagstat:
    flagCommand = ["samtools", "flagstat",  fileName + ".sort.bam", ">", fileName + ".flagstat"]
    
    flagCommand = " ".join(flagCommand)
    
    sys.stderr.write("\nFlagstat command:\n" + flagCommand + "\n")
    
    if not args.test:
        os.system(flagCommand)
        
        sys.stderr.write("\nFlagstat complete at {} {}\n".format(time.strftime("%d/%m/%Y"), time.strftime("%H:%M:%S")))




#rmdup

if not args.skipRmdup:
    rmdupCommand = ["picard", "MarkDuplicates",  "-Xmx"+args.javaXmx, "TMP_DIR=" + tmpDir, "MAX_RECORDS_IN_RAM=1000000", "REMOVE_DUPLICATES=true", "ASSUME_SORTED=true", "VALIDATION_STRINGENCY=SILENT", "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000", "CREATE_INDEX=true", "INPUT=" + fileName + ".sort.bam", "OUTPUT=" + fileName + ".rmdup.bam", "METRICS_FILE=" + fileName + ".rmdup.metrics", "2>", fileName + ".rmdup.log"]
    
    rmdupCommand = " ".join(rmdupCommand)
    
    sys.stderr.write("\nRemove Duplicates command:\n" + rmdupCommand + "\n")
    
    if not args.test:
        os.system(rmdupCommand)
        
        #check for success
        if subprocess.check_output(["tail", "-1", fileName + ".rmdup.log"]).split(b".")[0] != b"Runtime":
            sys.stderr.write("\nRemove duplicates failed - check log file.\n")
            sys.exit()
        
        sys.stderr.write("\nRemove duplicates complete at {} {}\n".format(time.strftime("%d/%m/%Y"), time.strftime("%H:%M:%S")))
        
        if args.cleanUpBams:
            os.remove(fileName + ".sort.bam")
        
        if args.cleanUpLogs:
            os.remove(fileName + ".rmdup.log")
            os.remove(fileName + ".rmdup.metrics")

#flagstat after duplicate removal

if not args.skipFlagstat:
    flagCommand = ["samtools", "flagstat",  fileName + ".rmdup.bam", ">", fileName + ".rmdup.flagstat"]
    
    flagCommand = " ".join(flagCommand)
    
    sys.stderr.write("\nFlagstat command:\n" + flagCommand + "\n")
    
    if not args.test:
        os.system(flagCommand)

        sys.stderr.write("\nFlagstat complete at {} {}\n".format(time.strftime("%d/%m/%Y"), time.strftime("%H:%M:%S")))

