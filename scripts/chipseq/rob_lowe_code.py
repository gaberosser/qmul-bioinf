import os
import gzip


def gff_to_tss_bed(fn, pad=2500):
    if os.path.splitext(fn)[1].lower() == '.gz':
        opener = gzip.open
    else:
        opener = open
    for k in opener(fn):
        if k.find("#")!=0:
            if k.find("protein_coding")!=-1 and k.split("\t")[2]=="gene":
                ids=k.split("ID=")[1].split(".")[0]
                if k.split("\t")[6]=="-":
                    print k.split("\t")[0].split("chr")[1]+"\t"+str(int(k.split("\t")[4])-pad)+"\t"+str(int(k.split("\t")[4])+pad)+"\t"+ids
                elif k.split("\t")[6]=="+":
                    print k.split("\t")[0].split("chr")[1]+"\t"+str(int(k.split("\t")[3])-pad)+"\t"+str(int(k.split("\t")[3])+pad)+"\t"+ids
                else:
                    print "ERROR",k
