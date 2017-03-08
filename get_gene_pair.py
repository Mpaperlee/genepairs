import os.path
import glob
import sys
def get_contig(GTF):
    with open(GTF) as IN:
        for i in IN.readlines():
            ELEMENT = i.strip().split('\t')
            contig, feature, direction, Tran_name = ELEMENT[0], ELEMENT[2], ELEMENT[6], ELEMENT[-1].split(';')[1].split('"')[-2]
            if os.path.isdir('./contig'):
                pass
            else:
                os.mkdir('./contig')
            FILEOUT = open('./contig/%s.txt' % contig,'a')
            if feature in ['start_codon','stop_codon']:
                FILEOUT.write(contig +'\t'+Tran_name+'\t'+feature+'\t'+ELEMENT[3]+'\t'+ELEMENT[4] +'\t' + direction + '\n')
def get_interval():
    allfile = glob.glob('./contig/*.txt')
    FILEOUT = open('results.txt','a')
    for i in allfile:
        if os.path.getsize(i):
            dic = {}
            DIC = {}
            print 'dealing with %s' % i
            with open(i) as IN:
                for i in IN.readlines():
                    ELEMENT = i.strip().split()
                    contig = ELEMENT[0]
                    if ELEMENT[-1] == '-':
                        if ELEMENT[1] not in dic:
                            dic[ELEMENT[1]] = [-int(ELEMENT[4])]
                        else:
                            dic[ELEMENT[1]].insert(0,-int(ELEMENT[3]))
                    if ELEMENT[-1] == '+':
                        if ELEMENT[1] not in dic:
                            dic[ELEMENT[1]] = [int(ELEMENT[3])]
                        else:
                            dic[ELEMENT[1]].append(int(ELEMENT[4]))
                            #print dic
                for i in dic.keys():
                    genename = i[:10]
                    if genename not in DIC:
                        DIC[genename] = [dic[i]]
                    else:
                        DIC[genename].append(dic[i])
                        #print DIC
                Dic  = {}
                for i in DIC:
                    if DIC[i][0][0] < 0 :
                        Dic[i] = sorted(DIC[i],key =lambda k: abs(DIC[i][0][0]))
                        DIC[i] = Dic[i][0]
                    else:
                        Dic[i] = sorted(DIC[i],key =lambda k: abs(DIC[i][0][-1]))
                        DIC[i] = Dic[i][-1]
                SORT =  sorted(DIC.items(), key=lambda d: abs(d[1][0]))   
                for k,v in zip(SORT[:-1],SORT[1:]):
                    if len(k[-1]) >1 and len(v[-1]) > 1:                
                        if k[-1][-1] * v[-1][0] < 0:
                            if k[-1][-1] < 0:
                                div = v[-1][0]-abs(k[-1][-1])
                                FILEOUT.write('\t'.join((contig,k[0],v[0],"conver","=> <=",str(div))) + '\n')
                            else:
                                div = abs(v[-1][0])-k[-1][-1]
                                FILEOUT.write('\t'.join((contig,k[0],v[0],"diver","<= =>",str(div))) + '\n')
                        else:pass
                    else:pass
        else:pass
    FILEOUT.close()
def gene_pairs_stastic():
    with open('results.txt') as IN:
        if not os.path.isdir('./conver'):
            os.mkdir('./conver')
        else:pass
        if not os.path.isdir('./diver'):
            os.mkdir('./diver')
        else:pass
        for i in IN.readlines():
            ALL = i.strip().split('\t')
            if ALL[-3] == "conver":
                with open('./conver/conver.txt','a') as OUT:
                    OUT.write(i)
                    if int(ALL[-1]) < 250:
                        with open('./conver/conver250.txt','a') as OUT:
                            OUT.write(i)
                    if int(ALL[-1]) < 500:
                        with open('./conver/conver500.txt','a') as OUT:
                            OUT.write(i)
                    if int(ALL[-1]) < 1000:
                        with open('./conver/conver1000.txt','a') as OUT:
                            OUT.write(i)
                    if int(ALL[-1]) > 1000:
                        with open('./conver/conver_over1000.txt','a') as OUT:
                            OUT.write(i)
                    else:pass
            else:
                with open('./diver/diver.txt','a') as OUT:
                    OUT.write(i)
                    if int(ALL[-1]) < 250:
                        with open('./diver/diver250.txt','a') as OUT:
                            OUT.write(i)
                    if int(ALL[-1]) < 500:
                        with open('./diver/diver500.txt','a') as OUT:
                            OUT.write(i)
                    if int(ALL[-1]) < 1000:
                        with open('./diver/diver1000.txt','a') as OUT:
                            OUT.write(i)
                    if int(ALL[-1]) > 1000:
                        with open('./diver/diver_over1000.txt','a') as OUT:
                            OUT.write(i)
                    else:pass
        os.popen('cat ./conver/conver1000.txt ./diver/diver1000.txt > results_below1000.txt')
if __name__=='__main__':
    gtf = sys.argv[1]
    get_contig(gtf)
    get_interval()
    gene_pairs_stastic()
    