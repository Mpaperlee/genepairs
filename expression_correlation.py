# -*- coding: utf-8 -*- 
import sys
import os.path
from numpy import *
import pandas as pd
from itertools import islice
#合并一堆基因的表达数据，并计算出每个基因不通条件下的总和，方便后期筛选。
def get_dataframe(file1,file2):
    gene_pairs = pd.read_csv(file1,sep = '\t',index_col=False)
    fpkm       = pd.read_csv(file2,sep = '\t',index_col=False)
    df = pd.merge(gene_pairs,fpkm,on = 'name1')
    newlist = list(df)
    newlist.remove('contig')
    newlist.remove('distance')
    df['sum1'] = df[newlist].sum(axis=1)
    fpkm.rename(columns = {'name1':'name2'},inplace =True)
    df2 = pd.merge(df,fpkm,on = 'name2')
    newlist = list(df2)
    newlist.remove('contig')
    newlist.remove('distance')
    df2['sum2'] = df2[newlist].sum(axis=1) - 2*df2['sum1']
    FINNAL = df2[(df2['sum1'] + df2['sum1'])>2]
    FINNAL.to_csv('tmp.txt',sep = '\t',index = False)    
#通过sum来筛选出同时表达和单边表达的类型
def deal_withdataframe():
    with open('tmp.txt') as IN:
        title = IN.readline()
        with open('single_inexp.txt','a') as OUT:
            OUT.write(title)
        with open('all_inexp.txt','a') as OUT:
            OUT.write(title)
        for i in islice(IN, 1, None):
            ALL = i.strip().split('\t')
            sum1, sum2 = ALL[12], ALL[19]
            if float(sum1) == 0.0 or float(sum2) == 0.0:
                with open('single_inexp.txt','a') as OUT:
                    OUT.write(i)
            else:
                if float(sum1) >2 and float(sum2) >2 :
                    with open('all_inexp.txt','a') as OUT:
                        OUT.write(i)
#计算向量的夹角余弦
def get_cos():
    if os.path.exists('COS.txt'):
        os.popen('mv COS.txt COS.txt.old')
    with open('all_inexp.txt') as IN:
        title = IN.readline()
        with open('COS.txt','a') as OUT:
            OUT.write('\t'.join(['contig','name1','name2','type','type_detail','distance','cos'])+'\n')
            for i in islice(IN, 1, None):
                ALL = i.strip().split('\t')
                v1 = mat(map(float,ALL[6:12]))
                v2 = mat(map(float,ALL[13:19]))
                cos = float(v1*v2.T)/(linalg.norm(v1)*linalg.norm(v2))
                sim = 0.5 + 0.5*cos
                OUT.write('\t'.join((ALL[0:6]))+'\t'+str(sim)+'\n')

if __name__=='__main__':
    file1 ,file2 = sys.argv[1:]
    get_dataframe(file1, file2)
    deal_withdataframe()
#get_cos()