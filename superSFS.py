# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 17:12:08 2021

@author: Chauvet
"""

import pandas as pd
import os
import sys
import matplotlib.pyplot as plt 

def for_sfs(outg,thr,i,o):
    with open(i) as f:
        for l in f:
            if '##' in l:
                with open(o,'a') as f1:
                    f1.write(l)
                continue
            elif '#CHROM' in l:
                l_o=l.split()
                # print(l)
                og_list=[]
                with open(outg) as f2:
                    for og in f2:
                        og=og.strip()
                        if len(og)<2:
                            print('wrong out group data')
                        else:
                            og_list+=[l_o.index(og)]
                with open(o,'a') as f1:
                    f1.write(l)
            # print(l[9:-6])
            else:
                l_o=l.split()
                om=''
                for n in og_list:
                    om+=l_o[n]
                    c_a=om.count('1')
                if c_a>int(thr):
                    l_n=l_o[:3]+[l_o[4]]+[l_o[3]]+l_o[5:9]
                    l_t=l_o[9:]
                    for gt in l_t:
                        if gt[0]=='0':
                            gt_n='1'+'|'
                        else:
                            gt_n='0'+'|'
                        if gt[2]=='0':
                            gt_n+='1'
                        else:
                            gt_n+='0'
                        l_n+=[gt_n]
                    l_n='\t'.join(l_n)+'\n'
                    # print(l_n)
                else:
                    l_n=l
                with open(o,'a') as f1:
                    f1.write(l_n)
                    
        # print(og_list)


    
def get_anno(a):
    pf=pd.read_csv(a,sep='\t')
    g_l={}
    for index in pf.index:
        a=pf.iloc[index].values
        if a[1] not in g_l.keys():
            g_l[a[1]]=[a[0]]
        else:
            g_l[a[1]]=g_l[a[1]]+[a[0]]
    # print(g_l['Wine'])
    return g_l
    # print(g_l)

def get_index(i,anno):
    with open(i) as f:
        g_i={}
        for l in f:
            if '#CHROM' in l:
                l=l.split()
                for key in anno.keys():
                    for x in anno[key]:
                        if key not in g_i.keys():
                            g_i[key]=[l.index(x)]
                        else:
                            g_i[key]=g_i[key]+[l.index(x)]
                return g_i
                break

def count_derall(i2,g_i,anno_list,o):        
    with open(i2) as f:
        for l in f:
            if '#' in l:
                continue
            else:
                l=l.split()
                g_d={}
                for key in g_i.keys():
                    gt=''
                    for x in g_i[key]:
                        gt+=l[x]
                        n_d=gt.count('1')
                    g_d[key]=n_d
                # print(g_d)
                g_da=''
                with open(o,'a') as f1:
                    for x in anno_list:
                        g_da+=str(g_d[x])+'\t'
                    g_da=g_da.strip()+'\n'
                    # print(g_da)
                    f1.write(g_da)
                
                  

def ind_sfs(i,o,gp,start):
    pf=pd.read_csv(i,sep='\t')
    # print(pf.head())
    pf_c=pf.groupby(gp).count()
    # print(pf_c)
    num=[]
    cou=[]
    for index in pf_c.index:
        num+=[index]
        cou+=[pf_c.iloc[index].values[0]]
    gp_n=pd.DataFrame()
    gp_n['NUM']=num
    gp_n['COUNT']=cou
    t_all=gp_n['COUNT'].sum()
    gp_n['PER']=gp_n['COUNT']/t_all
    len_x=len(gp_n['NUM'][1:])
    if start==0:
        with open(o+'%s.txt' %gp, 'a') as f1:
            f1.write('1 observation'+'\n')
            n1a=''
            for n1 in num:
                n1a+='d0_'+str(n1)+'\t'
            f1.write(n1a.strip()+'\n')
            n2a=''
            for n2 in cou:
                n2a+=str(n2)+'\t'
            f1.write(n2a.strip()+'\n')
    # print(len_x)
    # print(gp_n['COUNT'].values)
    plt.figure(figsize=(50, 30), dpi=300)
    ax=plt.gca()
    ax.set_ylabel("PCT. sites", labelpad=15, font = {'family': 'serif',
                                # 'style': 'italic',
                                'weight': 'bold', 
                                'size': 60,
                                })
    ax.set_xlabel("No. derived allele copies",labelpad=15, fontdict = {'family': 'serif',
                                          # 'style': 'italic',
                                          'weight': 'bold', 
                                          'size': 60,
                                          })
    ax.set_title(gp, y=1.02, fontdict={'family': 'serif',
                                    'weight': 'bold',
                                    'size': 70})
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.spines['top'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    plt.xlim(start-1,len_x+1)
    plt.tick_params(which='major', length=30,labelsize=50,width=5,pad=10)
    plt.tick_params(which='minor',length=20,width=3)
    x_la=[start]
    for x_l in range(1,len_x//10+1):
        x_la+=[x_l*10]
    plt.xticks(x_la)
    plt.bar(gp_n['NUM'][start:],gp_n['PER'][start:])
    plt.savefig(o+'%s_{}.pdf'.format(start) %gp)
    

print('\nThis is a tool for speculation of ancestral allel, calculation of sfs and drawing its bar plot. It is easy-to-use and runing fast. What you should prepare is the phased vcf file containg the data of populations you intrested and the outgroup, the outgroup name file, and the annotation file. Enjoy it!!!\n')
print('It has four models:')
print('    0：Using all function, from original vcf data to sfs barplot')
print('    1: Only speculate the ancestral allel and output new vcf file using speculated allel as reference')
print('    2: Only count the frequency of derived allel in each snp of each population')
print('    3: Only draw bar polt of sfs using data generated from the results of calutation of sfs\n')
print('Example:')
print('    Model 0: python superSFS 0 ogdir threshold vcfdir annodir modir coutdir plotdir group')
print('    Model 1: python superSFS 1 ogdir threshold vcfdir outdir')
print('    Model 2: python superSFS 2 annodir modir coutdir')
print('    Model 3: python superSFS 3 coutdir plotdir group\n')
print('Explation for each parameter:')
print('    ogdir: direction of outgroup names file')
print('    threshold: a number that if the sum of variant allel in outpgroup greater than it,the variant allel will be counted as ancestral allel')
print('    vcfdir: direction of vcf data')
print('    vannodir: direction of annotation file with sample names in first column and group name in second colum. This file should has header in first row')
print('    vmodir: assign the output direction of generated vcf file using speculated allel as reference')
print('    countdir: assign the output direction of calculation of derived allels for each snp in each group')
print('    plotdir: assign the output direction of bar plot of sfs')
print('    group: the group that you want to analysis')

if len(sys.argv) > 2:
    opt=sys.argv[1]
    if str(opt)=='0':
        ogdir, threshold, vcfdir, annodir, modir, coutdir, plotdir, group = sys.argv[2:]
        for_sfs(ogdir,threshold,vcfdir,modir)
        anno=get_anno(annodir)
        g_i=get_index(modir,anno)
        anno_list=list(anno.keys())
        if not os.path.isfile(coutdir):
            with open(coutdir,'w') as f:
                name=''
                for x in anno_list:
                    name+=x+'\t'
                name=name.strip()+'\n'
                f.write(name)
        count_derall(modir,g_i,anno_list,coutdir)
        ind_sfs(coutdir,plotdir,group,0)
    
    elif str(opt)=='1':
        ogdir, threshold, vcfdir, outdir = sys.argv[2:]
        for_sfs(ogdir,threshold,vcfdir,outdir)
    
    elif str(opt)=='2':
        annodir, modir, coutdir = sys.argv[2:]
        anno=get_anno(annodir)
        g_i=get_index(modir,anno)
        anno_list=list(anno.keys())
        if not os.path.isfile(coutdir):
            with open(coutdir,'w') as f:
                name=''
                for x in anno_list:
                    name+=x+'\t'
                name=name.strip()+'\n'
                f.write(name)
        count_derall(modir,g_i,anno_list,coutdir)
    elif str(opt)=='3':
        coutdir, plotdir, group = sys.argv[2:]
        ind_sfs(coutdir,plotdir,group,0)


