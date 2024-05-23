#!/usr/bin/env python
# coding: utf-8


import sys, subprocess, importlib, re
required_packages = ['jinja2==3.0.3', 'weasyprint==56.0']
for pack in required_packages:
    pack_name = pack.split('=')[0]
    if not importlib.util.find_spec(pack_name):
        stout = subprocess.check_call(f'{sys.executable} -m pip install {pack} -i https://pypi.tuna.tsinghua.edu.cn/simple', shell=True)
        print(f'{pack} successfully installed!')

import os, codecs, time
from argparse import ArgumentParser
#from subprogram import formattable, formatnumber, html2doc
from jinja2 import Environment, FileSystemLoader


def formatnumber(num):
    if re.match('^-?\d+\.\d+$', num):
        nums = re.search('^(-?\d+)\.(\d+)$', num).groups()
        if len(nums[1]) > 2:
            newnum = '%.2f' % float(num)
        else:
            newnum = num
    else:
        newnum = num
    return newnum

def formattable(file, format=False):
    if not format:
        newlines = []
        for line in  open(file):
            field = line.strip().split('\t')
            # field[0] = field[0].replace('*','')
            newlines.append(field)
    if format:
        newlines = []
        for line in  open(file):
            field = line.strip().split('\t')
            # field[0] = field[0].replace('*','')
            items = [formatnumber(i) for i in field]
            newlines.append(items)
        
    return(newlines)

def html2doc(htmlfile, docfile):
    import pypandoc
    output = pypandoc.convert_file(source_file=htmlfile, format='html', to='docx', outputfile=docfile, extra_args=['-RTS'])

paras = ArgumentParser(description='Neoantigen auto report program')
# paras.add_argument('-d', metavar='directory', help='Result directory for data summarizing.')
paras.add_argument('-t', metavar='TMB/TNB_file', help='File containing TMB/TNB calculating results.')
paras.add_argument('-a', metavar='HLA_file', help='File containing HLA analysis data.')
paras.add_argument('-n', metavar='Neo_file', help='File containing neoantigen analysis data.')
paras.add_argument('-p', metavar='project', default='project.txt', help='Tab delimited file containing project information. 1st line: project [project ID]; 2nd line: patient [patient ID]; 3rd line: species [species name]. DEFAULT=%(default)s.')
paras.add_argument('-o', metavar='outdir', default='./', help='Output directory for report.html and report.pdf. DEFAULT=[%(default)s].')
paras.add_argument('-c', metavar='cache_dir', help=' The directory containing font, img, src, and *html.')
args = paras.parse_args()

#dirname = os.path.dirname(os.path.realpath(__file__))
dirname = args.c
env = Environment(loader=FileSystemLoader(dirname), line_statement_prefix='#')
template = env.get_template('Neoantigen.html')

project = []
for line in open(args.p):
    field = line.strip().split('\t')
    project.append(field[1])
project.append(time.strftime('%Y-%m-%d', time.localtime()))


# TMB = open('fk01001_TP_TMB_TNB.txt').readlines()[-1].split('\t')[1]
# TNB = open('fk01001_TP_TMB_TNB.txt').readlines()[-1].split('\t')[2]
TMB = open(args.t).readlines()[-1].split('\t')[1]
TNB = open(args.t).readlines()[-1].split('\t')[2]


HLAtable = []
# for line in open('DASH.output.txt'):
for line in open(args.a):
    field = line.strip().split('\t')
    if line.startswith('Alleles'):
        items = ['Blood Sample', 'Tumor Sample', 'LOH']
    else:
        hla_type_ = field[0].split('_')
        hla_type = hla_type_[0].upper() + '-' + hla_type_[1].upper() + hla_type_[2] +':' + hla_type_[3]
        loh = 'Positive' if field[1] == 'True' else 'Negative'
        items = [hla_type, hla_type, loh]
    HLAtable.append(items)


# neohtml = formattable('pvacseq_fk01001_27aa_peptide_HLA_I_AF0.1.txt', format=False)
neohtml = formattable(args.n, format=False)
neotable = []
neotable.append(['HLA-I', 'Gene', 'MT Peptide', 'Allele Frequency', 'Expression (TPM)'])
for line in open(args.n):
    infor = []
    field = line.strip().split('\t')
    if line.startswith('chromosome'):
        header = field
        continue
    if len(field) < 2:
        continue
    HLA_infor = field[header.index('hla allele')]
    HLA_infor = HLA_infor.replace('*', '')
    infor.append(HLA_infor)
    infor.append(field[header.index('gene name')])
    infor.append(field[header.index('mt_peptide_sequence')])
    af_value = field[header.index('tumor dna vaf')]
    af_value = formatnumber(af_value)
    infor.append(af_value)
    tpm_value = field[header.index('gene expression')]
    tpm_value = formatnumber(tpm_value)
    infor.append(tpm_value)
    neotable.append(infor)
    
if not os.path.isdir(f'{args.o}/{project[1]}_report'):
    os.makedirs(f'{args.o}/{project[1]}_report')
    os.system(f'cp -r {dirname}/font {args.o}/{project[1]}_report')
    os.system(f'cp -r {dirname}/img {args.o}/{project[1]}_report')
    os.system(f'cp -r {dirname}/src {args.o}/{project[1]}_report')
fout = codecs.open(f'{args.o}/{project[1]}_report/{project[1]}_report.html', 'w', 'utf-8')
output = template.render(title='Neoantigen Sequencing', project=project, TMB=TMB, TNB=TNB, HLAtable=HLAtable, neotable=neotable, neohtml=neohtml)
fout.write(output)
fout.close()


os.system(f'weasyprint {args.o}/{project[1]}_report/{project[1]}_report.html {args.o}/{project[1]}_report/{project[1]}_report.pdf')

os.system(f'rm -rf {args.o}/{project[1]}_report/font')
os.system(f'rm -rf {args.o}/{project[1]}_report/img')
os.system(f'rm -rf {args.o}/{project[1]}_report/src')