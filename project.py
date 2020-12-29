#!/usr/bin/env python3
from __future__ import division

import gffutils
import os, itertools, re, gzip

def isheader(line):
	return line[0] == '>'

def aspairs(f):
	seq_id = ''
	sequence = ''
	for header,group in itertools.groupby(f, isheader):
		if header:
			line = next(group)
			seq_id = line[1:].split()[0]
		else:
			sequence = ''.join(line.strip() for line in group)
			yield seq_id, sequence

#url_gff = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/746/955/GCA_001746955.1_ASM174695v1/GCA_001746955.1_ASM174695v1_genomic.gff.gz"
url_fasta = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/746/955/GCA_001746955.1_ASM174695v1/GCA_001746955.1_ASM174695v1_genomic.fna.gz"

#file_gff = "GCA_001746955.1_ASM174695v1_genomic.gff.gz"
file_fasta = "GCA_001746955.1_ASM174695v1_genomic.fna.gz"

#if not os.path.exists(file_gff):
#	os.system("curl -O %s" %(url_gff))

if not os.path.exists(file_fasta):
	os.system("curl -O %s" %(url_fasta))

#Make a database from the chr1_gff3.gff3 file:
#force-True will rewrite the db
#merge_strategy="merge" will merge duplicate data
#from_string=True shows that the data is the actual data to use
with open("GS115_CRG.gff3", "rt") as file:
	data = file.read()
gffutils.create_db(data, dbfn='GS115_db.db', force=True, from_string=True,merge_strategy="merge")

#Defining a dictionary with chromosome IDs as dictionary keys and a list of dictionary values consisting of
#dict.value[0] = positive chromosome strand
#dict.value[1] = negative chromosome strand
chrs = {}
with gzip.open(file_fasta, "rt") as file:
	seqs = aspairs(file)
	j = 0
	for seq in seqs:
		if j == 4:
			break
		chr_id = seq[0]
		chr_str = seq[1]
#chr_str_for is the 5'- 3' strand:
		chr_str_for = chr_str.upper()
#chr_str_rev is the 3'- 5' strand:
		chr_str_r = (''.join(reversed(chr_str_for)))
		chr_str_rev = []
		for i in range(0, len(chr_str_r)):
			chr_str_rev.append(chr_str_r[i].translate(str.maketrans("ATCG", "TAGC")))
		chr_str_rev = (''.join(chr_str_rev))
		chr = []
		chr.append(chr_str_for)
		chr.append(chr_str_rev)
		chrs[str(chr_id)] = chr
		j += 1

#Choosing only "CDS" rows from the GFF file to avoid designing gRNAs located on introns:
db = gffutils.FeatureDB('GS115_db.db', keep_order=True)
CDS_features = db.features_of_type("CDS")
genes = list(CDS_features)

f=open("GS115_gRNA.txt","a+")
f.write("#This file contains all the gRNAs of K. phaffii GS115 strain that are located within the first 1/3 of each coding sequence of each gene.\n"\
"#the first set of gRNAs represent the gRNAs located on the positive strand of the chromosome, while the second set of gRNAs for each gene shows the occurence of all the gRNAs on the negative strand of the chromosome\n")

for chr_id in chrs.keys():
	print("%s\n" %(chr_id))
	chr_str_for = chrs[str(chr_id)][0]
	chr_str_rev = chrs[str(chr_id)][1]
	for a in range(0,len(genes)):
		if genes[a].seqid == chr_id:
			f.write("%s/%s\n\n" %(genes[a].id, chr_id))

#If genes are located on the positive strand, the codes takes the first 300 bp of each gene:
			if genes[a].strand == "+":
				start = genes[a].start
				end = genes[a].start + (genes[a].stop - genes[a].start)//3

				gRNA_gene = chr_str_for[start:end]
#Finds PAM sites (NGG) and gRNAs for the first 300 bp of each gene located at the positive strand:
				for i in range(22, len(gRNA_gene)):
					if gRNA_gene[i]=="G" and gRNA_gene[i-1]=="G":

						PAM_FWD = gRNA_gene[i-2 : i+1]
						gRNA_FWD = gRNA_gene[i-22 : i-2]
						gRNA_FWD_Core = gRNA_gene[i-14 : i+1]

						G = 0
						C = 0
						for n in gRNA_FWD:
							if n == "G":
								G += 1
							if n == "C":
								C += 1
							else:
								continue
						GC_Per_FWD = ("%f" %(((G + C)/20)*100))

#Applying some filters on found gRNAs before writing them out on the output file:
						if len(re.findall(r'TTT[T]+', gRNA_FWD)) >= 1:
							break
						else:
							for chr_id in chrs.keys():
							 canwrite = True
							 if genes[a].seqid == chr_id:
							  if not ( len(re.findall(gRNA_FWD_Core, chrs[str(chr_id)][0])) == 1 or len(re.findall(gRNA_FWD_Core, chrs[str(chr_id)][1])) == 0):
							   canwrite = False
							 if genes[a].seqid != chr_id:
							  if not ( len(re.findall(gRNA_FWD_Core, chrs[str(chr_id)][0])) == 0 or len(re.findall(gRNA_FWD_Core, chrs[str(chr_id)][1])) == 0):
							   canwrite = False

						if canwrite  == True:
							f.write("%s\t%s\t%s\n" %(PAM_FWD, gRNA_FWD, GC_Per_FWD))
#Finds PAM and gRNA in the negative strand:
				reversed_gRNA_gene = (''.join(reversed(gRNA_gene)))
				gRNA_gene_rev = []
				for i in range(0,len(reversed_gRNA_gene)):
					gRNA_gene_rev.append((reversed_gRNA_gene[i].translate(str.maketrans("ATCG", "TAGC"))))
				gRNA_gene_r = (''.join(gRNA_gene_rev))
				f.write("\n")

				for j in range(22, len(gRNA_gene_r)):
					if gRNA_gene_r[j]=="G" and gRNA_gene_r[j-1]=="G":

						PAM_RVS = gRNA_gene_r[j-2 : j+1]
						gRNA_RVS = gRNA_gene_r[j-22 : j-2]
						gRNA_RVS_Core = gRNA_gene_r[j-14 : j+1]

						C = 0
						G = 0

						for n in gRNA_RVS:
							if n == "G":
								G += 1
							if n == "C":
								C += 1
							else:
								continue
						GC_Per_RVS = ("%f" %(((G + C)/20)*100))

						if len(re.findall(r'TTT[T]+', gRNA_RVS)) >= 1:
							break
						else:
							for chr_id in chrs.keys():
							 canwrite = True
							 if genes[a].seqid == chr_id:
							  if not ( len(re.findall(gRNA_RVS_Core, chrs[str(chr_id)][0])) == 0 or len(re.findall(gRNA_RVS_Core, chrs[str(chr_id)][1])) == 1):
							   canwrite = False
							 if genes[a].seqid != chr_id:
							  if not ( len(re.findall(gRNA_RVS_Core, chrs[str(chr_id)][0])) == 0 or len(re.findall(gRNA_RVS_Core, chrs[str(chr_id)][1])) == 0):
							   canwrite = False
						if canwrite  == True:
							f.write("%s\t%s\t%s\n" %(PAM_RVS, gRNA_RVS, GC_Per_RVS))
				f.write("\n")
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			if genes[a].strand == "-":
				start = genes[a].stop
				end = start - (genes[a].stop - genes[a].start) // 3

				gRNA_gene = chr_str_for[end:start]

				for i in range(22, len(gRNA_gene)):
					if gRNA_gene[i] == "G" and gRNA_gene[i-1] == "G":

						PAM_FWD = gRNA_gene[i-2 : i+1]
						gRNA_FWD = gRNA_gene[i-22 : i-2]
						gRNA_FWD_Core = gRNA_gene[i-14 : i+1]

						G = 0
						C = 0
						for n in gRNA_FWD:
							if n == "G":
								G += 1
							if n == "C":
								C += 1
							else:
								continue
						GC_Per_FWD = ("%f" %(((G + C)/20)*100))

						if len(re.findall(r'TTT[T]+', gRNA_FWD)) >= 1:
							break
						else:
							for chr_id in chrs.keys():
							 canwrite = True
							 if genes[a].seqid == chr_id:
							  if not ( len(re.findall(gRNA_FWD_Core, chrs[str(chr_id)][0])) == 1 or len(re.findall(gRNA_FWD_Core, chrs[str(chr_id)][1])) == 0):
							   canwrite = False
							 if genes[a].seqid != chr_id:
							  if not ( len(re.findall(gRNA_FWD_Core, chrs[str(chr_id)][0])) == 0 or len(re.findall(gRNA_FWD_Core, chrs[str(chr_id)][1])) == 0):
							   canwrite = False
						if canwrite  == True:
							f.write("%s\t%s\t%s\n" %(PAM_FWD, gRNA_FWD, GC_Per_FWD))

#Finds the PAM and gRNAs in the negative strand:
				reversed_gRNA_gene = (''.join(reversed(gRNA_gene)))
				gRNA_gene_rev_str = []
				for i in range(0, len(reversed_gRNA_gene)):
					gRNA_gene_rev_str.append(reversed_gRNA_gene[i].translate(str.maketrans("ATCG", "TAGC")))
				gRNA_gene_r = (''.join(gRNA_gene_rev_str))
				f.write("\n")

				for j in range(22, len(gRNA_gene_r)):
					if gRNA_gene_r[j]=="G" and gRNA_gene_r[j-1]=="G":

						PAM_RVS = gRNA_gene_r[j-2 : j+1]
						gRNA_RVS = gRNA_gene_r[j-22 : j-2]
						gRNA_RVS_Core = gRNA_gene_r[j-14 : j+1]

						C = 0
						G = 0
						for n in gRNA_RVS:
							if n == "G":
								G += 1
							if n =="C":
								C += 1
							else:
								continue
						GC_Per_RVS = ("%f" %(((G + C)/20)*100))

						if len(re.findall(r'TTT[T]+', gRNA_RVS)) >= 1:
							break
						else:
							for chr_id in chrs.keys():
							 canwrite = True
							 if genes[a].seqid == chr_id:
							  if not ( len(re.findall(gRNA_RVS_Core, chrs[str(chr_id)][0])) == 1 or len(re.findall(gRNA_RVS_Core, chrs[str(chr_id)][1])) == 0):
							   canwrite = False
							 if genes[a].seqid != chr_id:
							  if not ( len(re.findall(gRNA_RVS_Core, chrs[str(chr_id)][0])) == 0 or len(re.findall(gRNA_RVS_Core, chrs[str(chr_id)][1])) == 0):
							   canwrite = False
						if canwrite  == True:
							f.write("%s\t%s\t%s\n" %(PAM_RVS, gRNA_RVS, GC_Per_RVS))

				f.write("\n")
