import HTSeq
from multiprocessing import Pool
import collections
import itertools
import numpy as np
import os
import re
import scipy as sp
import scipy.stats
from argparse import ArgumentParser,ArgumentTypeError
parser = ArgumentParser(description="Analysis of SNV-associated intronic polyadenylation events,designed by Zhaozz")
parser.add_argument("-bam",dest='bamfile',action="store",type=str,help="input bamfile")
parser.add_argument('-anno_txt',dest='anno_txt',action="store",type=str,help="input annotation file contains intron and flanking exons")
parser.add_argument("-snv_file",dest="snv_file",action="store",type=str,help="input file contains SNV information")
parser.add_argument("-output",dest="outfile",action="store",type=str,help="Output SNV-associated intronic polyadenylation information")
args = parser.parse_args()

cigar_char = ('M', '=', 'X')
minaqual = 10
def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2

def Get_label_information(label,annot,bam_reader):
    gas = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
    gene_count = {}
    for feature,rank,chrom,start,end,strand,length,exon_rank_left,exon_rank_right in annot[label]:
        iv=HTSeq.GenomicInterval(chrom,start,end,strand)
        gas[iv] += (feature,rank)
        gene_count[(feature,rank)] = 0
    boundary_left,boundary_right = min([i[3] for i in annot[label]]),max([i[4] for i in annot[label]])
    region_fetch = annot[label][0][2]+":"+str(int(boundary_left)-500)+"-"+str(int(boundary_right)+500)
    read_seq = bam_reader.fetch(region=region_fetch)
    read_seq_iter = iter(bam_reader.fetch())
    one_read = next(read_seq_iter)
    pe_mode = one_read.paired_end
    #pe_mode = True
    if pe_mode:
        read_seq = HTSeq.pair_SAM_alignments_with_buffer(read_seq)
    for a in read_seq:
        if not pe_mode:
            if not a.aligned:
                continue
            if a.optional_field('NH') > 1:
                continue
            iv_seq = (cigop.ref_iv for cigop in a.cigar if cigop.type == "M" and cigop.size >0)
        else:
            if ((a[0] and a[0].aQual<minaqual) or (a[1] and a[1].aQual<minaqual)):
                continue
            if ((a[0] and a[0].optional_field('NH') > 1) or (a[1] and a[1].optional_field('NH')>1)):
                continue
            if a[0] is not None and a[0].aligned:
                iv_seq = (cigop.ref_iv for cigop in a[0].cigar if cigop.type in cigar_char and cigop.size > 0)
            else:
                iv_seq=tuple()
            if a[1] is not None and a[1].aligned:
                iv_seq = itertools.chain(iv_seq,(invert_strand(cigop.ref_iv) for cigop in a[1].cigar if cigop.type in cigar_char and cigop.size > 0))
        feature_aligned=set()
        for iv in iv_seq:
            for iv2, val2 in gas[iv].steps():
                feature_aligned |= val2
                ga[iv] += 1  # for calculating coverage
        if len(feature_aligned) ==0:
            continue
        for f in [item for item in feature_aligned if item[0] == 'intron']:
            gene_count[f] +=1
        if 'intron' not in [x for x,y in feature_aligned]:
            for f in feature_aligned:
                gene_count[f] +=1
    return gas,ga,gene_count


def Estimation_abundance(Region_Coverage,break_point):
    downstream_cov_mean=np.mean(Region_Coverage[break_point:])
    upstream_cov_mean=np.mean(Region_Coverage[:break_point])
    Coverage_diff=Region_Coverage[:break_point]-upstream_cov_mean
    Coverage_diff=np.append(Coverage_diff,Region_Coverage[break_point:]-downstream_cov_mean)
    Mean_Squared_error=np.mean(Coverage_diff**2)
    return Mean_Squared_error


def Get_min_mseratio(cvg_region):
    search_result=[[],[]]
    global_mse=np.mean((cvg_region-np.mean(cvg_region))**2)
    cvg_result=[[],[]]
    for curr_search_point in range(50,len(cvg_region),int(len(cvg_region)/30)):
        Mean_Squared_error = Estimation_abundance(cvg_region, curr_search_point)
        mse_ratio=Mean_Squared_error/global_mse
        cvg_result[0].append(mse_ratio)
        cvg_result[1].append(curr_search_point)
    cvg_min_mse=min(cvg_result[0])
    cvg_min_mse_index=cvg_result[0].index(cvg_min_mse)
    cvg_min_mse_point=cvg_result[1][cvg_min_mse_index]
    return cvg_min_mse,cvg_min_mse_point


def Get_min_mseratio_refine(cvg_region,min_mse_point):
    research_result=[[],[]]
    global_mse=np.mean((cvg_region-np.mean(cvg_region))**2)
    for search_point in range(min_mse_point-int(len(cvg_region)/30),min_mse_point+int(len(cvg_region)/30),30):
        Mean_Squared_error = Estimation_abundance(cvg_region, search_point)
        mse_ratio=Mean_Squared_error/global_mse
        research_result[0].append(mse_ratio)
        research_result[1].append(search_point)
    min_mse_refine=min(research_result[0])
    min_mse_index=research_result[0].index(min_mse_refine)
    min_mse_point_refine=research_result[1][min_mse_index]
    return min_mse_refine,min_mse_point_refine


def Get_Skipend_dict(region_fetch,bamfile,strand):
    bam_reader=HTSeq.BAM_Reader(bamfile)
    read_seq = bam_reader.fetch(region=region_fetch)
    read_seq_iter = iter(bam_reader.fetch())
    one_read = next(read_seq_iter)
    skip_list=[]
    pe_mode = one_read.paired_end
    if pe_mode:
        read_seq = HTSeq.pair_SAM_alignments_with_buffer(read_seq)
    for a in read_seq:
        if not pe_mode:
            if not a.aligned:
                continue
            if a.optional_field('NH') > 1:
                continue
            if strand=="+":
                skip_list.extend([int(cigop.ref_iv.end) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
            else:
                skip_list.extend([int(cigop.ref_iv.start) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
        else:
            if ((a[0] and a[0].aQual<minaqual) or (a[1] and a[1].aQual<minaqual)):
                continue
            if ((a[0] and a[0].optional_field('NH') > 1) or (a[1] and a[1].optional_field('NH')>1)):
                continue
            if a[0] is not None and a[0].aligned:
                if strand=="+":
                    skip_list.extend([int(cigop.ref_iv.end) for cigop in a[0].cigar if cigop.type =="N" and cigop.size > 0])
                else:
                    skip_list.extend([int(cigop.ref_iv.start) for cigop in a[0].cigar if cigop.type =="N" and cigop.size > 0])
            if a[1] is not None and a[1].aligned:
                if strand=="+":
                    skip_list.extend([int(cigop.ref_iv.end) for cigop in a[1].cigar if cigop.type =="N" and cigop.size > 0])
                else:
                    skip_list.extend([int(cigop.ref_iv.start) for cigop in a[1].cigar if cigop.type =="N" and cigop.size > 0])
    skip_dict=dict(collections.Counter(skip_list))
    return skip_dict

def Get_Skipstart_dict(region_fetch,bamfile,strand):
    skip_list=[]
    bam_reader=HTSeq.BAM_Reader(bamfile)
    read_seq = bam_reader.fetch(region=region_fetch)
    read_seq_iter = iter(bam_reader.fetch())
    one_read = next(read_seq_iter)
    pe_mode = one_read.paired_end
    if pe_mode:
        read_seq = HTSeq.pair_SAM_alignments_with_buffer(read_seq)
    for a in read_seq:
        if not pe_mode:
            if not a.aligned:
                continue
            if a.optional_field('NH') > 1:
                continue
            if strand=="+":
                skip_list.extend([int(cigop.ref_iv.start) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
            else:
                skip_list.extend([int(cigop.ref_iv.end) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
        else:
            if ((a[0] and a[0].aQual<minaqual) or (a[1] and a[1].aQual<minaqual)):
                continue
            if ((a[0] and a[0].optional_field('NH') > 1) or (a[1] and a[1].optional_field('NH')>1)):
                continue
            if a[0] is not None and a[0].aligned:
                if strand=="+":
                    skip_list.extend([int(cigop.ref_iv.start) for cigop in a[0].cigar if cigop.type =="N" and cigop.size > 0])
                else:
                    skip_list.extend([int(cigop.ref_iv.end) for cigop in a[0].cigar if cigop.type =="N" and cigop.size > 0])
            if a[1] is not None and a[1].aligned:
                if strand=="+":
                    skip_list.extend([int(cigop.ref_iv.start) for cigop in a[1].cigar if cigop.type =="N" and cigop.size > 0])
                else:
                    skip_list.extend([int(cigop.ref_iv.end) for cigop in a[1].cigar if cigop.type =="N" and cigop.size > 0])
    skip_dict=dict(collections.Counter(skip_list))
    return skip_dict


def Get_Exonic_SNV_IPA(input_tuple):
    lastbase,bamfile = input_tuple
    label,region,SNV_position,strand,exon_rank,distance,intron_rank,mut_SS,patient_id,mut_change,mut_type = lastbase.strip().split("\t")
    exon_start,exon_end = [feature[3:5] for feature in annot[label] if feature[0]=="exon" and feature[1]== int(exon_rank.split("_")[1])][0]
    result = []
    exon_start = exon_start + 1
    if exon_start <= int(SNV_position) <= exon_end:
        bam_reader = HTSeq.BAM_Reader(bamfile)
        genome_base = mut_change.split(">")[0]
        chrom = label.split(":")[0]
        region_fetch = chrom + ":" + str(int(SNV_position)-500) + "-" + str(int(SNV_position)+500)
        read_seq = bam_reader.fetch(region=region_fetch)
        distance = min(abs(int(SNV_position) - exon_start),abs(int(SNV_position) - exon_end)) + 1
        splice_MUT = 0
        splice_WT = 0
        IR_MUT = 0
        IR_WT = 0
        for a in read_seq:
            for cigop in a.cigar:
                if strand == "+":
                    if cigop.type == "N" and ((cigop.ref_iv.start == exon_end and mut_SS == "5SS") or (cigop.ref_iv.end == exon_start-1 and mut_SS == "3SS")) :
                        if genome_base in a.optional_field('MD'):
                            splice_MUT += 1
                        else:
                            splice_WT += 1
                    elif (cigop.type == "M" and mut_SS == "5SS" and cigop.ref_iv.start < exon_end+3 < cigop.ref_iv.end and cigop.ref_iv.start < exon_end-3 < cigop.ref_iv.end) or (cigop.type == "M" and mut_SS =="3SS" and cigop.ref_iv.start < exon_start-3 < cigop.ref_iv.end and cigop.ref_iv.start < exon_start+3 < cigop.ref_iv.end):
                        if genome_base in a.optional_field('MD'):
                            IR_MUT +=1
                        else:
                            IR_WT +=1
                else:
                    if cigop.type == "N" and ((cigop.ref_iv.end == exon_start-1 and mut_SS == "5SS") or (cigop.ref_iv.start == exon_end and mut_SS == "3SS")) :
                        if genome_base in a.optional_field('MD'):
                            splice_MUT += 1
                        else:
                            splice_WT += 1
                    elif (cigop.type == "M" and mut_SS == "5SS" and cigop.ref_iv.start < exon_start+3 < cigop.ref_iv.end and cigop.ref_iv.start < exon_start-3 < cigop.ref_iv.end) or (cigop.type == "M" and mut_SS =="3SS" and cigop.ref_iv.start < exon_end-3 < cigop.ref_iv.end and cigop.ref_iv.start < exon_end+3 < cigop.ref_iv.end):
                        if genome_base in a.optional_field('MD'):
                            IR_MUT +=1
                        else:
                            IR_WT +=1
        ratio_val,Pvalue = sp.stats.fisher_exact([[splice_MUT,splice_WT],[IR_MUT,IR_WT]])
        gas,ga,gene_count = Get_label_information(label,annot,bam_reader)
        intronrank = int(intron_rank.split("_")[1])
        if intronrank != -1:
            if gene_count[('intron',intronrank)] > 30 :
                feature,rank,chrom,start,end,strand,length,exon_rank_left,exon_rank_right=[line for line in annot[label] if line[0]=="intron" and line[1] == intronrank][0]
                intron_start = start
                intron_end = end
                end_value = 10
                intron_region = chrom+":"+str(intron_start)+"-"+str(intron_end)
                IPAtype = "Composite"
                iv = HTSeq.GenomicInterval(chrom,start,end,strand)
                cvg_region = list(ga[iv])
                if strand == "-":
                    cvg_region = cvg_region[::-1]
                skipend_dict = Get_Skipend_dict(intron_region,bamfile,strand)
                for key,value in skipend_dict.items():
                    if int(start)+50<int(key)<int(end)-50 and int(value)>10:
                        if strand == "+":
                            skip_position = int(key)-int(start)
                        else:
                            skip_position = int(end)-int(key)
                        curr_label_all_cov = cvg_region[skip_position:]
                        IPAtype = "Skipped"
                        start = int(key)
                        end = int(key)
                        end_value = int(value)
                        break
                if len(cvg_region) > 250:
                    min_mse,min_mse_point = Get_min_mseratio(cvg_region)
                    if min_mse < 0.7 :
                        min_mse_refine,min_mse_point = Get_min_mseratio_refine(cvg_region,min_mse_point)
                        IPAabundance = np.mean(cvg_region[:min_mse_point])
                        up_down_diff = np.mean(cvg_region[:min_mse_point])-np.mean(cvg_region[min_mse_point:])
                        down_cov = len(list(filter(lambda x:x>5,cvg_region[min_mse_point:])))/(len(cvg_region)-min_mse_point)
                        up_cov = len(list(filter(lambda x:x>5,cvg_region[:min_mse_point])))/(len(cvg_region[:min_mse_point]))
                        IPA_point = min_mse_point
                        if min_mse_refine < 0.5 and up_down_diff > 2:
                            if strand == "+":
                                IPA_location = int(start)+IPA_point
                                IPA_inf = chrom+":"+str(start)+"-"+str(IPA_location)
                                exon_iv = tuple(i[0] for i in gas.steps() if i[1] == {('exon',int(exon_rank_left))})
                            else:
                                IPA_location = int(end)-IPA_point
                                IPA_inf = chrom+":"+str(IPA_location)+"-"+str(end)
                                exon_iv = tuple(i[0] for i in gas.steps() if i[1] == {('exon',int(exon_rank_right))})
                            skipstart_dict = Get_Skipstart_dict(intron_region,bamfile,strand)
                            for key,value in skipstart_dict.items():
                                if IPA_location-20< int(key) <IPA_location+20 and int(value) > 4:
                                    IPAtype = "A5SS"
                            if len(exon_iv) == 1:
                                exon_abundance = round(np.mean(sorted(list(ga[exon_iv[0]]),reverse=True)[:30]),3)
                                IPUI = round(IPAabundance/exon_abundance,3)
                                result.append([label,region,mut_SS,SNV_position,bamfile,patient_id,exon_rank,intron_rank,mut_change,mut_type,distance,splice_MUT,splice_WT,IR_MUT,IR_WT,Pvalue,intron_region,IPA_inf,IPAtype,min_mse_refine,up_cov,down_cov,IPUI])
    return result

refGene_txt = args.anno_txt
annot = collections.OrderedDict()
for line in open(refGene_txt):
    gene_label,feature,rank,position,length,exon_rank_left,exon_rank_right = line.strip().split('\t')
    chrom,iv_str,strand = position.strip().split(':')
    start,end = map(int,iv_str.strip().split('-'))
    annot.setdefault(gene_label,[]).append((feature,int(rank),chrom,start,end,strand,int(length),exon_rank_left,exon_rank_right))


SNV_list = []
for line in open(args.snv_file,"r"):
    SNV_list.append(line)

input_tuple = [SNV_list[0],args.bamfile]
result_list = Get_Exonic_SNV_IPA(input_tuple)

out = open(args.outfile,"w")
first_line = ["label","exon_region","mut_SS_type","SNV_position","bamfile","patient_id","exon_rank","intron_rank","mut_change","mut_type","distance","normal_splice_MUT","normal_splice_WT","abnormal_splice_MUT","abnormal_splice_WT","Pvalue","intron_region","IPA_inf","IPAtype","min_ratio_mse","up_cov","down_cov","IPUI"]
out.writelines("\t".join(first_line)+"\n")
for lst in result_list:
    out.writelines('\t'.join(list(map(str,lst)))+'\n')

out.close()
