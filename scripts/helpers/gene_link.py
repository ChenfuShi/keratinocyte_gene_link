#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# Helper function to link genes to snps

#########################################

import pandas as pd
import numpy as np
import pybedtools as pbed
import os

os.makedirs("/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools/", exist_ok = True)
pbed.helpers.set_tempdir("/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools/")
bed_genome_file = "/mnt/iusers01/jw01/mdefscs4/hg38.genome"



# gene linker
# new version that is so much faster using pybedtools
# add genes that are overlapping the snps that overlap active 27ac
# promoters are 1 bp and snps are extendend 1kb, equivalent to saying promoters are within 1 kb from active snps

def _find_overlapping_genes(snps_df, peaks_df, trans_pbed, SNP_filter_peak_overlapping):
    snps_pbed = pbed.BedTool.from_dataframe(snps_df)
    peaks_pbed = pbed.BedTool.from_dataframe(peaks_df)
    
    if SNP_filter_peak_overlapping:
        active_snps = snps_pbed.intersect(peaks_pbed, wa = True)
    else:
        active_snps = snps_pbed
    active_snps = active_snps.slop(g=bed_genome_file, b = 1000)
    
    found_intersect = trans_pbed.intersect(active_snps, wa=True, wb=True)
    found_intersect_df = found_intersect.to_dataframe(header=None,disable_auto_names=True).iloc[:,[3,4,8]]
    found_intersect_df.columns = ["gene_id","transcript_id","linked_SNP"]
    found_intersect_df = found_intersect_df.groupby('gene_id').agg({'transcript_id':lambda x: set(x), "linked_SNP":lambda x: set(x)})
    found_intersect_df.loc[:,"score"] = 1
    return found_intersect_df
    
    
def _find_all_overlapping_genes(snp_df,dict_of_dfs_peaks,transcripts_bed, SNP_filter_peak_overlapping):
    all_genes = []
    for sample in dict_of_dfs_peaks.keys():
        all_genes.append(_find_overlapping_genes(snp_df,dict_of_dfs_peaks[sample],transcripts_bed, SNP_filter_peak_overlapping)) 
    # merging all samples into one table
    genes_unique_df = pd.DataFrame(columns="gene_id transcript_id".split())
    for x,name in zip(all_genes,dict_of_dfs_peaks.keys()):
        genes_unique_df = genes_unique_df.merge(x.rename(columns={"score":name, "transcript_id":"transcript_" + name, "linked_SNP":"linked_SNP_" + name}), on=["gene_id"], how="outer")
    return genes_unique_df

def _join_set(x):
    l = []
    for i in x:
        l.extend(eval(i))
    return set(l)

def _call_OE_genes(SNPs,loops,peaks,transcripts_bed, SNP_filter_peak_OE):
    # get the genes from other ends
    peaks_bed = pbed.BedTool.from_dataframe(peaks)
    loops_bed = pbed.BedTool.from_dataframe(loops)
    inverted_loops_bed = pbed.BedTool.from_dataframe(loops.iloc[:,[3,4,5,0,1,2,6,7,8,9,10,11,12,13,14]])
    snps_bed = pbed.BedTool.from_dataframe(SNPs)
    # filter SNPs with the peaks
    if SNP_filter_peak_OE:
        snps_bed_filtered = snps_bed.intersect(peaks_bed,u=True)
    else:
        snps_bed_filtered = snps_bed

    # get other ends
    OE_a = loops_bed.intersect(b=snps_bed_filtered.slop(b=5000,g=bed_genome_file),wa=True,wb=True).to_dataframe(header=None,disable_auto_names=True).iloc[:,[3,4,5,8,18]]
    OE_b = inverted_loops_bed.intersect(b=snps_bed_filtered.slop(b=5000,g=bed_genome_file),wa=True,wb=True)
    if os.path.getsize(OE_b.fn) == 0:
        OEs = OE_a
    else:
        OE_b = OE_b.to_dataframe(header=None,disable_auto_names=True).iloc[:,[3,4,5,8,18]]
        OEs = pd.concat((OE_a,OE_b))
    OEs = OEs.groupby([3,4,5]).agg({8:min,18:set}).reset_index()
    # intersect genes 
    OE_pbed = pbed.BedTool.from_dataframe(OEs)
    found_intersect = transcripts_bed.intersect(OE_pbed.slop(b=5000,g=bed_genome_file), wa=True, wb=True)
    found_intersect_df = found_intersect.to_dataframe(header=None,disable_auto_names=True).iloc[:,[3,4,8,9]]
    found_intersect_df.columns = ["gene_id","transcript_id","score","linked_SNP"]
    found_intersect_df = found_intersect_df.groupby('gene_id').agg({'score': 'min', 'transcript_id':lambda x: set(x), "linked_SNP":lambda x: _join_set(x)}) 
    return found_intersect_df

def link_genes(SNPs, dict_loops, dict_peaks, gtf_transcripts, SNP_filter_peak_OE = False, SNP_filter_peak_overlapping = True):
    # wrapper function to call the genes from the SNPs
    """ implemented as follows: 
    SNPs must overlap peaks in relevant cell type to be used.
    call promoters within 1kb of the snps as overlapping;
    call genes within 5kb of a loop linking active snps within 5kb of the other end of the loop"""
    
    # prepare TSS
    transcripts_tss = gtf_transcripts[["seqname","TSS_start","gene_id","transcript_id"]].copy()
    transcripts_tss["TSS_end"] = transcripts_tss["TSS_start"] + 1
    transcripts_tss.rename(columns={"seqname":"chr","TSS_start":"start","TSS_end":"end"}, inplace=True)
    transcripts_tss = transcripts_tss[transcripts_tss["chr"].str.startswith("chr")]
    transcripts_bed = pbed.BedTool.from_dataframe(transcripts_tss[["chr","start","end","gene_id","transcript_id"]])  
    # call genes for each sample
    all_OE_genes = []
    for sample in dict_loops.keys():
        all_OE_genes.append(_call_OE_genes(SNPs,dict_loops[sample],dict_peaks[sample],transcripts_bed, SNP_filter_peak_OE))
    OE_merged = pd.DataFrame(columns="gene_id transcript_id".split())
    for x,name in zip(all_OE_genes,dict_loops.keys()):
        OE_merged = OE_merged.merge(x.rename(columns={"score":name, "transcript_id":"transcript_" + name, "linked_SNP":"linked_SNP_" + name}), on=["gene_id"], how="outer")
    Overlapping = _find_all_overlapping_genes(SNPs,dict_peaks,transcripts_bed, SNP_filter_peak_overlapping)
    
    return OE_merged , Overlapping
    
