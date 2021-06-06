
## Description
Somatic single nucleotide variants (SNVs) in cancer genome affect gene expression through various mechanisms depending on their genomic location. In this study, we found that somatic SNVs near splice site are associated with abnormal intronic polyadenylation (IpA). Here we show how to perform SNV-associated IpA analysis through an example (*GRHL2*).


#### RNAseq density plots and IGV browser screenshot showing that somatic variant could cause abnormal intronic polyadenylation in *GRHL2*.
<img src="https://github.com/ZhaozzReal/SNV_IPA/blob/main/Example_GRHL2_SNV_IPA.png" width="400" height="400"/>



**Command Example:**

```python Detect_SNV_IPA.py -bam GRHL2_TCGA-95-7043_test.bam -anno_txt hg38_annotation.txt -snv_file GRHL2_SNV_information.txt -output GRHL2_SNV_IPA.txt```
##### GRHL2_TCGA-95-7043_test.bam: a bamfile extracted reads mapping the locus *GRHL2* from RNAseq bamfile of patient TCGA-95-7043
##### hg38_annotation.txt: a file contains all information of introns and exons annotated by Refseq
##### GRHL2_SNV_information.txt: a file contains information about somatic variants of GRHL2 in patient TCGA-95-7043 


**The following python packages are necessary:**

HTSeq、itertools、numpy、collections、multiprocessing、scipy、argparse、re、os
