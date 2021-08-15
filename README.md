
## Description
Somatic single nucleotide variants (SNVs) in cancer genome affect gene expression through various mechanisms depending on their genomic location. In this study, we found that somatic SNVs near splice site are associated with abnormal intronic polyadenylation (IPA) . Here we give examples to show how to detect SNV-associated IPA events and SNV-associated intron retention (IR) events.

#### RNAseq density plots and IGV browser screenshots showing that somatic variant could cause intronic polyadenylation in *GRHL2* (A) and intron retention in *TP53* (B).

<img src="https://github.com/ZhaozzReal/SNV_IPA/blob/main/Examples.png" width="500" height="300"/>





### Detect SNV-associated IPA event in *GRHL2*

```python Detect_SNV_IPA.py -bam GRHL2_TCGA-95-7043_test.bam -anno_txt hg38_annotation.txt -snv_file GRHL2_SNV_information.txt -output GRHL2_SNV_IPA.txt```

**GRHL2_TCGA-95-7043_test.bam:** a bamfile extracted reads mapping the gene *GRHL2* from RNAseq bamfile of patient TCGA-95-7043

**hg38_annotation.txt:**  a file contains all information of introns and exons annotated by Refseq

**GRHL2_SNV_information.txt: ** a file contains information about somatic variants of *GRHL2* in patient TCGA-95-7043 





### Detect SNV-associated IR event in *TP53*

```python Detect_SNV_IR.py -bam TP53_TCGA-LA-A446_test.bam -anno_txt hg38_annotation.txt -snv_file TP53_SNV_information.txt -output TP53_SNV_IR.txt```

**TP53_TCGA-LA-A446_test.bam:**  a bamfile extracted reads mapping the gene *TP53* from RNAseq bamfile of patient TCGA-LA-A446

**GRHL2_SNV_information.txt:**  a file contains information about somatic variants of *GRHL2* in patient TCGA-95-7043 



**Statement:** The Cancer Genome Atlas (TCGA) characterizes a comprehensive list of genomic, epigenomic, and transcriptomic features in thousands of tumor samples. Here we give the examples to show how to perform SNV-associated IPA analysis and SNV-associated IR analysis. One can utilize our code functions in the python script file to further perform SNV-associated IPA analysis and SNV-associated IR analysis by intergrating more whole-exome sequencing and RNAseq datasets.

**The following python packages are necessary:**

HTSeq、itertools、numpy、collections、multiprocessing、scipy、argparse、re、os
