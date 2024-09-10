
## Description
Genetic variants could affect gene expression through various mechanisms depending on their genomic location. We found that genetic variants near splice sites are associated with abnormal intronic polyadenylation (IPA). Here we give examples to show how to detect variant-associated IPA events and variant-associated intron retention (IR) events.

#### RNAseq density plots and IGV browser screenshots showing that somatic single nucleotide variant (SNV) could cause intronic polyadenylation in *GRHL2* (A) and intron retention in *TP53* (B).

<img src="https://github.com/ZhaozzReal/SNV_IPA/blob/main/Examples.png" />



## Detect SNV-associated IPA event in *GRHL2*

```python Detect_SNV_IPA.py -bam GRHL2_TCGA-95-7043_test.bam -anno_txt hg38_annotation.txt -snv_file GRHL2_SNV_information.txt -output GRHL2_SNV_IPA.txt```

**GRHL2_TCGA-95-7043_test.bam:** a bamfile extracted reads mapping the gene *GRHL2* from RNAseq bamfile of patient TCGA-95-7043

**hg38_annotation.txt:**  a file contains all information of introns and exons annotated by Refseq

**GRHL2_SNV_information.txt:** a file contains information about somatic variants of *GRHL2* in patient TCGA-95-7043 



## Detect SNV-associated IR event in *TP53*

```python Detect_SNV_IR.py -bam TP53_TCGA-LA-A446_test.bam -anno_txt hg38_annotation.txt -snv_file TP53_SNV_information.txt -output TP53_SNV_IR.txt```

**TP53_TCGA-LA-A446_test.bam:**  a bamfile extracted reads mapping the gene *TP53* from RNAseq bamfile of patient TCGA-LA-A446

**GRHL2_SNV_information.txt:**  a file contains information about somatic variants of *GRHL2* in patient TCGA-95-7043 



**Statement:** The Cancer Genome Atlas (TCGA) characterizes a comprehensive list of genomic, epigenomic, and transcriptomic features in thousands of tumor samples. Here we give the examples to show how to perform variant-associated IPA analysis and variant-associated IR analysis. One can utilize our code functions in the python script file to further perform variant-associated IPA/IR analysis by intergrating more WGS/WES and RNAseq datasets.

**The following python packages are necessary:**

HTSeq、itertools、numpy、collections、multiprocessing、scipy、argparse、re、os



## Citation

*Please cite the following articles if you use IPAFinder in your research:*

* Zhao Z, Xu Q, Wei R, Wang W, Ding D, Yang Y, Yao J, Zhang L, Hu YQ, Wei G, Ni T. Cancer-associated dynamics and potential regulators of intronic polyadenylation revealed by IPAFinder using standard RNA-seq data. Genome Res. 2021 Sep 2. doi: 10.1101/gr.271627.120. PMID: 34475268.

* Zhao Z, Xu Q, Wei R, Huang L, Wang W, Wei G, Ni T. Comprehensive characterization of somatic variants associated with intronic polyadenylation in human cancers. Nucleic Acids Res. 2021 Oct 11. doi: 10.1093/nar/gkab772. PMID: 34508351.



## Contact

If you have any comments, suggestions, questions, bug reports, etc, feel free to contact Zhaozhao Zhao (zhaozzbio@126.com). 
And PLEASE attach your command line and log messages if possible.
