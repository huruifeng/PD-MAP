Last login: Fri Feb 17 16:41:20 on ttys006

The default interactive shell is now zsh.
To update your account to use zsh, please run `chsh -s /bin/zsh`.
For more details, please visit https://support.apple.com/kb/HT208050.
(base) rh460@PHS025804:~$ e2
Last login: Fri Feb 17 16:41:28 2023 from 10.22.205.238
--------------------------------------------------------------------------------
|         Welcome to ERISTwo, the HPC Supercomputer of the Scientific          |
|          Computing ecosystem managed by Enterprise Research IS (ERIS).       |
| USERS MUST COMPLY WITH MASS GENERAL BRIGHAM POLICY: All ACTIONS ARE RECORDED |
|------------------------------------------------------------------------------|
|    For information on how to use the cluster, see our knowledge base at:     |
|      https://rc.partners.org/kb/computational-resources/linux-cluster        |
|                                                                              |
|  Remember to always:                                                         |
|    1) Run jobs through the LSF job scheduler and reserve resources           |
|    2) Not run computational jobs on the login nodes (eris2n4 and eris2n5)    |
|    3) Keep a backup of essential data outside ERISTwo                        |
|    4) Remove patient names and other PII / PHI from filenames                |
|                                                                              |
|  Email questions, problems and software requests to hpcsupport@partners.org  |
--------------------------------------------------------------------------------

(base) [rh460@eris2n5 ~]$ bioinfo
(base) [rh460@eris2n5 bioinformatics]$ ll
total 59472
drwxrwx---.   7 root        4096 Jan 27 17:50 .
drwxr-xr-x. 308 root           0 Feb 13 18:27 ..
-r--r--r--.   1 5451112     4096 Aug  6  2022 ._bioinformatics_filestat_table.csv
-rw-r--r--.   1 root       38629 Jan 16  2021 bioinformatics_filestat_table.csv
-rw-r--r--.   1 4975670     4096 May 24  2021 ._.DS_Store
-rw-r--r--.   1 4975670    12292 Jul 31  2021 .DS_Store
drwxrws---.  11 xd010       4096 Jul 25  2022 external_data
-rw-r--r--.   1 root       16059 Jan 27 17:50 fsstats_report.csv
-rw-r--r--.   1 root       28980 Jan 27 17:50 fsstats_report.txt
drwxrws---.   4 xd010       4096 Jul  2  2022 projects
drwxrws---.   8 xd010       4096 Mar 14  2022 referenceGenome
drwxr-xr-x.   4 4975670     4096 Nov 26  2020 .Rproj.user
drwxrwx---.   4 tz949       4096 Dec 29 17:03 tools
-rw-r--r--.   1 root    23077537 Jan 22 23:54 usage_report_old.txt
-rw-r--r--.   1 root    24267863 Jan 22 23:54 usage_report.txt
(base) [rh460@eris2n5 bioinformatics]$ cd projects/
(base) [rh460@eris2n5 projects]$ ll
total 352
drwxrws---.  4 xd010 4096 Jul  2  2022 .
drwxrwx---.  7 root  4096 Jan 27 17:50 ..
drwxrws---. 15 xd010 4096 Dec 30 11:01 biohub
drwxrws---. 11 xd010 4096 Nov 26 11:48 donglab
-rw-r-Sr--.  1 xd010 4096 Feb  9  2021 ._.DS_Store
-rw-r-Sr--.  1 xd010 8196 Oct 23 22:37 .DS_Store
(base) [rh460@eris2n5 projects]$ cd donglab/
(base) [rh460@eris2n5 donglab]$ ll
total 992
drwxrws---. 11 xd010    4096 Nov 26 11:48 .
drwxrws---.  4 xd010    4096 Jul  2  2022 ..
drwxrws---.  3 xd010   20480 Jul 27  2022 acne_gwas_by_sex
drwxrws---.  4 rw552    4096 Dec 21 16:05 AMPPD_circRNA
drwxrws---.  8 xd010    4096 Dec 15 17:14 AMPPD_eRNA
drwxrws---.  7 xd010   40960 Apr 11  2021 ariela2021
-rw-rw----.  1 xd010    4096 Aug  8  2022 ._.DS_Store
-rw-rw----.  1 xd010   10244 Feb 19 19:58 .DS_Store
drwxrws---.  2 xd010    4096 Jul  2  2022 eccDNA
drwxrws---. 11 xd010    4096 Oct 21 22:06 fei2021
drwxrws---. 10 xd010    4096 Feb  1 14:56 intronAnalysis
drwxrwx---.  8 5401991  8192 Jun 17  2022 jhu2022
drwxrwx---.  4 rh460    4096 Mar 19  2022 mouse_receptor_2022
(base) [rh460@eris2n5 donglab]$ cd AMPPD_eRNA/
(base) [rh460@eris2n5 AMPPD_eRNA]$ ll
total 560
drwxrws---.  8 xd010 4096 Dec 15 17:14 .
drwxrws---. 11 xd010 4096 Nov 26 11:48 ..
drwxrwx---.  2 rw552 4096 Sep  2 17:30 bin
-rwxrwx---.  1 rw552  643 Aug 19  2022 config.txt
drwxr-sr-x.  4 rh460 4096 Nov 26 10:54 DE_eRNA
drwxr-sr-x.  3 rh460 4096 Nov 26 10:37 DE_genes
drwxrwx---.  7 rw552 4096 Feb 16 14:25 inputs
drwxrwx---. 11 rw552 4096 Jan 26 15:56 output
drwxrwx---.  2 rw552 8192 Feb  8 11:31 src
(base) [rh460@eris2n5 AMPPD_eRNA]$ cd src/
(base) [rh460@eris2n5 src]$ ll
total 272
drwxrwx---. 2 rw552  8192 Feb  8 11:31 .
drwxrws---. 8 xd010  4096 Dec 15 17:14 ..
-rwxrwx---. 1 rw552  1061 Sep 21 22:58 eRNA.characterize.merge.R
-rwxrwx---. 1 rw552 22301 Dec 14 16:25 eRNA.characterize.sh
(base) [rh460@eris2n5 src]$ cd ../inputs/
(base) [rh460@eris2n5 inputs]$ ll
total 11169600
drwxrwx---. 7 rw552       4096 Feb 16 14:25 .
drwxrws---. 8 xd010       4096 Dec 15 17:14 ..
drwxrwx---. 3 rw552       4096 Feb  8 16:46 bidirectional_pairs
-rwxrwx---. 1 rw552        375 Aug 16  2022 ChromInfo.txt
drwxrwx---. 2 rw552       4096 Dec  2 09:41 class
drwxrwx---. 2 rw552       4096 Dec 14 16:37 closest
-rw-r-----. 1 rw552   17191688 Jan  9 15:46 eRNA.characterize.feature.color.xls
-rwxrwx---. 1 rw552 2764730108 Nov 21 12:32 eRNA.merged.readCounts.v2.xls
-rwxrwx---. 1 rw552    5671288 Dec 12 09:31 eRNA_stranded_sorted.bed
-rwxrwx---. 1 rw552    1504244 Dec  6 16:46 gencode.genes.no_version.txt
-rw-rw----. 1 rh460 3354189895 Feb 16 13:41 gene_expr_matrix_tpm_row_genes.txt
-rw-rw----. 1 rh460 2995596732 Feb 16 13:41 gene_expr_matrix_tpm_row_samples.txt
-rw-rw----. 1 rw552      66773 Feb 16 14:25 LITAF_gene_exp.txt
drwxrwx---. 2 rw552       4096 Feb 16 14:51 minus
drwxrwx---. 2 rw552       4096 Feb 17 12:00 plus
-rwxrwx---. 1 rw552        374 Dec  6 16:46 README.txt
-rwxrwx---. 1 rw552    4637543 Aug 16  2022 toExclude.bed
(base) [rh460@eris2n5 inputs]$ cd ../
(base) [rh460@eris2n5 AMPPD_eRNA]$ ll
total 560
drwxrws---.  8 xd010 4096 Dec 15 17:14 .
drwxrws---. 11 xd010 4096 Nov 26 11:48 ..
drwxrwx---.  2 rw552 4096 Sep  2 17:30 bin
-rwxrwx---.  1 rw552  643 Aug 19  2022 config.txt
drwxr-sr-x.  4 rh460 4096 Nov 26 10:54 DE_eRNA
drwxr-sr-x.  3 rh460 4096 Nov 26 10:37 DE_genes
drwxrwx---.  7 rw552 4096 Feb 16 14:25 inputs
drwxrwx---. 11 rw552 4096 Jan 26 15:56 output
drwxrwx---.  2 rw552 8192 Feb  8 11:31 src
(base) [rh460@eris2n5 AMPPD_eRNA]$ cd bin/
(base) [rh460@eris2n5 bin]$ ll
total 176
drwxrwx---. 2 rw552 4096 Sep  2 17:30 .
drwxrws---. 8 xd010 4096 Dec 15 17:14 ..
-rwxrwx---. 1 rw552  889 Sep  2 17:30 getNormalizedCpGscore.awk
(base) [rh460@eris2n5 bin]$ cd ../
(base) [rh460@eris2n5 AMPPD_eRNA]$ ll
total 560
drwxrws---.  8 xd010 4096 Dec 15 17:14 .
drwxrws---. 11 xd010 4096 Nov 26 11:48 ..
drwxrwx---.  2 rw552 4096 Sep  2 17:30 bin
-rwxrwx---.  1 rw552  643 Aug 19  2022 config.txt
drwxr-sr-x.  4 rh460 4096 Nov 26 10:54 DE_eRNA
drwxr-sr-x.  3 rh460 4096 Nov 26 10:37 DE_genes
drwxrwx---.  7 rw552 4096 Feb 16 14:25 inputs
drwxrwx---. 11 rw552 4096 Jan 26 15:56 output
drwxrwx---.  2 rw552 8192 Feb  8 11:31 src
(base) [rh460@eris2n5 AMPPD_eRNA]$ cd src/
(base) [rh460@eris2n5 src]$ ll
total 272
drwxrwx---. 2 rw552  8192 Feb  8 11:31 .
drwxrws---. 8 xd010  4096 Dec 15 17:14 ..
-rwxrwx---. 1 rw552  1061 Sep 21 22:58 eRNA.characterize.merge.R
-rwxrwx---. 1 rw552 22301 Dec 14 16:25 eRNA.characterize.sh
(base) [rh460@eris2n5 src]$ vi eRNA.characterize.merge.R
(base) [rh460@eris2n5 src]$ vi eRNA.characterize.sh

# ============================
# merging data files for both strands
# ============================
### MOST LIKELY WILL BE DEPRECIATED AFTER WE FIND A BETTER METHOD
########## GWAS
cat $pipeline_path/output/minus/eRNA.enrichment/eRNA.minus.f16.GWASDisease.txt $pipeline_path/output/plus/eRNA.enrichment/eRNA.plus.f16.GWASDisease.txt | sort -k 9 | cut -f 5-9 | uniq | bedtools groupby -g 5 -c 5 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > $pipeline_path/output/merged/eRNA.f16.GWASDisease.counts
# just for testing
#cat $pipeline_path/output/minus/eRNA.enrichment/eRNA.minus.f16.GWASDisease.txt $pipeline_path/output/plus/eRNA.enrichment/eRNA.plus.f16.GWASDisease.txt | sort -k 9 | cut -f 5-9 | bedtools groupby -g 5 -c 4 -o concat |less

########## eQTL
cat $pipeline_path/output/minus/eRNA.enrichment/eRNA.minus.f18.eSNP.gtexCaviarDisease.txt $pipeline_path/output/plus/eRNA.enrichment/eRNA.plus.f18.eSNP.gtexCaviarDisease.txt | sort -k 9 | cut -f 5-9 | uniq | bedtools groupby -g 5 -c 5 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > $pipeline_path/output/merged/eRNA.f18.eSNP.gtexCaviarDisease.txt
cat $pipeline_path/output/minus/eRNA.enrichment/eRNA.minus.f18.eSNP.gtexDapgDisease.txt $pipeline_path/output/plus/eRNA.enrichment/eRNA.plus.f18.eSNP.gtexDapgDisease.txt | sort -k 9 | cut -f 5-9 | uniq | bedtools groupby -g 5 -c 5 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > $pipeline_path/output/merged/eRNA.f18.eSNP.gtexDapgDisease.txt
cat $pipeline_path/output/minus/eRNA.enrichment/eRNA.minus.f18.eSNP.pvalDisease.txt $pipeline_path/output/plus/eRNA.enrichment/eRNA.plus.f18.eSNP.pvalDisease.txt | sort -k 9 | cut -f 5-9 | uniq | bedtools groupby -g 5 -c 5 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > $pipeline_path/output/merged/eRNA.f18.eSNP.pvalDisease.txt

########## CAGE
cat $pipeline_path/output/minus/eRNA.minus.f08.CAGEenhtissue.txt $pipeline_path/output/plus/eRNA.plus.f08.CAGEenhtissue.txt | sort -k 10 | cut -f 5-9 | uniq | bedtools groupby -g 5 -c 5 -o count | sort -k 1 | bedtools groupby -g 1 -c 2 -o sum > $pipeline_path/output/merged/eRNA.f18.eSNP.gtexCaviarDisease.txt

### TODO
############################################################################################################################################################

# ====================================
# bidirectional_trans
# ====================================
#select only the columns you need (name, RPM, pair, RPM, dis)
#cut -f4,7,14,17,21 plus.minus.bed > cut.plus.minus.bed
#cut -f4,7,14,17,21 minus.plus.bed > cut.minus.plus.bed

#TODO test this
if [[ $STRAND == "minus" ]];
then
Rscript eRNA.bidir.R $EXTERNAL_FEATURE/bidir/cut.plus.minus.bed plus.minus
elif [[ $STRAND == "plus" ]];
then
Rscript eRNA.bidir.R $EXTERNAL_FEATURE/bidir/cut.minus.plus.bed minus.plus
fi
############################################################################################################################################################
