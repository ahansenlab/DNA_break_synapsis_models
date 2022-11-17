#!/usr/bin/env python
# coding: utf-8

# In[31]:


import subprocess
import os
import sys
#ResultFolderName = "PsiblastRepresentatives_aligned"
#with open("DGR_result_cluster_checked_combined", "r") as f:
#    lines = f.readlines()

output_folder = sys.argv[1]
FileName_list = [filename for filename in os.listdir(output_folder) if filename.startswith("LEFcount") and filename.endswith("bed")]

for FileName in FileName_list:
    
    BedgraphFile = FileName.replace("bed","bedgraph")

    p = subprocess.Popen("""awk '{printf  "%s\\t%d\\t%d\\t%2.3f\\n" , $1,$2,$3,$5}' """ + output_folder+"/"+FileName + " > "+ output_folder+"/"+BedgraphFile, shell = True)
    p.communicate(timeout=9999999)
    p.wait(timeout=9999999)
    
    SortedBedgraphFile = FileName.replace(".bed","_sorted.bedgraph")
    
    p = subprocess.Popen("sort -k1,1 -k2,2n " + output_folder+"/"+BedgraphFile + " > "+ output_folder+"/"+SortedBedgraphFile, shell = True)
    p.communicate(timeout=9999999)
    p.wait(timeout=9999999)
    
    BigwigFile = FileName.replace("bed","bw")
    p = subprocess.Popen("bedGraphToBigWig " + output_folder+"/"+SortedBedgraphFile + " " + output_folder+"/DSBchrom.sizes "  + output_folder+"/"+BigwigFile, shell = True)
    p.communicate(timeout=9999999)
    p.wait(timeout=9999999)

# In[ ]:




