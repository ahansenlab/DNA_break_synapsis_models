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
FileName_list = [filename for filename in os.listdir(output_folder) if filename.endswith("bw")]

for FileName in FileName_list:
    
    MatFile = FileName.replace("bw","mat.gz")
    p = subprocess.Popen("computeMatrix reference-point --referencePoint TSS -b 500 -a 500 -S " + output_folder+"/"+FileName + " -R " + output_folder+"/"+"DSBcoordinates.bed -o "+ output_folder+"/"+MatFile, shell = True)
    p.communicate(timeout=9999999)
    p.wait(timeout=9999999)
    
    pngFile = FileName.replace("bw","pdf")
    p = subprocess.Popen("plotHeatmap -m " + output_folder+"/"+MatFile + """ --colorMap 'Blues' --plotFileFormat pdf --yMin -1 --yMax 3 --zMin 0 --zMax 5 --missingDataColor 1 --refPointLabel "DSB" -out """ + output_folder+"/"+pngFile, shell = True)
    p.communicate(timeout=9999999)
    p.wait(timeout=9999999)

# In[ ]:




