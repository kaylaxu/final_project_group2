import cptac
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
cptac.download(dataset="Gbm")
gbm = cptac.Gbm()
protein_data = gbm.get_proteomics()
rna_data = gbm.get_transcriptomics()
clinical_data = gbm.get_clinical()
ok = np.intersect1d(protein_data.index, clinical_data.index)
ok_protein = protein_data.loc[ok,:]
ok_clinic = clinical_data.loc[ok,:]
alive_boolmask = (ok_clinic.loc[:, "cause_of_death"]).isnull()
dead_boolmask = ok_clinic.loc[:,"cause_of_death"].notnull()
protein_alive = ok_protein.loc[alive_boolmask,:]
protein_dead = ok_protein.loc[dead_boolmask,:]

###PTEN
pten_alive = protein_alive.loc[:, "PTEN"]
pten_dead = protein_dead.loc[:, "PTEN"]
pten = [pten_alive, pten_dead]

# Creating plot
plt.boxplot(pten)
plt.xticks([1, 2], ["Alive", "Dead"])
plt.title("PTEN")
plt.savefig('ABSOLUTE/PATH/TO/qbio490_name/final_project_group2/outputs/PTEN_protein_boxplot')

###TP53
tp53_alive = protein_alive.loc[:, "TP53"]
tp53_dead = protein_dead.loc[:, "TP53"]
mask = tp53_alive.notnull()
tp53_alive = tp53_alive.loc[mask]
mask = tp53_dead.notnull()
tp53_dead = tp53_dead.loc[mask]
tp53 = [tp53_alive, tp53_dead]
# Creating plot
plt.boxplot(tp53)
plt.xticks([1, 2], ["Alive", "Dead"])
plt.title("TP53")
plt.savefig('ABSOLUTE/PATH/TO/qbio490_name/final_project_group2/outputs/TP53_protein_boxplot')

###MME
MME_alive = protein_alive.loc[:, "MME"]
MME_dead = protein_dead.loc[:, "MME"]
mask = MME_alive.notnull()
MME_alive = MME_alive.loc[mask]
mask = MME_dead.notnull()
MME_dead = MME_dead.loc[mask]
MME = [MME_alive, MME_dead]
plt.boxplot(MME)
plt.xticks([1, 2], ["Alive", "Dead"])
plt.title("MME")
plt.savefig('ABSOLUTE/PATH/TO/qbio490_name/final_project_group2/outputs/MME_protein_boxplot')

### GLI1
GLI1_alive = protein_alive.loc[:, "GLI1"]
GLI1_dead = protein_dead.loc[:, "GLI1"]
mask = GLI1_dead.notnull()
GLI1_dead = GLI1_dead.loc[mask]
mask = GLI1_alive.notnull()
GLI1_alive = GLI1_alive.loc[mask]
GLI1 = [GLI1_alive, GLI1_dead]
plt.boxplot(GLI1)
plt.xticks([1,2],["Alive", "Dead"])
plt.title("GLI1")
plt.savefig('ABSOLUTE/PATH/TO/qbio490_name/final_project_group2/outputs/GLI1_protein_boxplot')


