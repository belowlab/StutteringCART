import sys
import warnings
import math
import statsmodels
import numpy as np
from scipy import stats
import statsmodels.formula.api as smf

phe_enrich = open("stuttering_extended_phewas_enrichment.csv","r")#ENTER ENRICHMENT RESULTS TABLE
phe_list = []
p_value_threshold = 0.0001
count_threshold = 19
phecode_index = {}
phe_ind_dic = {}
xin = 0
for line in phe_enrich:
    if(xin>1):
        spline = line.split(",")
        if(float(spline[3])<p_value_threshold and int(spline[2])>count_threshold):
            phe_list.append(spline[0])
            phe_ind_dic[xin-1] = spline[0]
    xin+=1
phe_list.sort()
ind_x = 0
for phe in phe_list:
    phecode_index[ind_x] = phe
    ind_x+=1
#print(phecode_index)

firth_table = open("Extended_Stuttering_CART_Table.csv","r")
status_list = []
parameter_list = []
for line in firth_table:
    spline = (line.rstrip()).split(",")
    status_list.append(int(spline[0]))
    indp_list = []
    for phe in spline[2:]:#Skipping the first 2 componants (first is status, second is intercept term used for regression ML model from before
        indp_list.append(int(phe))
    parameter_list.append(indp_list)
X = np.array(parameter_list)
y =np.array(status_list)


from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=1)
from sklearn import tree

model = tree.DecisionTreeClassifier(min_samples_leaf=5)
print(model)
model.fit(X_train, y_train)

y_predict = model.predict(X_test)

from sklearn.metrics import accuracy_score

print("Overall model accuracy:",accuracy_score(y_test, y_predict))

from sklearn.metrics import confusion_matrix

print(confusion_matrix(y_test, y_predict))
tn, fp, fn, tp = confusion_matrix(y_test,y_predict).ravel()
print("tn:",tn)
print("fp:",fp)
print("fn:",fn)
print("tp:",tp)
print("PPV = ",(tp/(tp+fp)))

print(y_predict)


BioVU_phe_file = open("inputfilehere","r")#PHECODE FILE FOR PATIENT PREDICTION #FORMATING# #ID,PHECODE
#phe_list_file = open("./SP1_Regression/SP1_Stuttering_PhecodesxEffectSize.csv","r")

#Get the phecodes for all the GRIDs so we can construct out firthtalbe
GRID_phe_dic = {}
x = 0
print("Getting phecodes for BioVU patients... this will take awhile")
for line in BioVU_phe_file:
    if(x>0):
        spline = (line.rstrip()).split(",")
        if(spline[0] not in GRID_phe_dic):
            GRID_phe_dic[spline[0]] = [spline[1]]
        else:
            GRID_phe_dic[spline[0]].append(spline[1])
    x+=1
    #if(x>100000):
        #break
####Making the table to feed into the prediction

total_BV_table = []
GRID_index_reference_list = []
for GRID in GRID_phe_dic:
    GRID_index_reference_list.append(GRID)
    ind_list_tab = []
    for inx in phe_ind_dic:
        if(phe_ind_dic[inx] in GRID_phe_dic[GRID]):
            ind_list_tab.append(1)
        else:
            ind_list_tab.append(0)
    total_BV_table.append(ind_list_tab)
BioVU_X = np.array(total_BV_table)
BioVU_y_predict = model.predict(BioVU_X)

pred_ind = 0
positive_hits = []
for person in BioVU_y_predict:
    if(person==1):
        positive_hits.append(GRID_index_reference_list[pred_ind])
    pred_ind+=1

print("Total identified patients: ",len(positive_hits))

f = open("BioVU_patients_predicted_to_stutter.txt","w+")

for biovu_pat in positive_hits:
    f.write(biovu_pat)
    f.write("\n")

f.close()
