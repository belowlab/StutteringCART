import numpy as np



pair_file = open("","r")#Case_control pair file #FORMATING# #CASE\tPAIRED_CONTROL_1\tPARIED_CONTROL_2\t...PAIRED_CONTROL_N
phe_enrich = open("stuttering_extended_phewas_enrichment.csv","r")###Stuttering Enrichment Restuls
GRID_x_phe_file = open("Extended_Stuttering_Binary_PheCode_Counts.csv","r")#File produced from script 1
case_count = 0
for line in pair_file:
    case_count+=1
pair_file = open("","r")#Case_control pair file #CASE\tPAIRED_CONTROL_1\tPARIED_CONTROL_2\t...PAIRED_CONTROL_N #same file as line one
phe_list = []

p_value_threshold = 0.00000001
count_threshold = .035*case_count #frequency limit, adjust as needed





phe_list = []

x = 0
for line in phe_enrich:
    if(x>1):
        spline = line.split(",")
        if(float(spline[3])<p_value_threshold and int(spline[2])>count_threshold):
            phe_list.append(spline[0])
    x+=1

phe_list.sort()
#phecode_index[phecode] = index
#GPD[GRID] = [phecode1,phecode2,...]
phecode_index = {}
GPD = {}
#List of case and controls used in the SD pair set
case_list = []
control_list = []

unp = 0
for line in pair_file:
    line = line.rstrip()
    spline = line.split("\t")
    if(len(spline)>1):
        case_list.append(spline[0])
        for control in spline[1:]:
            control_list.append(control)

ind_x = 0
for phe in phe_list:
    phecode_index[ind_x] = phe
    ind_x+=1

rows = (len(case_list)+len(control_list))
cols = (len(phe_list)+1)
Reg_Table = {}



#Create an empty dictrionary table for regression
#RegTable[GRID] = [status,intercept term,phecode 1,phecode 2, phecode 3, phecode 4, phecode 5...]
#
for case in case_list:
    Reg_Table[case] = []
    (Reg_Table[case]).append("1")#status
    (Reg_Table[case]).append("1")#intercept term
    for thing in phe_list:
        Reg_Table[case].append("0")
for control in control_list:
    Reg_Table[control] = []
    Reg_Table[control].append("0")#status
    Reg_Table[control].append("1")#intercept term
    for thing in phe_list:
        Reg_Table[control].append("0")


#Read file prodcued by 1_script
for line in GRID_x_phe_file:
    spline = (line.rstrip()).split(",")
    GPD[spline[0]] = spline[1:]

missing_case_list = []

m_case = 0
m_con = 0
for case in case_list:
    if(case not in GPD):
        m_case+=1
        missing_case_list.append(case)
for case in control_list:
    if(case not in GPD):
        m_con+=1

print(m_case)
print(m_con)
print(len(Reg_Table))

for thing in missing_case_list:
    del(Reg_Table[thing])

#Now we need to fill the table
#YOU will need phecode_index (maps the index to the code) and GPD (shows what phecodes each GRID has)


for GRID in Reg_Table:
    #(Reg_Table[GRID])[0] = status
    #(Reg_Table[GRID])[1] = intercept term
    x = 2
    for phec in (Reg_Table[GRID])[2:]:
        if(phecode_index[x-2] in GPD[GRID]):
            (Reg_Table[GRID])[x] = "1"
            #print("found one: ", GRID," has phecode ",phecode_index[x])
        x+=1

#Write out table for CART regression scropt
f = open("Extended_Stuttering_CART_Table.csv","w+")
for person in Reg_Table:
    f.write((Reg_Table[person])[0])
    for val in (Reg_Table[person])[1:]:
        f.write(",")
        f.write(val)
    f.write("\n")









