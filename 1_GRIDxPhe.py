pair_file = open("","r")#Case_control pair file #FORMATING# #CASE\tPAIRED_CONTROL_1\tPARIED_CONTROL_2\t...PAIRED_CONTROL_N
phe_enrich = open("stuttering_extended_phewas_enrichment.csv","r")###Stuttering Enrichment Restuls
phe_file = open("","r")#Phecode file, #FORMATING# #ID,PHECODE,PHECODENAME,PHECODEDATE
case_count = 0
for line in pair_file:
    case_count+=1
pair_file = open("","r")#Case_control pair file #CASE\tPAIRED_CONTROL_1\tPARIED_CONTROL_2\t...PAIRED_CONTROL_N #same file as line one
phe_list = []

p_value_threshold = 0.00000001
count_threshold = .035*case_count #frequency limit, adjust as needed


phe_enrich.readline()#HEADER
for line in phe_enrich:
    spline = line.split(",")
    if(float(spline[3])<p_value_threshold and int(spline[2])>count_threshold):
        print(spline[3])
        phe_list.append(spline[0])

phe_list.sort()
#GRID_dic[GRID] = [status,[phecode status]]
#phecode_index[phecode] = index
#GPD[GRID] = [phecode1,phecode2,...]
phecode_index = {}
GRID_dic = {}
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
    phecode_index[phe] = ind_x
    ind_x+=1







line_n = 0
for line in phe_file:
    spline = line.split(",")
    if(spline[0] in case_list or spline[0] in control_list):
        if(spline[0] not in GPD):
            GPD[spline[0]] = []
            GPD[spline[0]].append(spline[1])
        else:
            GPD[spline[0]].append(spline[1])
    line_n+=1

f = open("Extended_Stuttering_Binary_PheCode_Counts.csv","w+")#WRITE OUT BINARY PHECODE STATUS
for grid in GPD:
    f.write(grid)
    for phe in GPD[grid]:
        f.write(",")
        f.write(phe)
    f.write("\n")
    

