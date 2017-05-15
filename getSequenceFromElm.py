f = open("phosphoELM_all_2015-04.dump", "r")
strTemp = f.readline()
result = []

stringResult = ""
pre_accID = ""
while True:
    strTemp = f.readline()
    if strTemp == "":
        break

    strList = strTemp.split("\t")
    # print strList
    strList = strList[0:3]

    if pre_accID != strList[0]:
        pre_accID = strList[0]
        # print(strList[1])
        f_seq_out = open("./Sequences_from_elm/" + strList[0], "w")

        strTemp = ">" + strList[0] + "\t" + strList[2] + "\n" + strList[1] + "\n"
        f_seq_out.write(strTemp)

        # result.append(strList)
        result.append(strTemp)

        stringResult += strTemp