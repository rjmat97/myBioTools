import json
import sys
import re

def document():
    print("---------------------------------------------------------------------\n"+
    "documentation:\n"+
    "---------------------------------------------------------------------\n"+
    "this program is meant to compare and find common sequences\n"+
    "from multi fasta files.\n"+
    "It can accept two files as srguments.\n"+
    "the format of the files is very crucial.\n" +
    "accession Id should be first line and consecutive line should be\n"+
    "sequence string. The sequence string should only be one line long. \n"+

    "---------------------------------------------------------------------\n"+
    "if no files are passed as arguments, the program searches\n"+
    "for the files \"small.fasta\" and \"big.fasta\"\n"+
    "---------------------------------------------------------------------\n")

def file2JSON(file):
    json_ret = {}
    try:
        with open(file) as f:
            content = f.readlines()
            for i in range(0, len(content), 2):
                json_ret[content[i]] = content[i+1]
    except:
        print("an ERROR occured---------------------------")
        print(file + "file does not exist!")
        print("exiting progream")
        print("-------------------------------------------")
        exit()
    return json_ret

def save2JSON(file):
    json_in = file2JSON(file)
    save_name = re.sub("\..*", ".json", file)
    with open(save_name, 'w') as f:
        json.dump(json_in, f)

def compJSONs(file1, file2):
    json1 = file2JSON(file1)
    json2 = file2JSON(file2)
    comKeys = 0
    for key in json1.keys():
        seq1 = json1[key]
        for inKey in json2.keys():
            if(json2[inKey] == seq1): 
                print(str(key)+" and "+str(inKey)+" have the same sequence")
                comKeys+=1
    print("-------------------------------------------")
    print("comparison complete!")
    if(comKeys == 0):
        print("no common sequences found")
    else:
        print(str(comKeys)+" number of common sequences found")
    print("-------------------------------------------")


if(len(sys.argv) == 3):
    small_file = sys.argv[1]
    big_file   = sys.argv[2]
    compJSONs(small_file, big_file)
elif(len(sys.argv) == 2):
    if(sys.argv[1] == "-h" or sys.argv[1] == "help"):
        document()
    else:
        print("-------------------------------------------")
        print("save 2 JSON mode activated:")
        save2JSON(sys.argv[1])
        print("file saved successfully")
        print("-------------------------------------------")
else:
    small_file = small.fasta
    big_file   = big.fasta
    compJSONs(small_file, big_file)