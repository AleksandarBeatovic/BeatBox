#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os.path
def csv_read(path_to_csv_file, delim=','):
    if not os.path.exists(path_to_csv_file):
        print("Error, such file doesn't exist")
        return []
    else:
        csv = open(path_to_csv_file, "r")
        reads = csv.readlines()
        lines = []
        for line in reads:
            lines.append(line[1:len(line)-2])
               
        lines_1 = []
        for line in lines:
            lines_1.append(line.split(delim))
            
        
        for line in lines_1:
            for i in line:
                if i[0] == '"' and i[-1] != '"':
                    line[line.index(i)] = i + delim + line[line.index(i)+1]
             
            for i in line:
                if i[0] != '"' and i[-1] == '"':
                    del line[line.index(i)]
            
        return (lines_1)
        csv.close()
def write_csv(path_to_csv_file, data, delim =','):
    file = open(path_to_csv_file,"w")
    
    lines_str = []
    for line in data:
        x = delim.join(line)
        lines_str.append(x)
             
            
        
           
    for i in lines_str:
        file.write(i + "\n")
    file.close()


# In[ ]:




