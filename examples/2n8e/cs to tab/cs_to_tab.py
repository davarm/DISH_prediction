"""Converting the RCSB .cs file into .tab format"""

"""Maybe just convert to a talos.tab so can then write out proper pipeline?"""

get=open('2n8e_cs.str','r')
get2=open('2n8e.tab','w')
for line in get:
   # line=line.replace('?','')
    #line=line.replace(' . ','')
    line=" ".join(line.split())
    line=line.split(' ')
    new=str(line[3]+' '+line[5]+' '+line[6]+' '+line[9])
    get2.write(new)
    get2.write('\n')