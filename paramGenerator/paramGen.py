import os
import sys
import re
import string

def setupVariableList(fileName,tagTableList,variableList):

	closepatt = re.compile('</.*>.*\n')
	patt = re.compile('<.*>.*\n')
	
	errorflag = 0
	openbracket = 0
	paramPos = 0
	tag = ''

	print 'Processing: ',fileName
	f = open(fileName)
	for line in f:
		ll = line.split()

		if(ll != []):
		
			cm = closepatt.match(line)
			if (cm and openbracket==0):
				errorflag = 1
				break	
			elif (cm and openbracket==1):
				openbracket = 0
				paramPos = 0
				tag = ''
			else:
				m = patt.match(line)
				if(m and openbracket==1):
					errorflag = 1
					break

				elif (m and openbracket==0):	 
					#Extract inner string
					endpos = line.find('>')
					innerset = line[m.start()+1:endpos]
					innerset = innerset.strip()

					openbracket = 1
					pieces = innerset.split()

					#obtain tag name
					pos = 0
					parameterList = []
					for piece in pieces:
						if(pos == 0):
							tag = pieces[0]
						else:	
							parameter = ''
							value = ''	
							#Find field parameters and values
							pos = string.find(piece,'=')
							if(pos != -1):
								parameter = piece[:pos]
								value = piece[pos+1:]
								parameterList.append((parameter,value))
						pos += 1	
										
#					print '##',tag,value		
					tagTableList.append((tag,parameterList))
						
				else:
					#if not illegal i.e openbracket==true and (!m || !cm)
					#then it MUST be variable list
#					print ll
					variableName = tag + '.' + ll[0]
					variableValue = ll[1]
					variableList.append((variableName,[paramPos,variableValue]))
					paramPos += 1	
	f.close()

#########################################################################################################
	
targetDir = sys.argv[1]
parameterFile = sys.argv[2]
print 'Target Dir: ',targetDir

#Create <tag,field value> and <parameter,value> 
tagTableList = []
variableList = []
setupVariableList(parameterFile,tagTableList,variableList)

tagDic = dict(tagTableList)
varDic = dict(variableList)

if(len(varDic)!=len(variableList)):
	print "Error! Duplicate variable name found!"
	sys.exit()

#Now process

fileTable = []
p = 0
#First do substitutions and create array list to match locations
for l in tagDic.keys():
	
	paramlist = tagDic[l]
	
#	print l,paramlist
	empytlist = []
	neededsize = 0
	none = 0
	pos = 0
	for params in paramlist:
	
		if(params[1][0]=='$'):
			varname = params[1][1:]
			subvalue = varDic[varname][1]
		else:
			subvalue = params[1]
	
		tagDic[l][pos] = (tagDic[l][pos][0],subvalue)	
		pos += 1
		
		if(params[0]=='size'):
			neededsize = int(params[1])
		
		if(params[0]=='file' and params[1]=='none'):
			none = 1

	if(none == 1):
		neededsize = 0

#	print l,paramlist,neededsize	
	if(neededsize>0):
		fileTable.append((l,[]))
		for i in range(0,neededsize):
			fileTable[p][1].append('')
	
#		print fileTable[p][0],fileTable[p][1]
		p += 1		

#initialize parameter array list
#print tagDic

fileTableDic = dict(fileTable)

#Now setup parameter at correct locations and matchup sizes etc.
for l in varDic.keys():

	getDotPos = string.find(l,'.')
	tag = l[:getDotPos]	
	parameter = l[getDotPos+1:]

	pos = int(varDic[l][0])
	value = varDic[l][1]

#	print tag,parameter,pos,value

	if (tag in fileTableDic) and (pos<len(fileTableDic[tag])):
			
		fileTableDic[tag][pos] = value
	elif tag in fileTableDic and pos>=len(fileTableDic[tag]):
		print 'Warning!',tag,'defines max-size',len(fileTableDic[tag]),' but pos =',pos,'@ tag =',tag,' found!'
 
#Now do writes
#print fileTableDic
			
for l in fileTableDic:
	
	filename = ''
	outfileFound = 0
	for parampair in tagDic[l]:
		if(parampair[0] == 'file'):
			outfileFound = 1
			filename = parampair[1]
			break

	if(outfileFound == 1):
		wf = open(targetDir + '/' + filename,'w')
		for values in fileTableDic[l]:
			wf.write(values+'\n')		
		wf.close()

