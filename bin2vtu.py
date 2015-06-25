#!/usr/bin/python

# ux*iHat+uy*jHat+uz*kHat

import sys
import struct as st #
#cvw = custom vtu writer
import makevtu as cvw
# Import subprocess library, so system calls can be made
# in this specific case we want to use it to delete old .vtu files
# in case partition size changes etc
import subprocess
import binascii

#"Global constants":
FIN = "output.dat"	#Std. input file format
FOUT = "output"	#Std. output file format

def main(itr):
	print("iter: " + str(itr))
	# Delete all vtu files of same name
	subprocess.call("rm " + FOUT + "_" + str(itr).zfill(5) + ".vtu", shell=True);
	# Try to open the file
	try:
		fin = open(FIN,'rb')
	except:
		exit("Could not open file.. exiting")

	# The file always starts with a user defined string
	# Here it is discarded, but if you have some information in it,
	# it can be saved.
	readInString(fin)

	print("Reading in mesh information")
	#Load information from header
	nDom,nPointsT,nCellsT,nPFields,nCFields,nodesPerElement=readHeader(fin)

	# The cell and point field names:
	pointFieldNames = readInString(fin)
	cellFieldNames  = readInString(fin)
	# Convert to tuples with, assume names are comma separated. Also strip for whitespaces/trailing characters
	pointFieldNames = [x.strip() for x in pointFieldNames.split(',')]
	cellFieldNames = [x.strip() for x in cellFieldNames.split(',')]

	print("Read/write in mesh data: nodes")
	# Open output stream
	try:
		fout = open(FOUT+ "_" + str(itr).zfill(5) + ".vtu",'wb')
	except:
		exit("Cannot create output file... exiting")

	# Node list - simple take and put
	rawP = ""
	for i in range(nDom):
		#Multiply by 3 because we have 3D
		#Multiply by 4 == length of float32
		rawP += fin.read(3*4*nPointsT[i])
	cvw.writeHeader(fout,sum(nPointsT),sum(nCellsT))
	cvw.writeRawPoints(fout,rawP)
	rawP = None # delete rawP from memory

	print("Read/write in mesh data: element connectivity")

	#Cell data needs to be processed and cannot be read in raw
	rawP = ""
	for i in range(nDom):
		#Multiply by 8 because we have 8 nodes per element
                #Multiply by 8 == length of unsigned long int
		rawP += fin.read(8*8*nCellsT[i])
	cvw.writeRawCellsConn(fout,rawP)
	#print st.unpack('Q'*128*8,rawP[0:8*128*8])

	rawP = ""
	for i in range(nDom):
                #Multiply by 8 == length of unsigned long int
                rawP += fin.read(8*nCellsT[i])
	cvw.writeRawCellsOffset(fout,rawP)

	# Convert from binary to e.g. floats and return in a tuple:
	#print st.unpack('Q'*128,rawP[0:8*128])

	rawP = ""
        for i in range(nDom):
                #Multiply by 8 == length of unsigned long int
                rawP += fin.read(8*nCellsT[i])

        cvw.writeRawCellsType(fout,rawP)
	#print st.unpack('Q'*128,rawP[0:8*128])
	rawP=None
	print("Done writing in mesh")


	#Write out a vtu file for each time step
	dataset = 0
	foundRequestedDataset = False
	while(1):
		try:
			iteration = readdata(fin,'Q')
			iteration = iteration[0]
			print("Optimization iter. " + str(iteration) + " = dataset " + str(dataset) + ", you requested dataset " + str(itr))
		except:
			break #break loop

		if int(dataset)==int(itr):
			foundRequestedDataset = True
			print("Processing dataset " + str(dataset))
			lPFieldNames = []
			lCFieldNames = []
			lrawPFields = []
			lrawCFields = []
			for j in range(nPFields[i]):
				lrawPFields.append("")
			for j in range(nCFields[i]):
				lrawCFields.append("")

			for i in range(nDom):
				for j in range(nPFields[i]):
					if(i==0):
						try:
							lPFieldNames.append(pointFieldNames[j])
						except:
							lPFieldNames.append("Point Field " + str(j))
					#Multiply by 4 == length of float32
					lrawPFields[j] += fin.read(4*nPointsT[i])

				for j in range(nCFields[i]):
					if(i==0):
						try:
							lCFieldNames.append(cellFieldNames[j])
						except:
							lCFieldNames.append("Cell Field " + str(j))
					#Multiply by 4 == length of float32
					lrawCFields[j] += fin.read(4*nCellsT[i])
			cvw.writeRawScalarPointData(fout,lrawPFields,lPFieldNames)
			cvw.writeRawScalarCellData(fout,lrawCFields,lCFieldNames)
			cvw.writeFooter(fout)
			fout.close()
		else:
			tmp1 = 0
			for i in range(nDom):
                                for j in range(nPFields[i]):
					#fin.read(4*nPointsT[i])
					tmp1 += 4*nPointsT[i]
				for j in range(nCFields[i]):
					#fin.read(4*nCellsT[i])
					tmp1 += 4*nCellsT[i]
			fin.seek(tmp1,1)
		dataset += 1

	fin.close()
	if foundRequestedDataset:
		print("Done")
	else:
		print("!! The requested dataset was NOT found!! ")
		subprocess.call("rm " + FOUT + "*.vtu", shell=True);




def getNoNodes(i):
	if(i==10):
		return 4
	if(i==12):
		return 8
	if(i==1000):
		return 8
	exit("Sorry, but the element type " + str(i) + " is not defined. You may add it to getNoNodes() yourself... exiting")

#If the input file has specified a custom cell number format they can be added here
#and in the getNoNodes()
def convertToVtkCell(i):
	if(i==1000):
		return 12
	return i


def readdata(fin,inpformat):
	# fin = file input
	# inpformat = type of character/string/number to read - see python manual.
	#  sequence of datatypes

	# How many bytes should we read when format is inpformat:
	bytecount = st.calcsize(inpformat)
	# Read the bytes into tmp
	tmp = fin.read(bytecount)
	# Convert from binary to e.g. floats and return in a tuple:
	return st.unpack(inpformat,tmp)


def readHeader(fin):
	#Should be called right after fin.open()
	try:
		nDom = readdata(fin,'Q')[0]

		tmp = readdata(fin,'Q'*nDom*4)
		nPointsT = list(tmp[0:nDom])
		nCellsT  = list(tmp[nDom:2*nDom])
		nPFields = list(tmp[2*nDom:3*nDom])
		nCFields = list(tmp[3*nDom:4*nDom])

		nodesPerElement = readdata(fin,'Q')[0]
	except:
		exit("Could not read header format... exiting")

	return nDom,nPointsT,nCellsT,nPFields,nCFields,nodesPerElement


def readInString(fin):
# Reads in a string until an end line symbol is detected
	string = ''
	while(1):
		try:
			tmp = readdata(fin,'c')[0] # The c means data of character type (length 1)
			string += tmp
		except:
			exit("File ended while scanning for string. String not present or properly terminated?... exiting")
		# Break if end character detected
		if(tmp == '\x01'):
			string = string[0:-1]; # Dont want to save last character
			break
	return string

# Make sure main is only called when the file is executed
if __name__ == "__main__":
	itr = 0
	if len(sys.argv) > 1:
	        itr = sys.argv[1]

	main(itr)
