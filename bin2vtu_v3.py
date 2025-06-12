#!/usr/bin/env python3

import sys
import struct as st
import makevtu_v3 as cvw
import subprocess
import binascii

FIN = "output_00000.dat"
FOUT = "output"

def main(itr):
    print(f"iter: {itr}")
    
    # Clean up old .vtu files
    subprocess.call(f"rm {FOUT}*.vtu", shell=True)

    try:
        fin = open(FIN, 'rb')
    except IOError:
        sys.exit("Could not open file.. exiting")

    readInString(fin)  # Discard user-defined string

    print("Reading in mesh information")
    nDom, nPointsT, nCellsT, nPFields, nCFields, nodesPerElement = readHeader(fin)

    pointFieldNames = readInString(fin).split(',')
    cellFieldNames = readInString(fin).split(',')
    pointFieldNames = [x.strip() for x in pointFieldNames]
    cellFieldNames = [x.strip() for x in cellFieldNames]

    print("Read/write in mesh data: nodes")
    try:
        fout = open(f"{FOUT}_{str(itr).zfill(5)}.vtu", 'wb')
    except IOError:
        sys.exit("Cannot create output file... exiting")

    rawP = b""
    for i in range(nDom):
        rawP += fin.read(3 * 4 * nPointsT[i])
    cvw.writeHeader(fout, sum(nPointsT), sum(nCellsT))
    cvw.writeRawPoints(fout, rawP)

    rawP = b""
    for i in range(nDom):
        rawP += fin.read(8 * 8 * nCellsT[i])
    cvw.writeRawCellsConn(fout, rawP)

    rawP = b""
    for i in range(nDom):
        rawP += fin.read(8 * nCellsT[i])
    cvw.writeRawCellsOffset(fout, rawP)

    rawP = b""
    for i in range(nDom):
        rawP += fin.read(8 * nCellsT[i])
    cvw.writeRawCellsType(fout, rawP)

    print("Done writing in mesh")

    dataset = 0
    foundRequestedDataset = False

    while True:
        try:
            iteration = readdata(fin, 'Q')[0]
            print(f"Optimization iter. {iteration} = dataset {dataset}, you requested dataset {itr}")
        except Exception:
            break

        if int(dataset) == int(itr):
            foundRequestedDataset = True
            print(f"Processing dataset {dataset}")
            lPFieldNames, lCFieldNames = [], []
            lrawPFields = [b"" for _ in range(nPFields[0])]
            lrawCFields = [b"" for _ in range(nCFields[0])]

            for i in range(nDom):
                for j in range(nPFields[i]):
                    if i == 0:
                        try:
                            lPFieldNames.append(pointFieldNames[j])
                        except IndexError:
                            lPFieldNames.append(f"Point Field {j}")
                    lrawPFields[j] += fin.read(4 * nPointsT[i])

                for j in range(nCFields[i]):
                    if i == 0:
                        try:
                            lCFieldNames.append(cellFieldNames[j])
                        except IndexError:
                            lCFieldNames.append(f"Cell Field {j}")
                    lrawCFields[j] += fin.read(4 * nCellsT[i])

            cvw.writeRawScalarPointData(fout, lrawPFields, lPFieldNames)
            cvw.writeRawScalarCellData(fout, lrawCFields, lCFieldNames)
            cvw.writeFooter(fout)
            fout.close()
        else:
            offset = sum(4 * nPointsT[i] * nPFields[i] + 4 * nCellsT[i] * nCFields[i] for i in range(nDom))
            fin.seek(offset, 1)

        dataset += 1

    fin.close()
    if foundRequestedDataset:
        print("Done")
    else:
        print("!! The requested dataset was NOT found!! ")
        subprocess.call(f"rm {FOUT}*.vtu", shell=True)


def getNoNodes(i):
    if i == 10:
        return 4
    if i == 12 or i == 1000:
        return 8
    sys.exit(f"Unsupported element type {i}. Extend getNoNodes() to support it.")


def convertToVtkCell(i):
    return 12 if i == 1000 else i


def readdata(fin, inpformat):
    bytecount = st.calcsize(inpformat)
    tmp = fin.read(bytecount)
    return st.unpack(inpformat, tmp)


def readHeader(fin):
    try:
        nDom = readdata(fin, 'Q')[0]
        tmp = readdata(fin, 'Q' * nDom * 4)
        nPointsT = list(tmp[0:nDom])
        nCellsT = list(tmp[nDom:2 * nDom])
        nPFields = list(tmp[2 * nDom:3 * nDom])
        nCFields = list(tmp[3 * nDom:4 * nDom])
        nodesPerElement = readdata(fin, 'Q')[0]
    except Exception:
        sys.exit("Could not read header format... exiting")

    return nDom, nPointsT, nCellsT, nPFields, nCFields, nodesPerElement


def readInString(fin):
    string = b''
    while True:
        try:
            tmp = readdata(fin, 'c')[0]
            if tmp == b'\x01':
                return string[:-1].decode('utf-8', errors='replace')
            string += tmp
        except Exception:
            sys.exit("File ended while scanning for string. Improper termination... exiting")

# Run only if executed directly
if __name__ == "__main__":
    itr = 0
    if len(sys.argv) > 1:
        itr = int(sys.argv[1])
    main(itr)

