#!/usr/bin/env python3

import struct as st
import base64

FOUT = "vtutest.vtu"

def main():
    try:
        fout = open(FOUT, 'wb')
    except IOError:
        exit("Couldn't create file...")

    writeHeader(fout, 4, 1)
    lP = [0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 1, 1]
    writePoints(fout, lP)
    lConn = [0, 1, 2, 3]
    lOffset = [4]
    lTypes = [10]
    writeCells(fout, lConn, lOffset, lTypes)

    lldata = []
    lldata.append(["Scalars", 1, 2, -1, 0.5])
    lldata.append(["OtherSet", 1.0, 2, 3, 127.3])
    writeScalarPointData(fout, lldata)

    lldata = []
    lldata.append(["Scalars", 1.0])
    lldata.append(["Scalars2", 27.0])
    writeScalarCellData(fout, lldata)

    writeFooter(fout)
    fout.close()


def writeHeader(fout, nP, nC):
    fout.write(b'<?xml version="1.0"?>\n')
    fout.write(b'<VTKFile type="UnstructuredGrid" version="1.0" header_type="UInt64" byte_order="LittleEndian">\n')
    fout.write(b'<UnstructuredGrid>\n')
    fout.write(f'\t<Piece NumberOfPoints="{nP}" NumberOfCells="{nC}">\n'.encode('utf-8'))


def writeFooter(fout):
    fout.write(b'\t</Piece>\n')
    fout.write(b'</UnstructuredGrid>\n')
    fout.write(b'</VTKFile>')


def writePoints(fout, lP):
    fout.write(b'\t\t<Points>\n')
    fout.write(b'\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="binary">\n')
    data = b''.join(st.pack('f', number) for number in lP)
    fout.write(writeBin64(data))
    fout.write(b'\t\t\t</DataArray>\n')
    fout.write(b'\t\t</Points>\n')


def writeCells(fout, lConn, lOffset, lTypes):
    fout.write(b'\t\t<Cells>\n')

    fout.write(b'\t\t\t<DataArray type="UInt32" Name="connectivity" format="binary">\n')
    data = b''.join(st.pack('I', number) for number in lConn)
    fout.write(writeBin64(data))
    fout.write(b'\t\t\t</DataArray>\n')

    fout.write(b'\t\t\t<DataArray type="UInt32" Name="offsets" format="binary">\n')
    data = b''.join(st.pack('I', number) for number in lOffset)
    fout.write(writeBin64(data))
    fout.write(b'\t\t\t</DataArray>\n')

    fout.write(b'\t\t\t<DataArray type="UInt32" Name="types" format="binary">\n')
    data = b''.join(st.pack('I', number) for number in lTypes)
    fout.write(writeBin64(data))
    fout.write(b'\t\t\t</DataArray>\n')

    fout.write(b'\t\t</Cells>\n')


def writeBin64(bindata):
    datalen = st.pack('Q', len(bindata))
    return base64.b64encode(datalen) + base64.b64encode(bindata)


def writeScalarPointData(fout, lldata):
    fout.write(b'\t\t<PointData Scalars="scalars">\n')
    for ldata in lldata:
        setname = ldata[0]
        values = ldata[1:]
        data = b''.join(st.pack('f', v) for v in values)
        fout.write(f'\t\t\t<DataArray type="Float32" Name="{setname}" format="binary">\n'.encode('utf-8'))
        fout.write(writeBin64(data))
        fout.write(b'\t\t\t</DataArray>\n')
    fout.write(b'\t\t</PointData>\n')


def writeRawScalarPointData(fout, ldata, setnames):
    fout.write(b'\t\t<PointData Scalars="scalars">\n')
    for name, data in zip(setnames, ldata):
        fout.write(f'\t\t\t<DataArray type="Float32" Name="{name}" format="binary">\n'.encode('utf-8'))
        fout.write(writeBin64(data))
        fout.write(b'\t\t\t</DataArray>\n')
    fout.write(b'\t\t</PointData>\n')


def writeScalarCellData(fout, lldata):
    fout.write(b'\t\t<CellData Scalars="scalars">\n')
    for ldata in lldata:
        setname = ldata[0]
        values = ldata[1:]
        data = b''.join(st.pack('f', v) for v in values)
        fout.write(f'\t\t\t<DataArray type="Float32" Name="{setname}" format="binary">\n'.encode('utf-8'))
        fout.write(writeBin64(data))
        fout.write(b'\t\t\t</DataArray>\n')
    fout.write(b'\t\t</CellData>\n')


def writeRawScalarCellData(fout, ldata, setnames):
    fout.write(b'\t\t<CellData Scalars="scalars">\n')
    for name, data in zip(setnames, ldata):
        fout.write(f'\t\t\t<DataArray type="Float32" Name="{name}" format="binary">\n'.encode('utf-8'))
        fout.write(writeBin64(data))
        fout.write(b'\t\t\t</DataArray>\n')
    fout.write(b'\t\t</CellData>\n')


# Optional raw functions for compressed bulk data:
def writeRawPoints(fout, data):
    fout.write(b'\t\t<Points>\n')
    fout.write(b'\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="binary">\n')
    fout.write(writeBin64(data))
    fout.write(b'\t\t\t</DataArray>\n')
    fout.write(b'\t\t</Points>\n')


def writeRawCellsConn(fout, data):
    fout.write(b'\t\t<Cells>\n')
    fout.write(b'\t\t\t<DataArray type="UInt64" Name="connectivity" format="binary">\n')
    fout.write(writeBin64(data))
    fout.write(b'\t\t\t</DataArray>\n')


def writeRawCellsOffset(fout, data):
    fout.write(b'\t\t\t<DataArray type="UInt64" Name="offsets" format="binary">\n')
    fout.write(writeBin64(data))
    fout.write(b'\t\t\t</DataArray>\n')


def writeRawCellsType(fout, data):
    fout.write(b'\t\t\t<DataArray type="UInt64" Name="types" format="binary">\n')
    fout.write(writeBin64(data))
    fout.write(b'\t\t\t</DataArray>\n')
    fout.write(b'\t\t</Cells>\n')


if __name__ == "__main__":
    main()

