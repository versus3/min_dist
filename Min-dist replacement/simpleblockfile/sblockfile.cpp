#pragma warning(disable:4996)
// sblockfile.cpp
// Simple Block File Class, Implementation
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 16/06/2010
//		Last Modified: 24/09/2012

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sblockfile.h"


SBlockFile::SBlockFile(char* szFileName, int nBS, int nIS,	int nPP, int nIT) {
	if((m_fp = fopen(szFileName, "rb")) == NULL){
		printf("Can not open file %s\n", szFileName);
		exit(0);
	}

	m_nBlockSize = nBS;
	m_nItemSize = nIS;
	m_nItemPerPage = nPP;
	m_nItemTotal = nIT;

	m_szBlock = new char[m_nBlockSize];
	m_nPageAccesses = 0;
}

SBlockFile::~SBlockFile() {
	delete[] m_szBlock;
	m_szBlock = NULL;
	fclose(m_fp);
}


int SBlockFile::getItem(void* dataItem, int index) {
	if(index < 0 || index > m_nItemTotal){
		return 0;
	}

	int nOffset = m_nBlockSize * (index / m_nItemPerPage) + sizeof(int) + (index % m_nItemPerPage)*m_nItemSize;
	fseek(m_fp, nOffset, SEEK_SET);
	fread(dataItem, m_nItemSize, 1, m_fp);
	m_nPageAccesses++;
	return 1;
}

int SBlockFile::readBlock(void* data) {
	int nRecordCount = 0;

	fread((void*)m_szBlock, sizeof(char), m_nBlockSize, m_fp);
	memcpy((void*)&nRecordCount, m_szBlock, sizeof(int));
	memcpy(data, m_szBlock + sizeof(int), nRecordCount*m_nItemSize);

	m_nPageAccesses++;
	return nRecordCount;
}

void SBlockFile::rewind() {
	fseek(m_fp, 0, SEEK_SET);
}
