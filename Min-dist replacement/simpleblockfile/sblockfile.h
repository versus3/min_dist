// sblockfile.h
// Simple Block File Class, Declaration
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 16/06/2010
//		Last Modified: 24/09/2012

#ifndef __SBLOCKFILE_H
#define __SBLOCKFILE_H

/* manages a sequentially accessable file */
class SBlockFile{
public:
	SBlockFile(char* szFileName, int nBS, int nIS,	int nPP, int nIT);

	~SBlockFile();

	int getItem(void* dataItem, int index);
	int readBlock(void* data);
	void rewind();

	FILE* m_fp;
	int m_nBlockSize;
	int m_nItemSize;
	int m_nItemPerPage;
	int m_nItemTotal;

	char* m_szBlock;
	long m_nPageAccesses;
};

#endif