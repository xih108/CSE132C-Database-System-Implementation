/**
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include <vector>
#include "btree.h"
#include "page.h"
#include "filescan.h"
#include "page_iterator.h"
#include "file_iterator.h"
#include "exceptions/insufficient_space_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/end_of_file_exception.h"
#define extra_test 1

#define checkPassFail(a, b) 																				\
{																																		\
	if(a == b)																												\
		std::cout << "\nTest passed at line no:" << __LINE__ << "\n";		\
	else																															\
	{																																	\
		std::cout << "\nTest FAILS at line no:" << __LINE__;						\
		std::cout << "\nExpected no of records:" << b;									\
		std::cout << "\nActual no of records found:" << a;							\
		std::cout << std::endl;																					\
		exit(1);																												\
	}																																	\
}

using namespace badgerdb;

// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------
int testNum = 1;
const std::string relationName = "relA";
//If the relation size is changed then the second parameter 2 chechPassFail may need to be changed to number of record that are expected to be found during the scan, else tests will erroneously be reported to have failed.
const int	relationSize = 5000;
std::string intIndexName, doubleIndexName, stringIndexName;

// This is the structure for tuples in the base relation

typedef struct tuple {
	int i;
	double d;
	char s[64];
} RECORD;

PageFile* file1;
RecordId rid;
RECORD record1;
std::string dbRecord1;

BufMgr * bufMgr = new BufMgr(100);

// -----------------------------------------------------------------------------
// Forward declarations
// -----------------------------------------------------------------------------

void createRelationForward();
void createRelationBackward();
void createRelationRandom();
void createRelationSparse();
void createRelationNegative();
void createRelationEmpty();

void indexTests();
void indexTests_sparse();
void indexTests_negative();
void indexTests_empty();

void intTests();
void intTests_sparse();
void intTests_negative();
void intTests_empty();
int intScan(BTreeIndex *index, int lowVal, Operator lowOp, int highVal, Operator highOp);
void doubleTests();
void doubleTests_sparse();
void doubleTests_negative();
void doubleTests_empty();
int doubleScan(BTreeIndex *index, double lowVal, Operator lowOp, double highVal, Operator highOp);
void stringTests();
void stringTests_sparse();
void stringTests_negative();
void stringTests_empty();
int stringScan(BTreeIndex *index, int lowVal, Operator lowOp, int highVal, Operator highOp);
void test1();
void test2();
void test3();
//extra tests
void test4();
void test5();
void test6();
void errorTests();
void deleteRelation();

int main(int argc, char **argv)
{
	if( argc != 2 )
	{
		std::cout << "Expects one argument as a number between 1 to 3 to choose datatype of key.\n";
		std::cout << "For INTEGER keys run as: ./badgerdb_main 1\n";
		std::cout << "For DOUBLE keys run as: ./badgerdb_main 2\n";
		std::cout << "For STRING keys run as: ./badgerdb_main 3\n";
		return 0;
	}

	sscanf(argv[1],"%d",&testNum);

	switch(testNum)
	{
		case 1:
			std::cout << "leaf size:" << INTARRAYLEAFSIZE << " non-leaf size:" << INTARRAYNONLEAFSIZE << std::endl;
			break;
		case 2:
			std::cout << "leaf size:" << DOUBLEARRAYLEAFSIZE << " non-leaf size:" << DOUBLEARRAYNONLEAFSIZE << std::endl;
			break;
		case 3:
			std::cout << "leaf size:" << STRINGARRAYLEAFSIZE << " non-leaf size:" << STRINGARRAYNONLEAFSIZE << std::endl;
			break;
	}


  // Clean up from any previous runs that crashed.
  try
	{
    File::remove(relationName);
  }
	catch(FileNotFoundException)
	{
  }

	{
		// Create a new database file.
		PageFile new_file = PageFile::create(relationName);

		// Allocate some pages and put data on them.
		for (int i = 0; i < 20; ++i)
		{
			PageId new_page_number;
			Page new_page = new_file.allocatePage(new_page_number);

    	sprintf(record1.s, "%05d string record", i);
    	record1.i = i;
    	record1.d = (double)i;
    	std::string new_data(reinterpret_cast<char*>(&record1), sizeof(record1));

			new_page.insertRecord(new_data);
			new_file.writePage(new_page_number, new_page);
		}

	}
	// new_file goes out of scope here, so file is automatically closed.

	{
		FileScan fscan(relationName, bufMgr);

		try
		{
			RecordId scanRid;
			while(1)
			{
				fscan.scanNext(scanRid);
				//Assuming RECORD.i is our key, lets extract the key, which we know is INTEGER and whose byte offset is also know inside the record. 
				std::string recordStr = fscan.getRecord();
				const char *record = recordStr.c_str();
				int key = *((int *)(record + offsetof (RECORD, i)));
				std::cout << "Extracted : " << key << std::endl;
			}
		}
		catch(EndOfFileException e)
		{
			std::cout << "Read all records" << std::endl;
		}
	}
	// filescan goes out of scope here, so relation file gets closed.

	File::remove(relationName);

	test1();
	test2();
	test3();

	#ifdef extra_test
		test4();
		test5();
		test6();
	#endif
	errorTests();

  return 1;
}

void test1()
{
	// Create a relation with tuples valued 0 to relationSize and perform index tests 
	// on attributes of all three types (int, double, string)
	std::cout << "---------------------" << std::endl;
	std::cout << "createRelationForward" << std::endl;
	createRelationForward();
	indexTests();
	deleteRelation();
}

void test2()
{
	// Create a relation with tuples valued 0 to relationSize in reverse order and perform index tests 
	// on attributes of all three types (int, double, string)
	std::cout << "----------------------" << std::endl;
	std::cout << "createRelationBackward" << std::endl;
	createRelationBackward();
	indexTests();
	deleteRelation();
}

void test3()
{
	// Create a relation with tuples valued 0 to relationSize in random order and perform index tests 
	// on attributes of all three types (int, double, string)
	std::cout << "--------------------" << std::endl;
	std::cout << "createRelationRandom" << std::endl;
	createRelationRandom();
	indexTests();
	deleteRelation();
}

void test4()
{
	std::cout << "--------------------" << std::endl;
	std::cout << "createRelationSparse" << std::endl;
	createRelationSparse();
	indexTests_sparse();
	deleteRelation();
}

void test5()
{
	std::cout << "--------------------" << std::endl;
	std::cout << "createRelationNegative" << std::endl;
	createRelationNegative();
	indexTests_negative();
	deleteRelation();
}

void test6()
{
	std::cout << "--------------------" << std::endl;
	std::cout << "createRelationEmpty" << std::endl;
	createRelationEmpty();
	indexTests_empty();
	deleteRelation();
}

// -----------------------------------------------------------------------------
// createRelationForward
// -----------------------------------------------------------------------------

void createRelationForward()
{
	std::vector<RecordId> ridVec;
  // destroy any old copies of relation file
	try
	{
		File::remove(relationName);
	}
	catch(FileNotFoundException e)
	{
	}

  file1 = new PageFile(relationName, true);

  // initialize all of record1.s to keep purify happy
  memset(record1.s, ' ', sizeof(record1.s));
	PageId new_page_number;
  Page new_page = file1->allocatePage(new_page_number);

  // Insert a bunch of tuples into the relation.
  for(int i = 0; i < relationSize; i++ )
	{
    sprintf(record1.s, "%05d string record", i);
    record1.i = i;
    record1.d = (double)i;
    std::string new_data(reinterpret_cast<char*>(&record1), sizeof(record1));

		while(1)
		{
			try
			{
    		new_page.insertRecord(new_data);
				break;
			}
			catch(InsufficientSpaceException e)
			{
				file1->writePage(new_page_number, new_page);
  			new_page = file1->allocatePage(new_page_number);
			}
		}
  }

	file1->writePage(new_page_number, new_page);
}

// -----------------------------------------------------------------------------
// createRelationBackward
// -----------------------------------------------------------------------------

void createRelationBackward()
{
  // destroy any old copies of relation file
	try
	{
		File::remove(relationName);
	}
	catch(FileNotFoundException e)
	{
	}
  file1 = new PageFile(relationName, true);

  // initialize all of record1.s to keep purify happy
  memset(record1.s, ' ', sizeof(record1.s));
	PageId new_page_number;
  Page new_page = file1->allocatePage(new_page_number);

  // Insert a bunch of tuples into the relation.
  for(int i = relationSize - 1; i >= 0; i-- )
	{
    sprintf(record1.s, "%05d string record", i);
    record1.i = i;
    record1.d = i;

    std::string new_data(reinterpret_cast<char*>(&record1), sizeof(RECORD));

		while(1)
		{
			try
			{
    		new_page.insertRecord(new_data);
				break;
			}
			catch(InsufficientSpaceException e)
			{
				file1->writePage(new_page_number, new_page);
  			new_page = file1->allocatePage(new_page_number);
			}
		}
  }

	file1->writePage(new_page_number, new_page);
}

// -----------------------------------------------------------------------------
// createRelationRandom
// -----------------------------------------------------------------------------

void createRelationRandom()
{
  // destroy any old copies of relation file
	try
	{
		File::remove(relationName);
	}
	catch(FileNotFoundException e)
	{
	}
  file1 = new PageFile(relationName, true);

  // initialize all of record1.s to keep purify happy
  memset(record1.s, ' ', sizeof(record1.s));
	PageId new_page_number;
  Page new_page = file1->allocatePage(new_page_number);

  // insert records in random order

  std::vector<int> intvec(relationSize);
  for( int i = 0; i < relationSize; i++ )
  {
    intvec[i] = i;
  }

  long pos;
  int val;
	int i = 0;
  while( i < relationSize )
  {
    pos = random() % (relationSize-i);
    val = intvec[pos];
    sprintf(record1.s, "%05d string record", val);
    record1.i = val;
    record1.d = val;

    std::string new_data(reinterpret_cast<char*>(&record1), sizeof(RECORD));

		while(1)
		{
			try
			{
    		new_page.insertRecord(new_data);
				break;
			}
			catch(InsufficientSpaceException e)
			{
      	file1->writePage(new_page_number, new_page);
  			new_page = file1->allocatePage(new_page_number);
			}
		}

		int temp = intvec[relationSize-1-i];
		intvec[relationSize-1-i] = intvec[pos];
		intvec[pos] = temp;
		i++;
  }
  
	file1->writePage(new_page_number, new_page);
}


// -----------------------------------------------------------------------------
// createRelationSparse
// -----------------------------------------------------------------------------

void createRelationSparse()
{
  // destroy any old copies of relation file
	try
	{
		File::remove(relationName);
	}
	catch(FileNotFoundException e)
	{
	}
  file1 = new PageFile(relationName, true);

  // initialize all of record1.s to keep purify happy
  memset(record1.s, ' ', sizeof(record1.s));
  PageId new_page_number;
  Page new_page = file1->allocatePage(new_page_number);

  // insert records in random order

  std::vector<int> intvec(relationSize);
  for( int i = 0; i < relationSize; i++ )
  {
    intvec[i] = 5*i;
  }

  long pos;
  int val;
	int i = 0;
  while( i < relationSize )
  {
    pos = random() % (relationSize-i);
    val = intvec[pos];
    sprintf(record1.s, "%05d string record", val);
    record1.i = val;
    record1.d = val;

    std::string new_data(reinterpret_cast<char*>(&record1), sizeof(RECORD));

		while(1)
		{
			try
			{
    		new_page.insertRecord(new_data);
				break;
			}
			catch(InsufficientSpaceException e)
			{
      	file1->writePage(new_page_number, new_page);
  			new_page = file1->allocatePage(new_page_number);
			}
		}

		int temp = intvec[relationSize-1-i];
		intvec[relationSize-1-i] = intvec[pos];
		intvec[pos] = temp;
		i++;
  }
  
	file1->writePage(new_page_number, new_page);
}


// -----------------------------------------------------------------------------
// createRelationNNegative
// -----------------------------------------------------------------------------

void createRelationNegative()
{
  // destroy any old copies of relation file
	try
	{
		File::remove(relationName);
	}
	catch(FileNotFoundException e)
	{
	}
  file1 = new PageFile(relationName, true);

  // initialize all of record1.s to keep purify happy
  memset(record1.s, ' ', sizeof(record1.s));
  PageId new_page_number;
  Page new_page = file1->allocatePage(new_page_number);

  // insert records in random order

  std::vector<int> intvec(relationSize);
  for( int i = 0; i < relationSize; i++ )
  {
	if (i % 2 == 0){
		intvec[i] = i / 2;
	}
    else{
		intvec[i] = - (i+1) / 2;
	}
  }

  long pos;
  int val;
	int i = 0;
  while( i < relationSize )
  {
    pos = random() % (relationSize-i);
    val = intvec[pos];
    sprintf(record1.s, "%05d string record", val);
    record1.i = val;
    record1.d = val;

    std::string new_data(reinterpret_cast<char*>(&record1), sizeof(RECORD));

		while(1)
		{
			try
			{
    		new_page.insertRecord(new_data);
				break;
			}
			catch(InsufficientSpaceException e)
			{
      	file1->writePage(new_page_number, new_page);
  			new_page = file1->allocatePage(new_page_number);
			}
		}

		int temp = intvec[relationSize-1-i];
		intvec[relationSize-1-i] = intvec[pos];
		intvec[pos] = temp;
		i++;
  }
  
	file1->writePage(new_page_number, new_page);
}

void createRelationEmpty()
{
	try
	{
		File::remove(relationName);
	}
	catch(FileNotFoundException e)
	{
	}
  file1 = new PageFile(relationName, true);	
}

// -----------------------------------------------------------------------------
// indexTests
// -----------------------------------------------------------------------------

void indexTests()
{
  if(testNum == 1)
  {
    intTests();
		try
		{
			File::remove(intIndexName);
		}
  	catch(FileNotFoundException e)
  	{
  	}
  }
  else if(testNum == 2)
  {
    doubleTests();
		try
		{
			File::remove(doubleIndexName);
		}
  	catch(FileNotFoundException e)
  	{
  	}
  }
  else if(testNum == 3)
  {
    stringTests();
		try
		{
			File::remove(stringIndexName);
		}
  	catch(FileNotFoundException e)
  	{
  	}
  }
}

// -----------------------------------------------------------------------------
// indexTestsSparse
// -----------------------------------------------------------------------------
void indexTests_sparse()
{
  if(testNum == 1)
  {
    intTests_sparse();
		try
		{
			File::remove(intIndexName);
		}
  	catch(FileNotFoundException e)
  	{
  	}
  }
  else if(testNum == 2)
  {
    doubleTests_sparse();
		try
		{
			File::remove(doubleIndexName);
		}
  	catch(FileNotFoundException e)
  	{
  	}
  }
  else if(testNum == 3)
  {
    stringTests_sparse();
		try
		{
			File::remove(stringIndexName);
		}
  	catch(FileNotFoundException e)
  	{
  	}
  }
}

// -----------------------------------------------------------------------------
// indexTestsNegative
// -----------------------------------------------------------------------------
void indexTests_negative()
{
  if(testNum == 1)
  {
    intTests_negative();
		try
		{
			File::remove(intIndexName);
		}
  	catch(FileNotFoundException e)
  	{
  	}
  }
  else if(testNum == 2)
  {
    doubleTests_negative();
		try
		{
			File::remove(doubleIndexName);
		}
  	catch(FileNotFoundException e)
  	{
  	}
  }
  else if(testNum == 3)
  {
    stringTests_negative();
		try
		{
			File::remove(stringIndexName);
		}
  	catch(FileNotFoundException e)
  	{
  	}
  }
}

void indexTests_empty() 
{
	if (testNum == 1) 
	{
		intTests_empty();
		try {
			File::remove(intIndexName);
		}
		catch(FileNotFoundException e) {}
	}
	else if (testNum == 2) 
	{
		doubleTests_empty();
		try {
			File::remove(doubleIndexName);
		}
		catch(FileNotFoundException e) {}
	}
	else if (testNum == 3) 
	{
		stringTests_empty();
		try {
			File::remove(stringIndexName);
		}
		catch(FileNotFoundException e) {}
	}
}
// -----------------------------------------------------------------------------
// intTests
// -----------------------------------------------------------------------------

void intTests()
{
  std::cout << "Create a B+ Tree index on the integer field" << std::endl;
  BTreeIndex index(relationName, intIndexName, bufMgr, offsetof(tuple,i), INTEGER);

	// run some tests
	checkPassFail(intScan(&index,25,GT,40,LT), 14)
	checkPassFail(intScan(&index,20,GTE,35,LTE), 16)
	checkPassFail(intScan(&index,-3,GT,3,LT), 3)
	checkPassFail(intScan(&index,996,GT,1001,LT), 4)
	checkPassFail(intScan(&index,-1000,GT,6000,LT), 5000)
	checkPassFail(intScan(&index,300,GT,400,LT), 99)
	checkPassFail(intScan(&index,3000,GTE,4000,LT), 1000)
}

void intTests_sparse()
{
  std::cout << "Create a B+ Tree index on the integer field" << std::endl;
  BTreeIndex index(relationName, intIndexName, bufMgr, offsetof(tuple,i), INTEGER);

	// run some tests
	checkPassFail(intScan(&index,25,GT,40,LT), 2)
	checkPassFail(intScan(&index,20,GTE,35,LTE), 4)
	checkPassFail(intScan(&index,-3,GT,3,LT), 1)
	checkPassFail(intScan(&index,996,GT,1001,LT), 1)
	checkPassFail(intScan(&index,-1000,GT,30000,LT), 5000)
	checkPassFail(intScan(&index,300,GT,400,LT), 19)
	checkPassFail(intScan(&index,3000,GTE,4000,LT), 200)
}

void intTests_negative()
{
  std::cout << "Create a B+ Tree index on the integer field" << std::endl;
  BTreeIndex index(relationName, intIndexName, bufMgr, offsetof(tuple,i), INTEGER);

	// run some tests
	checkPassFail(intScan(&index,-25,GT,40,LT), 64)
	checkPassFail(intScan(&index,-20,GTE,35,LTE), 56)
	checkPassFail(intScan(&index,-3,GT,3,LT), 5)
	checkPassFail(intScan(&index,996,GT,1001,LT), 4)
	checkPassFail(intScan(&index,-1000,GT,6000,LT), 3499)
	checkPassFail(intScan(&index,300,GT,400,LT), 99)
	checkPassFail(intScan(&index,-4000,GTE,-1000,LT), 1500)
}

void intTests_empty() 
{
	std::cout << "Create a emptry B+ Tree index on the integer field" << std::endl;
	BTreeIndex index(relationName, intIndexName, bufMgr, offsetof(tuple,i), INTEGER);

	// run some tests
	checkPassFail(intScan(&index,25,GT,40,LT), 0)
	checkPassFail(intScan(&index,20,GTE,35,LTE), 0)
	checkPassFail(intScan(&index,-3,GT,3,LT), 0)
	checkPassFail(intScan(&index,996,GT,1001,LT), 0)
	checkPassFail(intScan(&index,0,GT,1,LT), 0)
	checkPassFail(intScan(&index,300,GT,400,LT), 0)
	checkPassFail(intScan(&index,3000,GTE,4000,LT), 0)

}

int intScan(BTreeIndex * index, int lowVal, Operator lowOp, int highVal, Operator highOp)
{
  RecordId scanRid;
	Page *curPage;

  std::cout << "Scan for ";
  if( lowOp == GT ) { std::cout << "("; } else { std::cout << "["; }
  std::cout << lowVal << "," << highVal;
  if( highOp == LT ) { std::cout << ")"; } else { std::cout << "]"; }
  std::cout << std::endl;

  int numResults = 0;
	
	try
	{
	bufMgr->flushFile(file1);
  	index->startScan(&lowVal, lowOp, &highVal, highOp);
	}
	catch(NoSuchKeyFoundException e)
	{
    std::cout << "No Key Found satisfying the scan criteria." << std::endl;
		return 0;
	}

	while(1)
	{
		try
		{
			index->scanNext(scanRid);
			bufMgr->readPage(file1, scanRid.page_number, curPage);
			RECORD myRec = *(reinterpret_cast<const RECORD*>(curPage->getRecord(scanRid).data()));
			bufMgr->unPinPage(file1, scanRid.page_number, false);

			if( numResults < 5 )
			{
				std::cout << "at:" << scanRid.page_number << "," << scanRid.slot_number;
				std::cout << " -->:" << myRec.i << ":" << myRec.d << ":" << myRec.s << ":" <<std::endl;
			}
			else if( numResults == 5 )
			{
				std::cout << "..." << std::endl;
			}
		}
		catch(IndexScanCompletedException e)
		{
			break;
		}

		numResults++;
	}

  if( numResults >= 5 )
  {
    std::cout << "Number of results: " << numResults << std::endl;
  }
  index->endScan();
  std::cout << std::endl;

	return numResults;
}

// -----------------------------------------------------------------------------
// doubleTests
// -----------------------------------------------------------------------------

void doubleTests()
{
  std::cout << "Create a B+ Tree index on the double field" << std::endl;
  BTreeIndex index(relationName, doubleIndexName, bufMgr, offsetof(tuple,d), DOUBLE);

	// run some tests
	checkPassFail(doubleScan(&index,25,GT,40,LT), 14)
	checkPassFail(doubleScan(&index,20,GTE,35,LTE), 16)
	checkPassFail(doubleScan(&index,-3,GT,3,LT), 3)
	checkPassFail(doubleScan(&index,996,GT,1001,LT), 4)
	checkPassFail(doubleScan(&index,0,GT,1,LT), 0)
	checkPassFail(doubleScan(&index,300,GT,400,LT), 99)
	checkPassFail(doubleScan(&index,3000,GTE,4000,LT), 1000)
}

void doubleTests_sparse()
{
	std::cout << "Create a B+ Tree index on the double field" << std::endl;
	BTreeIndex index(relationName, doubleIndexName, bufMgr, offsetof(tuple,d), DOUBLE);

	checkPassFail(doubleScan(&index,25,GT,40,LT), 2)
	checkPassFail(doubleScan(&index,20,GTE,35,LTE), 4)
	checkPassFail(doubleScan(&index,-3,GT,3,LT), 1)
	checkPassFail(doubleScan(&index,996,GT,1001,LT), 1)
	checkPassFail(doubleScan(&index,-1000,GT,30000,LT), 5000)
	checkPassFail(doubleScan(&index,300,GT,400,LT), 19)
	checkPassFail(doubleScan(&index,3000,GTE,4000,LT), 200)
}

void doubleTests_negative()
{
	std::cout << "Create a B+ Tree index on the double field" << std::endl;
	BTreeIndex index(relationName, doubleIndexName, bufMgr, offsetof(tuple,d), DOUBLE);

	checkPassFail(doubleScan(&index,-25,GT,40,LT), 64)
	checkPassFail(doubleScan(&index,-20,GTE,35,LTE), 56)
	checkPassFail(doubleScan(&index,-3,GT,3,LT), 5)
	checkPassFail(doubleScan(&index,996,GT,1001,LT), 4)
	checkPassFail(doubleScan(&index,-1000,GT,6000,LT), 3499)
	checkPassFail(doubleScan(&index,300,GT,400,LT), 99)
	checkPassFail(doubleScan(&index,-4000,GTE,-1000,LT), 1500)
}

void doubleTests_empty()
{
  std::cout << "Create a B+ Tree index on the double field" << std::endl;
  BTreeIndex index(relationName, doubleIndexName, bufMgr, offsetof(tuple,d), DOUBLE);

	// run some tests
	checkPassFail(doubleScan(&index,25,GT,40,LT), 0)
	checkPassFail(doubleScan(&index,20,GTE,35,LTE), 0)
	checkPassFail(doubleScan(&index,-3,GT,3,LT), 0)
	checkPassFail(doubleScan(&index,996,GT,1001,LT), 0)
	checkPassFail(doubleScan(&index,0,GT,1,LT), 0)
	checkPassFail(doubleScan(&index,300,GT,400,LT), 0)
	checkPassFail(doubleScan(&index,3000,GTE,4000,LT), 0)
}

int doubleScan(BTreeIndex * index, double lowVal, Operator lowOp, double highVal, Operator highOp)
{
  RecordId scanRid;
	Page *curPage;

  std::cout << "Scan for ";
  if( lowOp == GT ) { std::cout << "("; } else { std::cout << "["; }
  std::cout << lowVal << "," << highVal;
  if( highOp == LT ) { std::cout << ")"; } else { std::cout << "]"; }
  std::cout << std::endl;

  int numResults = 0;

	try
	{
  	index->startScan(&lowVal, lowOp, &highVal, highOp);
	}
	catch(NoSuchKeyFoundException e)
	{
    std::cout << "No Key Found satisfying the scan criteria." << std::endl;
		return 0;
	}

	while(1)
	{
		try
		{
			index->scanNext(scanRid);
			bufMgr->readPage(file1, scanRid.page_number, curPage);
			RECORD myRec = *(reinterpret_cast<const RECORD*>(curPage->getRecord(scanRid).data()));
			bufMgr->unPinPage(file1, scanRid.page_number, false);

			if( numResults < 5 )
			{
				std::cout << "rid:" << scanRid.page_number << "," << scanRid.slot_number;
				std::cout << " -->:" << myRec.i << ":" << myRec.d << ":" << myRec.s << ":" <<std::endl;
			}
			else if( numResults == 5 )
			{
				std::cout << "..." << std::endl;
			}
		}
		catch(IndexScanCompletedException e)
		{
			break;
		}

		numResults++;
	}

  if( numResults >= 5 )
  {
    std::cout << "Number of results: " << numResults << std::endl;
  }
  index->endScan();
  std::cout << std::endl;

	return numResults;
}

// -----------------------------------------------------------------------------
// stringTests
// -----------------------------------------------------------------------------

void stringTests()
{
  std::cout << "Create a B+ Tree index on the string field" << std::endl;
  BTreeIndex index(relationName, stringIndexName, bufMgr, offsetof(tuple,s), STRING);

	// run some tests
	checkPassFail(stringScan(&index,25,GT,40,LT), 14)
	checkPassFail(stringScan(&index,20,GTE,35,LTE), 16)
	checkPassFail(stringScan(&index,-3,GT,3,LT), 3)
	checkPassFail(stringScan(&index,996,GT,1001,LT), 4)
	checkPassFail(stringScan(&index,0,GT,1,LT), 0)
	checkPassFail(stringScan(&index,300,GT,400,LT), 99)
	checkPassFail(stringScan(&index,3000,GTE,4000,LT), 1000)
}

void stringTests_sparse()
{
  std::cout << "Create a B+ Tree index on the string field" << std::endl;
  BTreeIndex index(relationName, stringIndexName, bufMgr, offsetof(tuple,s), STRING);

	// run some tests
	checkPassFail(stringScan(&index,25,GT,40,LT), 2)
	checkPassFail(stringScan(&index,20,GTE,35,LTE), 4)
	checkPassFail(stringScan(&index,-3,GT,3,LT), 1)
	checkPassFail(stringScan(&index,996,GT,1001,LT), 1)
	checkPassFail(stringScan(&index,0,GT,1,LT), 0)
	checkPassFail(stringScan(&index,300,GT,400,LT), 19)
	checkPassFail(stringScan(&index,3000,GTE,4000,LT), 200)
}

void stringTests_negative()
{
  std::cout << "Create a B+ Tree index on the string field" << std::endl;
  BTreeIndex index(relationName, stringIndexName, bufMgr, offsetof(tuple,s), STRING);

	// run some tests
	checkPassFail(stringScan(&index,0,GT,40,LT), 39)
	checkPassFail(stringScan(&index,-2000,GTE,35,LTE), 537)
	checkPassFail(stringScan(&index,-2497,GT,3,LT), 6)
	checkPassFail(stringScan(&index,996,GT,1001,LT), 4)
	checkPassFail(stringScan(&index,-1000,GT,6000,LT), 4000)
	checkPassFail(stringScan(&index,300,GT,400,LT), 99)
	checkPassFail(stringScan(&index,-1000,GTE,-4000,LT), 1501)

}

void stringTests_empty()
{
  std::cout << "Create a B+ Tree index on the string field" << std::endl;
  BTreeIndex index(relationName, stringIndexName, bufMgr, offsetof(tuple,s), STRING);

	// run some tests
	checkPassFail(stringScan(&index,25,GT,40,LT), 0)
	checkPassFail(stringScan(&index,20,GTE,35,LTE), 0)
	checkPassFail(stringScan(&index,-3,GT,3,LT), 0)
	checkPassFail(stringScan(&index,996,GT,1001,LT), 0)
	checkPassFail(stringScan(&index,0,GT,1,LT), 0)
	checkPassFail(stringScan(&index,300,GT,400,LT), 0)
	checkPassFail(stringScan(&index,3000,GTE,4000,LT), 0)
}

int stringScan(BTreeIndex * index, int lowVal, Operator lowOp, int highVal, Operator highOp)
{
  RecordId scanRid;
	Page *curPage;

  std::cout << "Scan for ";
  if( lowOp == GT ) { std::cout << "("; } else { std::cout << "["; }
  std::cout << lowVal << "," << highVal;
  if( highOp == LT ) { std::cout << ")"; } else { std::cout << "]"; }
  std::cout << std::endl;

  char lowValStr[100];
  sprintf(lowValStr,"%05d string record",lowVal);
  char highValStr[100];
  sprintf(highValStr,"%05d string record",highVal);

  int numResults = 0;

	try
	{
  	index->startScan(lowValStr, lowOp, highValStr, highOp);
	}
	catch(NoSuchKeyFoundException e)
	{
    std::cout << "No Key Found satisfying the scan criteria." << std::endl;
		return 0;
	}

	while(1)
	{
		try
		{
			index->scanNext(scanRid);
			bufMgr->readPage(file1, scanRid.page_number, curPage);
			RECORD myRec = *(reinterpret_cast<const RECORD*>(curPage->getRecord(scanRid).data()));
			bufMgr->unPinPage(file1, scanRid.page_number, false);

			if( numResults < 5 )
			{
				std::cout << "rid:" << scanRid.page_number << "," << scanRid.slot_number;
				std::cout << " -->:" << myRec.i << ":" << myRec.d << ":" << myRec.s << ":" <<std::endl;
			}
			else if( numResults == 5 )
			{
				std::cout << "..." << std::endl;
			}
		}
		catch(IndexScanCompletedException e)
		{
			break;
		}

		numResults++;
	}

  if( numResults >= 5 )
  {
    std::cout << "Number of results: " << numResults << std::endl;
  }
  index->endScan();
  std::cout << std::endl;

	return numResults;
}

// -----------------------------------------------------------------------------
// errorTests
// -----------------------------------------------------------------------------

void errorTests()
{
	std::cout << "Error handling tests" << std::endl;
	std::cout << "--------------------" << std::endl;
	// Given error test

	try
	{
		File::remove(relationName);
	}
	catch(FileNotFoundException e)
	{
	}

  file1 = new PageFile(relationName, true);
	
  // initialize all of record1.s to keep purify happy
  memset(record1.s, ' ', sizeof(record1.s));
	PageId new_page_number;
  Page new_page = file1->allocatePage(new_page_number);

  // Insert a bunch of tuples into the relation.
	for(int i = 0; i <10; i++ ) 
	{
    sprintf(record1.s, "%05d string record", i);
    record1.i = i;
    record1.d = (double)i;
    std::string new_data(reinterpret_cast<char*>(&record1), sizeof(record1));

		while(1)
		{
			try
			{
    		new_page.insertRecord(new_data);
				break;
			}
			catch(InsufficientSpaceException e)
			{
				file1->writePage(new_page_number, new_page);
  			new_page = file1->allocatePage(new_page_number);
			}
		}
  }

	file1->writePage(new_page_number, new_page);

  BTreeIndex index(relationName, intIndexName, bufMgr, offsetof(tuple,i), INTEGER);
	
	int int2 = 2;
	int int5 = 5;

	// Scan Tests
	std::cout << "Call endScan before startScan" << std::endl;
	try
	{
		index.endScan();
		std::cout << "ScanNotInitialized Test 1 Failed." << std::endl;
	}
	catch(ScanNotInitializedException e)
	{
		std::cout << "ScanNotInitialized Test 1 Passed." << std::endl;
	}
	
	std::cout << "Call scanNext before startScan" << std::endl;
	try
	{
		RecordId foo;
		index.scanNext(foo);
		std::cout << "ScanNotInitialized Test 2 Failed." << std::endl;
	}
	catch(ScanNotInitializedException e)
	{
		std::cout << "ScanNotInitialized Test 2 Passed." << std::endl;
	}
	
	std::cout << "Scan with bad lowOp" << std::endl;
	try
	{
  	index.startScan(&int2, LTE, &int5, LTE);
		std::cout << "BadOpcodesException Test 1 Failed." << std::endl;
	}
	catch(BadOpcodesException e)
	{
		std::cout << "BadOpcodesException Test 1 Passed." << std::endl;
	}
	
	std::cout << "Scan with bad highOp" << std::endl;
	try
	{
  	index.startScan(&int2, GTE, &int5, GTE);
		std::cout << "BadOpcodesException Test 2 Failed." << std::endl;
	}
	catch(BadOpcodesException e)
	{
		std::cout << "BadOpcodesException Test 2 Passed." << std::endl;
	}


	std::cout << "Scan with bad range" << std::endl;
	try
	{
  	index.startScan(&int5, GTE, &int2, LTE);
		std::cout << "BadScanrangeException Test 1 Failed." << std::endl;
	}
	catch(BadScanrangeException e)
	{
		std::cout << "BadScanrangeException Test 1 Passed." << std::endl;
	}

	deleteRelation();
}

void deleteRelation()
{
	if(file1)
	{
		bufMgr->flushFile(file1);
		delete file1;
		file1 = NULL;
	}
	try
	{
		File::remove(relationName);
	}
	catch(FileNotFoundException e)
	{
	}
}
