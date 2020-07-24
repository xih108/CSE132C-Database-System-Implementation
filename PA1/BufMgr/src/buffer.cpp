/**
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include <memory>
#include <iostream>
#include "buffer.h"
#include "exceptions/buffer_exceeded_exception.h"
#include "exceptions/page_not_pinned_exception.h"
#include "exceptions/page_pinned_exception.h"
#include "exceptions/bad_buffer_exception.h"
#include "exceptions/hash_not_found_exception.h"

namespace badgerdb { 

BufMgr::BufMgr(std::uint32_t bufs)
	: numBufs(bufs) {
	bufDescTable = new BufDesc[bufs];

  for (FrameId i = 0; i < bufs; i++) 
  {
  	bufDescTable[i].frameNo = i;
  	bufDescTable[i].valid = false;
  }

  bufPool = new Page[bufs];

	int htsize = ((((int) (bufs * 1.2))*2)/2)+1;
  hashTable = new BufHashTbl (htsize);  // allocate the buffer hash table

  clockHand = bufs - 1;
}


BufMgr::~BufMgr() {
	
	for (FrameId i = 0; i < numBufs; i++) 
	{
		BufDesc* temp = &bufDescTable[i];
		// check valid and dirty bit
		if (temp->valid && temp->dirty)
		{
			temp->file->writePage(bufPool[i]);
		}
	}

	// deallocate buffer pool, BufDesc table
	delete[] bufDescTable;
	delete[] bufPool;

	// deallocate hashtable
	delete hashTable;
}

void BufMgr::advanceClock()
{ 
	// advance clockHand clockwise by one
	clockHand  = (clockHand + 1) % numBufs;	
}

void BufMgr::allocBuf(FrameId & frame) 
{
	for (FrameId i = 0; i < 2 * numBufs; i++){
		advanceClock();

		// check for valid 
		if (!bufDescTable[clockHand].valid)
		{
			frame = clockHand;
			bufDescTable[clockHand].Clear();
			return;
		}

		// check for refbit
		if (!bufDescTable[clockHand].refbit)
		{
			// check for page pinned
			if (bufDescTable[clockHand].pinCnt == 0)
			{
				frame = clockHand;
				// writing a dirty page back to disk
				if (bufDescTable[clockHand].dirty){
					bufDescTable[clockHand].file->writePage(bufPool[clockHand]);
				}
				// remove the appropriate entry from the hash table
				hashTable->remove(bufDescTable[clockHand].file, bufDescTable[clockHand].pageNo);
				bufDescTable[clockHand].Clear();

				return;
			}
		}
		else
		{
			// clear refbit
			bufDescTable[clockHand].refbit = false;
		}
		
	}
	// all buffer frames are pinned
	throw BufferExceededException();
}

	
void BufMgr::readPage(File* file, const PageId pageNo, Page*& page)
{
	FrameId frameNo;
	try
	{	
		// check whether the page is already in the buffer pool 
		hashTable->lookup(file, pageNo, frameNo);
		// set refbit
		bufDescTable[frameNo].refbit = true;
		// increment pinCnt
		bufDescTable[frameNo].pinCnt += 1;
		page = &bufPool[frameNo];

	}
	catch(HashNotFoundException e)
	{
		// allocate a buffer frame
		allocBuf(frameNo);
		// read page from disk into buffer pool
		bufPool[frameNo] = file->readPage(pageNo);

		// insert page into hashtable
		hashTable->insert(file, pageNo, frameNo);
		// set up the frame 
		bufDescTable[frameNo].Set(file, pageNo);
		page = &bufPool[frameNo];
	}

	

}


void BufMgr::unPinPage(File* file, const PageId pageNo, const bool dirty) 
{
	FrameId frameNo;
	try 
	{
		// check whether the page is in the buffer pool 
		hashTable->lookup(file, pageNo, frameNo);
	}
	catch(HashNotFoundException e){
		return;
	}
	
	// set dirty bit if dirty is true
	if (dirty)
	{
		bufDescTable[frameNo].dirty = dirty;
	}

	// check the pinCnt of the frame
	if (bufDescTable[frameNo].pinCnt == 0)
	{
		throw PageNotPinnedException(file->filename(), pageNo, frameNo);
	}
	else
	{	
		// decrements the pinCnt of the frame
		bufDescTable[frameNo].pinCnt --;
	}

}

void BufMgr::flushFile(const File* file) 
{
	// scan the bufTable
	for (FrameId i = 0; i < numBufs; i++)
	{
		BufDesc* temp = &bufDescTable[i];
		// check if the page belongs to the file
		
		if (temp->file == file)
		{
			
			// check if the page is valid
			if (!temp->valid)
			{
				throw BadBufferException(i, temp->dirty, temp->valid, temp->refbit);	
			}

			// check if the page is pinned
			if (temp->pinCnt > 0)
			{
				throw PagePinnedException(file->filename(), temp->pageNo, i);
			}

			// check the dirty bit
			if (temp->dirty)
			{
				// flush the dirty page to disk 
				temp->file->writePage(bufPool[i]);

				//set dirty bit to false
				temp->dirty = false;
			}
			// remove the page from the hashtable 
			hashTable->remove(file, temp->pageNo);
			// Clear the page frame
			temp->Clear();
		
		}

	}
}

void BufMgr::allocPage(File* file, PageId &pageNo, Page*& page) 
{
	// allocate an empty page in the file
	Page newpage = file->allocatePage();
	pageNo = newpage.page_number();
	// allocate buf, insert into hash table, set up the frame
	readPage(file, pageNo, page);
	
}

void BufMgr::disposePage(File* file, const PageId PageNo)
{
	FrameId frameNo;
	try
	{
		// check whether the page is already in the buffer pool 
		hashTable->lookup(file, PageNo, frameNo);
		// free frame
		bufDescTable[frameNo].Clear(); 
		// remove the entry from the hash table
		hashTable->remove(file, PageNo);
		// delete the page from file
		file->deletePage(PageNo);
	}
	catch (HashNotFoundException e)
	{
		// delete the page from file
		file->deletePage(PageNo);
	}
	
}

void BufMgr::printSelf(void) 
{
  BufDesc* tmpbuf;
	int validFrames = 0;
  
  for (std::uint32_t i = 0; i < numBufs; i++)
	{
  	tmpbuf = &(bufDescTable[i]);
		std::cout << "FrameNo:" << i << " ";
		tmpbuf->Print();

  	if (tmpbuf->valid == true)
    	validFrames++;
  }

	std::cout << "Total Number of Valid Frames:" << validFrames << "\n";
}

}
