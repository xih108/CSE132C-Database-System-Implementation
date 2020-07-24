/**
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

/**
 * Xinyi He A13561164
 * Zhixian Chen A13680919
 * @brief The purpose of this file is to build a b+ index tree which could perform insertion and range 
 * search opertaion. The tree node will be stored into a file as a page and we could use buffer manager 
 * to get the page from the file.
 */

#include "btree.h"
#include "filescan.h"
#include "exceptions/bad_index_info_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/end_of_file_exception.h"


//#define DEBUG
using namespace std;
namespace badgerdb
{

// -----------------------------------------------------------------------------
// BTreeIndex::BTreeIndex -- Constructor
// -----------------------------------------------------------------------------
/**
 * The constructor of the BTreeIndex.
 * Build the btree into a index file based on the relation file.
 * If the index file exists, open the file. Othewise, create a new file
*/
BTreeIndex::BTreeIndex(const std::string & relationName,
		std::string & outIndexName,
		BufMgr *bufMgrIn,
		const int attrByteOffset,
		const Datatype attrType)
{
	//build the index file
	ostringstream idxStr;
	idxStr << relationName << '.' << attrByteOffset;
	string indexName = idxStr.str();

	outIndexName = indexName;

	//initialize the field for the btree index
	this->attributeType = attrType;
	this->attrByteOffset = attrByteOffset;
	this->bufMgr = bufMgrIn;
	this->scanExecuting = false;

	if (this->attributeType == INTEGER){
		this->nodeOccupancy = INTARRAYNONLEAFSIZE;
		this->leafOccupancy = INTARRAYLEAFSIZE;
	}
	else if (this->attributeType == DOUBLE){
		this->nodeOccupancy = DOUBLEARRAYNONLEAFSIZE;
		this->leafOccupancy = DOUBLEARRAYLEAFSIZE;
	}
	else {
		this->nodeOccupancy = STRINGARRAYNONLEAFSIZE;
		this->leafOccupancy = STRINGARRAYLEAFSIZE;
	}
	
	//if the index file exists
	if (File::exists(indexName)){
		this->file = new BlobFile(indexName, false);
		this->headerPageNum = this->file->getFirstPageNo();
		Page* headerPage;
		this->bufMgr->readPage(this->file, headerPageNum, headerPage);
		IndexMetaInfo* metaInfo = (IndexMetaInfo*) headerPage;
		this->bufMgr->unPinPage(this->file, this->headerPageNum, false); 

		//check whether the meta page matches
		if (strcmp(metaInfo->relationName, relationName.c_str()) != 0 
		|| metaInfo->attrByteOffset != attrByteOffset
		|| metaInfo->attrType != attrType)
			throw BadIndexInfoException("Relation does not match");

		this->rootPageNum = metaInfo->rootPageNo;
	}
	//create a new index file to build the btree
	else {
		this->file = new BlobFile(outIndexName, true);

		//allocate the header page and root page
		Page* headerPage;
		this->bufMgr->allocPage(this->file, this->headerPageNum, headerPage);
		Page* rootPage;
		this->bufMgr->allocPage(this->file, this->rootPageNum, rootPage);
		
		//initalize the meta info page, header page
		FileScan fscan(relationName, this->bufMgr);
		IndexMetaInfo* metaInfo = (IndexMetaInfo*) headerPage;
		strcpy(metaInfo->relationName, relationName.c_str());
		metaInfo->attrByteOffset = attrByteOffset;
		metaInfo->attrType = attrType;
		metaInfo->rootPageNo = this->rootPageNum;

		//create a dummy root with two leaf nodes so that the root cannot be nonleaf
		Page* leaf1;
		PageId leafNo1;
		this->bufMgr->allocPage(this->file, leafNo1, leaf1);

		Page* leaf2;
		PageId leafNo2;
		this->bufMgr->allocPage(this->file, leafNo2, leaf2);
		
		//for the int type, insert dummy root node with key 1000
		if (this->attributeType == INTEGER ) {
			NonLeafNodeInt* rootNode = (NonLeafNodeInt*) rootPage;
			rootNode->level = 1;

			rootNode->keyArray[0] = 1000;

			rootNode->pageNoArray[0] = leafNo1;
			rootNode->pageNoArray[1] = leafNo2;
			
			LeafNodeInt* leafNode1 = (LeafNodeInt*) leaf1;
			LeafNodeInt* leafNode2 = (LeafNodeInt*) leaf2;
				
			leafNode1->rightSibPageNo = leafNo2;
		}
		//for the double type, insert dummy root node with key 1000
		else if (this->attributeType == DOUBLE) {
			NonLeafNodeDouble* rootNode = (NonLeafNodeDouble*) rootPage;
			rootNode->level = 1;

			rootNode->keyArray[0] = 1000;

			rootNode->pageNoArray[0] = leafNo1;
			rootNode->pageNoArray[1] = leafNo2;
			
			LeafNodeDouble* leafNode1 = (LeafNodeDouble*) leaf1;
			LeafNodeDouble* leafNode2 = (LeafNodeDouble*) leaf2;
				
			leafNode1->rightSibPageNo = leafNo2;
		}
		//for the string type, insert dummy root node with key 011111111
		else {
			NonLeafNodeString* rootNode = (NonLeafNodeString*) rootPage;
			rootNode->level = 1;

			std::string dummykey = "011111111\0";
			strncpy(rootNode->keyArray[0], dummykey.c_str(), 10);

			rootNode->pageNoArray[0] = leafNo1;
			rootNode->pageNoArray[1] = leafNo2;
			
			LeafNodeString* leafNode1 = (LeafNodeString*) leaf1;
			LeafNodeString* leafNode2 = (LeafNodeString*) leaf2;
				
			leafNode1->rightSibPageNo = leafNo2;
		}
		
		this->bufMgr->unPinPage(this->file, leafNo1, true); 
		this->bufMgr->unPinPage(this->file, leafNo2, true); 

		this->bufMgr->unPinPage(this->file, this->headerPageNum, true); 
		this->bufMgr->unPinPage(this->file, this->rootPageNum, true); 

		//build the btree based on the relational file using file scan
		try
		{
			RecordId scanRid;
			while(1)
			{
				fscan.scanNext(scanRid);
				std::string recordStr = fscan.getRecord();
				const char *record = recordStr.c_str();
				void* key = (void *)(record + attrByteOffset);
				
				//insert the key into the tree
				insertEntry(key, scanRid);
			}
		}
		catch(EndOfFileException e)
		{
			// std::cout << "Read all records" << std::endl;
		}

	}

}

// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------
/**
 * The destructor of b+ index tree.
 * Finish the scan operation. Flush the file and delete the index file
*/
BTreeIndex::~BTreeIndex()
{
	//clear up state variable
	if (this->scanExecuting)
		this->endScan();

	//flush the file
	this->bufMgr->flushFile(this->file);
	delete this->file;
}

template <class T, class nonleaf_T>
const void BTreeIndex::updateRoot(T midpagekey){
	
	//allocate page for the new root page
	Page* newroot;
	PageId newrootNo;
	this->bufMgr->allocPage(this->file, newrootNo, newroot);
	nonleaf_T* newrootnode = (nonleaf_T*) newroot;
	
	//insert key into the new root node and link the children to the new root
	assign(newrootnode->keyArray[0], midpagekey.key);
	newrootnode->pageNoArray[0] = this->rootPageNum;
	newrootnode->pageNoArray[1] = midpagekey.pageNo;
	newrootnode->level = 0;
	
	//update the rootpagenum fielf
	this->rootPageNum = newrootNo;

	//update the meta info for the root page number
	Page* headerPage;
	this->bufMgr->readPage(this->file, this->headerPageNum, headerPage);
	IndexMetaInfo* metaInfo = (IndexMetaInfo*) headerPage;
	metaInfo->rootPageNo = this->rootPageNum;
	this->bufMgr->unPinPage(this->file, this->headerPageNum, true); 
	this->bufMgr->unPinPage(this->file, newrootNo, true);
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------
/**
 * Insert a new entry into the tree. Start from the root node to find the place to insert the entry
 * into the leaf node. The value returned from the insertRecursive indicates whether the root node 
 * needs to be updated
 * @param key			Key to insert, pointer to integer/double/ string
 * @param rid			Record ID of a record whose entry is getting inserted into the index.
**/
const void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{	

	//insert the key into the tree based on the attribute type
	if (this->attributeType == INTEGER){
		RIDKeyPair<int> ridkey;
		ridkey.set(rid, *(int *)key);
		PageKeyPair<int> midpagekey = this->insertRecursive<int, NonLeafNodeInt, LeafNodeInt>(this->rootPageNum, ridkey);
		
		//check whether the root node needs to be updated
		if (midpagekey.pageNo != 0){
			updateRoot<PageKeyPair<int>, NonLeafNodeInt>(midpagekey);
		}
	}
	else if (this->attributeType == DOUBLE){
		RIDKeyPair<double> ridkey;
		ridkey.set(rid, *(double *)key);
		PageKeyPair<double> midpagekey = this->insertRecursive<double, NonLeafNodeDouble, LeafNodeDouble>(this->rootPageNum, ridkey);
		
		//check whether the root node needs to be updated
		if (midpagekey.pageNo != 0){
			updateRoot<PageKeyPair<double>, NonLeafNodeDouble>(midpagekey);
		}
	}
	else {
		
		//get the first 10 character of the inserted key
		RIDKeyPair<string> ridkey;
		char a[10];
		strncpy(a, (char *)key, 9);
		a[9] = '\0';
		string str = a;
		ridkey.set(rid, str);
		
		PageKeyPair<string> midpagekey = this->insertRecursive<string, NonLeafNodeString, LeafNodeString>(this->rootPageNum, ridkey);
		
		//check whehter the root node needs to be updated
		if (midpagekey.pageNo != 0){
			updateRoot<PageKeyPair<string>, NonLeafNodeString>(midpagekey);
		}
	} 	
}

/**
 * The assign function to give one double to another double
*/
const void BTreeIndex::assign(double& left, double right) {
	left = right;
}

/**
 * The assign function to give one int to another int
*/
const void BTreeIndex::assign(int& left, int right) {
	left = right;
}

/**
 * The assign function to give one char array to another char array
*/
const void BTreeIndex::assign(char left[], char right[]) {
	strncpy(left, right, 10);
}

/**
 * The assign function to give one string to another string
*/
const void BTreeIndex::assign(string& left, string right) {
	left = right;
}

/**
 * The assign function to give one string to another char array
*/
const void BTreeIndex::assign(char left[], string right) {
	strncpy(left, right.c_str(), 10);
}

/**
 * Traverse the tree from the current node to find a place insert the entry
 * if the current node's child is leafnode, we insert the entry
 * If the current node's child is nonleafnode, we recursively call insertrecursive to find the place
 * The return pagekeypair indicate that if the leaf node got split, then mid key to insert into the parent
 * @param root			current node page id 
 * @param ridkey		the inserted entry's ridkey pair 
**/
template <class T, class nonleaf_T, class leaf_T>
PageKeyPair<T> BTreeIndex::insertRecursive(PageId root,  RIDKeyPair<T> ridkey){

	Page* page;
	this->bufMgr->readPage(this->file, root, page);
	nonleaf_T* rootNode = (nonleaf_T*) page;
	this->bufMgr->unPinPage(this->file, root, true);

	//if the child of the current node is nonleaf node 
	if (rootNode -> level == 0){

		//traverse down to find the leafnode
		PageId nextNode = this->findPageNoInNonLeaf<T, nonleaf_T>(root, ridkey);
		PageKeyPair<T> midpagekey;
		midpagekey = this->insertRecursive<T, nonleaf_T, leaf_T>(nextNode, ridkey);

		//check whether the child node got split, if not do nothing
		if (midpagekey.pageNo == 0){
			return midpagekey;
		}
		else 
		{
			// if the current page node is not full, insert the key into the current node
			if (rootNode -> pageNoArray[this->nodeOccupancy] == 0){
				insertMidintoParent<T, nonleaf_T>(root, midpagekey);
				midpagekey.set(0, midpagekey.key);
				return midpagekey;
			}
			//if the node is full, split nonleaf node
			else {
				PageKeyPair<T> leftpagekey;
				leftpagekey.set(root, midpagekey.key);
				return splitNonLeafNode<T, nonleaf_T>(root, leftpagekey, midpagekey);	
			}
			
		}
	}
	//the child of the current node is leaf node, then insert
	else
	{
		//get the leaf node
		PageId leafNo = this->findPageNoInNonLeaf<T, nonleaf_T>(root, ridkey);
		Page* leafpage ;
		this->bufMgr->readPage(this->file, leafNo, leafpage);
		leaf_T* leafNode = (leaf_T*) leafpage;
		this->bufMgr->unPinPage(this->file, leafNo, true); 
		T findkey = ridkey.key;
		
		//if the node is not full, insert directly
		if (leafNode->ridArray[this->leafOccupancy-1].page_number == 0){
			bool insert = false;

			//loop the array to find the right place to insert the key	
			for (int i = this->leafOccupancy-1; i >= 0; i--){
				if (leafNode->ridArray[i].page_number != 0){
					assign(leafNode->keyArray[i+1], leafNode->keyArray[i]);
					leafNode->ridArray[i+1] = leafNode->ridArray[i];
					if (findkey >= leafNode-> keyArray[i]){
						assign(leafNode->keyArray[i+1], findkey);
						leafNode->ridArray[i+1] = ridkey.rid;
						insert = true;
						break;
					}
				}
			}

			//insert at the end 
			if (!insert){
				assign(leafNode->keyArray[0], findkey);
				leafNode->ridArray[0] = ridkey.rid;
			}	

			//no need to insert key to parent node
			PageKeyPair<T> midpagekey;
			midpagekey.set(0, findkey);
			return midpagekey;
		}
		//if the child node is full, split the child node
		else
		{
			PageKeyPair<T> midpagekey = splitLeafNode<T, leaf_T>(leafNo, ridkey);

			//if the parent node is not full, insert the mid key into the parent node directly
			if (rootNode -> pageNoArray[this->nodeOccupancy] == 0){
				insertMidintoParent<T, nonleaf_T>(root, midpagekey);
				midpagekey.set(0, midpagekey.key);
				return midpagekey;
			}
			//if the parent node is full, split the parent node to insert
			else {
				PageKeyPair<T> leftpagekey;
				leftpagekey.set(leafNo, midpagekey.key);
				return splitNonLeafNode<T, nonleaf_T>(root, leftpagekey, midpagekey);	
			}
		}

	}
}


/**
 * Insert the mid key to the parent node and parent node is not full
 * Inser the mid key to the right position and shift the original value accrordingly
 * @param rootNo		Parent node in order to insert the entry into
 * @param midpagekey	The entry to be inserted
**/
template <class T, class nonleaf_T>
const void BTreeIndex:: insertMidintoParent(PageId rootNo, PageKeyPair<T> midpagekey){
	Page* rootpage;
	this->bufMgr->readPage(this->file, rootNo, rootpage);
	nonleaf_T* rootNode = (nonleaf_T*) rootpage;
	this->bufMgr->unPinPage(this->file, rootNo, true); 
	bool insert = false;
	
	//insert the mid key into parent node directly from and shift the original keys
	for (int i = this->nodeOccupancy-1; i >= 0; i--){
		if (rootNode->pageNoArray[i+1] != 0){
			assign(rootNode->keyArray[i+1], rootNode->keyArray[i]);
			rootNode->pageNoArray[i+2] = rootNode->pageNoArray[i+1];

			//find the place to insert
			if (midpagekey.key > rootNode->keyArray[i]){
				assign(rootNode->keyArray[i+1], midpagekey.key);
				rootNode->pageNoArray[i+2] = midpagekey.pageNo;
				insert = true;
				break;
			}
		}
	}
	
	//insert at the front if the key is smallest
	if (!insert){
		assign(rootNode->keyArray[0], midpagekey.key);
		rootNode->pageNoArray[1] = midpagekey.pageNo;
	}
}


/**
 * The function to split the leafnode if the leaf node is already full
 * Allocate a new page to split the leafnode and assign the original entry based on key value
 * The return value of this function is the midkey entry of the leaf node (include the inserted key).
 * The return value is needed to insert to the parent node
 * @param leafNo		The leafnode which is needed to insert the entry and is full
 * @param ridkey		The ridkey pair to be inserted
**/
template <class T, class leaf_T>
PageKeyPair<T> BTreeIndex::splitLeafNode(PageId leafNo, RIDKeyPair<T> ridkey){

	//get the leaf node
	Page* leafPage; 
	this->bufMgr->readPage(this->file, leafNo, leafPage);
	leaf_T* leafNode = (leaf_T *) leafPage;

	//allocate new page for the splited page
	Page* leaf2;
	PageId leafNo2;
	this->bufMgr->allocPage(this->file, leafNo2, leaf2);
	leaf_T* leafNode2 = (leaf_T *) leaf2;

	//store key and original keys in the leaf node into a tmparray sorted
	T tmpkeyArray[ this->leafOccupancy + 1];
	RecordId tmpridArray[ this->leafOccupancy + 1 ];
	int index = 0;
	bool insert = false;

	//find the right place to insert the inserted record
	for (int i = 0; i < this->leafOccupancy ; i++){
		if (ridkey.key < leafNode->keyArray[i] && !insert){
			assign(tmpkeyArray[index], ridkey.key);
			tmpridArray[index] = ridkey.rid;
			insert = true;
			index ++;	
		}
		
		assign(tmpkeyArray[index], leafNode->keyArray[i]);
		tmpridArray[index] = leafNode->ridArray[i];
		index ++;
		
	}

	//insert at the end if the inserted key is largest
	if (!insert){
		assign(tmpkeyArray[index], ridkey.key);
		tmpridArray[index] = ridkey.rid;
	}
	
	//find the mid key to insert into the parent node
	int mid = (this->leafOccupancy + 1) / 2;
	T key = tmpkeyArray[mid];

	for (int i = 0; i < this->leafOccupancy + 1; i++){

		//assign left half of the tmpArray to original node (left node)
		if (i < mid){

			assign(leafNode->keyArray[i], tmpkeyArray[i]);
			leafNode->ridArray[i] = tmpridArray[i];
		}
		//assign right half of the tmparray to right node (allocated node)
		if (i >= mid){
			assign(leafNode2->keyArray[i - mid], tmpkeyArray[i]);
			leafNode2->ridArray[i - mid] = tmpridArray[i];
		}
		//set the rest of original node to empty
		if (i >= mid && i < this->leafOccupancy ){
			leafNode->ridArray[i].page_number = 0;
		}
	}

	//connect the right node and left node
	leafNode2->rightSibPageNo = leafNode->rightSibPageNo;
	leafNode->rightSibPageNo = leafNo2;
	this->bufMgr->unPinPage(this->file, leafNo, true); 
	this->bufMgr->unPinPage(this->file, leafNo2, true); 

	//return the mid key in order to insert to parent node
	PageKeyPair<T> midpagekey;
	midpagekey.set(leafNo2, key);
	return midpagekey;
}

/**
 * The function to split the nonleafnode if the nonleaf node is already full
 * Same as the splitleafNode, we sort the original entries with the inserted entry
 * Split these into two nodes. The return value is needed to insert the midkey into the parent node
 * @param nonLeafNo		The pageid of the nonleaf node to insert entry
 * @param leftpagekey	Contains the key value and the midkey's left pagekey pair
 * @param rightpagekey	The midkey's right pagekey pair
**/
template <class T, class nonleaf_T>
PageKeyPair<T> BTreeIndex::splitNonLeafNode(PageId nonLeafNo, PageKeyPair<T> leftpagekey, PageKeyPair<T> rightpagekey){

	//get the full nonleaf node and allocate a new page as the node at its right at the same level
	Page* nonLeafPage;
	this->bufMgr->readPage(this->file, nonLeafNo, nonLeafPage);
	nonleaf_T* nonleafNode = (nonleaf_T*) nonLeafPage;

	Page* nonleaf2;
	PageId nonleafNo2;
	this->bufMgr->allocPage(this->file, nonleafNo2, nonleaf2);
	nonleaf_T* nonleafNode2 = (nonleaf_T*) nonleaf2;
	nonleafNode2->level = nonleafNode->level;

	//create a tmp key to sort the inserted key from child and full nonleaf node
	T tmpkeyArray[ this->nodeOccupancy + 1];
	PageId tmppageNoArray[ this->nodeOccupancy + 2 ];
	int index = 0;
	int index2 = 0;
	bool insert = false;

	//insert all key and pageNo into the tmp array we create
	for (int i = 0; i < this->nodeOccupancy ; i++){
		if (leftpagekey.key < nonleafNode->keyArray[i] && !insert){
			assign(tmpkeyArray[index], leftpagekey.key);
			tmppageNoArray[index2] = leftpagekey.pageNo;
			index2 = index + 1;
			tmppageNoArray[index2] = rightpagekey.pageNo;
			index ++;	
			index2 ++; 
			insert = true;
		}
		
		//insert the original value if inserted key is not inserted
		if (!insert){
			assign(tmpkeyArray[index], nonleafNode->keyArray[i]);
			tmppageNoArray[index2] = nonleafNode->pageNoArray[i];
			index ++;	
			index2 ++;
		}
		//insert the original value if inserted key is inserted
		if (insert){
			assign(tmpkeyArray[index], nonleafNode->keyArray[i]);
			tmppageNoArray[index2] = nonleafNode->pageNoArray[i+1];
			index ++;	
			index2 ++;
		}
		
	}
	//if key is the largest among the original keys in nonleaf node 
	if (!insert){
		assign(tmpkeyArray[index], leftpagekey.key);
		index2 = index + 1;
		tmppageNoArray[index2] = rightpagekey.pageNo;
	}

	//find the mid point to split the tmp array
	int mid = (this->nodeOccupancy + 1) / 2;
	T key = tmpkeyArray[mid];

	for (int i = 0; i < this->nodeOccupancy + 1; i++){

		//assign the left half to the original node
		if (i < mid){
			assign(nonleafNode->keyArray[i], tmpkeyArray[i]);
			nonleafNode->pageNoArray[i] = tmppageNoArray[i];
		}
		//insert the mid key into the allocated node
		if (i == mid){
			nonleafNode->pageNoArray[i] = tmppageNoArray[i];
			nonleafNode2->pageNoArray[0] = tmppageNoArray[i+1];
		}
		//insert the right half to the allocated node
		if (i > mid){
			assign(nonleafNode2->keyArray[i - mid - 1], tmpkeyArray[i]);
			nonleafNode2->pageNoArray[i - mid] = tmppageNoArray[i];
		}
		//clear the original node rest half
		if (i > mid && i <= this->nodeOccupancy ){
			nonleafNode->pageNoArray[i] = 0;
		}
	}
	
	this->bufMgr->unPinPage(this->file, nonLeafNo, true); 
	this->bufMgr->unPinPage(this->file, nonleafNo2, true); 
	
	//return the allocated pagekey pair
	PageKeyPair<T> pagekey;
	pagekey.set(nonleafNo2, key);
	return pagekey;
	
}

/**
 * Find the next pageNo that might contains the node to insert the entry in order to traverse the tree
 * @param currentNo		The pageId we need to loop to find the next node
 * @param ridkey		The ridkey pair to be inserted
**/
template<class T, class nonleaf_T>
PageId BTreeIndex::findPageNoInNonLeaf(PageId currentNode, RIDKeyPair<T> ridkey ){
	
	//get the current node based on page number
	Page* currentPage;
	this->bufMgr->readPage(this->file, currentNode, currentPage);

	PageId nextNode;
	nonleaf_T*	curnode = (nonleaf_T*) currentPage;
	this->bufMgr->unPinPage(this->file, currentNode, false); 
	
	T findkey = ridkey.key;
	bool find = false;	

	//loop through the current node	
	for (int i = this->nodeOccupancy-1; i >= 0; i--){
		
		if (findkey > curnode->keyArray[i] && curnode->pageNoArray[i+1] != 0){
			nextNode = curnode->pageNoArray[i+1];
			find = true;
			break;
		}
	}
	//key is the smallest
	if (!find){
		nextNode = curnode->pageNoArray[0];
	}

	return nextNode;

}

// -----------------------------------------------------------------------------
// BTreeIndex::startScan
// -----------------------------------------------------------------------------
/**
 * Initialize the field to start scanning process
 * Traverse down the tree to find the first leaf node that matching the range 
 * If find the leaf node, set the next entry and currentPageNo
 * @param lowVal	Low value of range, pointer to integer / double / char string
 * @param lowOp		Low operator (GT/GTE)
 * @param highVal	High value of range, pointer to integer / double / char string
 * @param highOp	High operator (LT/LTE)
 * @throws  BadOpcodesException If lowOp and highOp do not contain one of their their expected values 
 * @throws  BadScanrangeException If lowVal > highval
 * @throws  NoSuchKeyFoundException If there is no key in the B+ tree that satisfies the scan criteria.
**/
const void BTreeIndex::startScan(const void* lowValParm,
				   const Operator lowOpParm,
				   const void* highValParm,
				   const Operator highOpParm)

{
	//initialize the field
	this->lowOp = lowOpParm;
	this->highOp = highOpParm;

	//check whether the scan is valid. Throw exception if not
	if ((lowOp != GTE && lowOp != GT) || (highOp != LTE && highOp != LT)){
		throw  BadOpcodesException();
	}

	//call the typescan template function based on attribute
	if (this->attributeType == INTEGER){
		this->lowValInt = *(int*)lowValParm;
		this->highValInt = *(int*)highValParm;
		if (this->lowValInt > this->highValInt){
			throw BadScanrangeException();
		}
		typeScan<int, NonLeafNodeInt, LeafNodeInt>(this->lowValInt, lowOpParm, this->highValInt, highOpParm);
	}
	else if (this->attributeType == DOUBLE) {
		this->lowValDouble = *(double*) lowValParm;
		this->highValDouble = *(double*)highValParm;
		if (this->lowValDouble > this->highValDouble){
			throw BadScanrangeException();
		}	
		typeScan<double, NonLeafNodeDouble, LeafNodeDouble>(this->lowValDouble, lowOpParm, this->highValDouble, highOpParm);
	}
	else {
	
		//get the first ten character of the input parameter
		char tmplow[10];
		char tmphigh[10];
		strncpy(tmplow, (char*) lowValParm, 9);
		strncpy(tmphigh, (char*) highValParm, 9);
		tmplow[9] = '\0';
		tmphigh[9] = '\0';
		this->lowValString = tmplow;
		this->highValString = tmphigh;

		if (this->lowValString > this->highValString){
			throw BadScanrangeException();
		}
		typeScan<string, NonLeafNodeString, LeafNodeString>( this->lowValString, lowOpParm, this->highValString,  highOpParm);
	}
	
	this->scanExecuting = true;
}


/**
 * The compare function to compare the relation of the two values based on op
 * @param op		GT/GTE/LT/LTE
 * @param left		The right value to be compared (int/double/string)
 * @param right		The right value to be compared (int/double/string)
**/
template<class T> 
bool BTreeIndex::compare(const Operator op, T left, T right){

	if (op == GT){
		return left < right;
	}
	else if (op == GTE){
		return left <= right;
	}
	else if (op == LT){
		return left < right;
	}
	else {
		return left <= right;
	}
}


/**
 * Called by startScan. Set the field to start scan process
 * The function is implemented in order to work for all three attributes (int/double/string)
 * The template is used to meet that goal
 * @param lowValParm	Low value of range, pointer to integer / double / char string
 * @param lowOpParm		Low operator (GT/GTE)
 * @param highValParm	High value of range, pointer to integer / double / char string
 * @param highOpParm	High operator (LT/LTE)
**/
template <class T, class nonleaf_T, class leaf_T>
const void BTreeIndex::typeScan(T lowValParm,
				   const Operator lowOpParm,
				   T highValParm,
				   const Operator highOpParm)
{		
	Page* rootpage;
	PageId pageNum = this->rootPageNum;
	
	//traverse the tree to get the leaf node
	while (true){
		
		this->bufMgr->readPage(this->file, pageNum, rootpage);
		nonleaf_T* rootnode = (nonleaf_T*) rootpage;
		this->bufMgr->unPinPage(this->file, pageNum, false);
		bool parentleaf = rootnode->level == 1;

		pageNum = rootnode->pageNoArray[0];

		if (parentleaf){
			break;
		}		
	}
	
	leaf_T* leafnode;

	//keep going right using rightsibling number
	while (true){
		this->bufMgr->readPage(this->file, pageNum, rootpage);
		leafnode = (leaf_T*) rootpage;
	
		for (int i = 0; i < this->leafOccupancy; i++){
			if (leafnode->ridArray[i].page_number == 0){
				break;
			}
			//find the first record that matches the range
			if (compare<T>(lowOpParm, lowValParm, (T) leafnode->keyArray[i]) && compare(highOpParm, (T) leafnode->keyArray[i], highValParm)){
				this->currentPageNum = pageNum;
				this->nextEntry = i;
				this->currentPageData = rootpage;
				return;
			}
		}

		//if all the key are smaller than the range minimal
		this->bufMgr->unPinPage(this->file, pageNum, false);
		if (leafnode->rightSibPageNo == 0){
			throw NoSuchKeyFoundException();
		}
		else
		{	
			pageNum = leafnode->rightSibPageNo;
		}
	}
	
}
// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------
/**
 * Set the next entry field that matches the scan range
 * Return the next record from current page being scanned. If current page has been scanned to its entirety, move on to the right sibling of current page, if any exists, to start scanning that page. Make sure to unpin any pages that are no longer required.
 * @param outRid	RecordId of next record found that satisfies the scan criteria returned in this
 * @throws ScanNotInitializedException If no scan has been initialized.
 * @throws IndexScanCompletedException If no more records, satisfying the scan criteria, are left to be scanned.
**/
const void BTreeIndex::scanNext(RecordId& outRid) 
{
	//call typescannNext based on attribute type
	if (this->scanExecuting == false) {
		throw ScanNotInitializedException();
	}

	if (this->attributeType == INTEGER){
		typescanNext<int, LeafNodeInt>(this->lowValInt, this->highValInt, outRid);
	}
	else if (this->attributeType == DOUBLE){
		typescanNext<double, LeafNodeDouble>(this->lowValDouble, this->highValDouble, outRid);
	}
	else {
		typescanNext<std::string, LeafNodeString>(this->lowValString, this->highValString, outRid);
	}

}

/**
 * Same as the scanNext. This function works for all three attribute types by using template function
 * Return the next record from current page being scanned. 
 * If current page has been scanned to its entirety, move on to the right sibling of current page, if any exists, 
 * to start scanning that page. Make sure to unpin any pages that are no longer required.
 * @param lowValParm	The low value in the range
 * @param highValParm 	The high value in the range
 * @param outRid		RecordId of next record found that satisfies the scan criteria returned in this
**/
template <class T, class leaf_T>
const void BTreeIndex::typescanNext(T lowValParm, T highValParm, RecordId& outRid) {
	leaf_T* leafnode = (leaf_T*) currentPageData;
	
	//go to the next node (rightsibling of the current node) if possible
	if (this->nextEntry >= this->leafOccupancy || leafnode->ridArray[this->nextEntry].page_number == 0){
		this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
		this->nextEntry = 0;
		this->currentPageNum = leafnode->rightSibPageNo;
		if (this->currentPageNum == 0){
			throw IndexScanCompletedException();
		} 
		else{
			this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
		}
	}
	
	leafnode = (leaf_T*) this->currentPageData;

	//check wehther the record meets the range
	if (compare<T>(this->lowOp, lowValParm, leafnode->keyArray[this->nextEntry]) && compare<T>(this->highOp, leafnode->keyArray[this->nextEntry], highValParm)){	
		outRid = leafnode->ridArray[this->nextEntry];
		this->nextEntry = nextEntry+1;
	}
	else {
		throw IndexScanCompletedException();	
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//
/**
 * Terminate the scan process by unpining used page and clear the field
 **/
const void BTreeIndex::endScan() 
{
	//terminate the current scan
	if (this->scanExecuting == false) {
		throw ScanNotInitializedException();
	}
	this->scanExecuting = false;

	//unpin page that is used during scan process
	if (this->currentPageNum != 0)
		this->bufMgr->unPinPage(this-> file, this->currentPageNum, false);
}

}
