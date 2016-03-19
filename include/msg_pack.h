#ifndef MSG_PACK_H
#define MSG_PACK_H
#include <iostream>
#include <queue>
#define NULLTAG 1
#define CELLTAG 2
#define CHILDCELLTAG 3
#define LEVELTAG 7
#define BODYTAG 8
#define FLUSHTAG 9
#define MAXTAG  15
#define LEVELSHIFT 5
#define REQUESTSHIFT 4
#define DIRECTIONSHIFT 1
#define GRAINSHIFT 16
#define LEVELMASK 0x1F
#define REQUESTMASK 0xF
#define DIRECTIONMASK 0x1
#define GRAINMASK 0xFFFF
#define RAWPTR(vec) &vec[0]
#define TOGGLEDIRECTION(tag) tag ^= DIRECTIONMASK;

namespace exafmm {
	Cell __cell___; 
	Body __body___;
	const int CELLWORD = sizeof(__cell___) / 4;
	const int BODYWORD = sizeof(__body___) / 4;
	const uint8_t SENDBIT = 1;
	const uint8_t RECEIVEBIT = 0;

	inline int encryptMessage(uint16_t grainSize, uint8_t requestType, uint8_t level, uint8_t direction) {
		int tag = int(grainSize);
		tag<<=REQUESTSHIFT; tag|=requestType;        
	  tag<<=LEVELSHIFT; tag|=level;
	  tag<<=DIRECTIONSHIFT; tag|= direction;
	  return tag;  
	}

	inline void decryptMessage(int tag, uint16_t& grainSize, uint8_t& requestType, uint8_t& level, uint8_t& direction) {
		direction = tag & DIRECTIONMASK;
		tag >>= DIRECTIONSHIFT;
		level = tag & LEVELMASK;
		tag >>= LEVELSHIFT;
		requestType = tag & REQUESTMASK;	
		tag >>= REQUESTSHIFT;
		grainSize = tag & GRAINMASK;
	}

	inline uint16_t getGrainSize(int const& tag) {
		return ((tag >> (REQUESTSHIFT + DIRECTIONSHIFT + LEVELSHIFT)) & GRAINMASK);
	}

	inline uint8_t getMessageType(int const& tag) {
		return ((tag >> (DIRECTIONSHIFT + LEVELSHIFT)) & REQUESTMASK);
	}

	inline uint8_t getMessageLevel(int const& tag) {
		return ((tag >> DIRECTIONSHIFT) & LEVELMASK);
	}

	inline uint8_t getMessageDirection(int const& tag) {
		return (tag & DIRECTIONMASK);
	}


#if DFS
	inline void packSubtree(Cells& cells, C_iter const& C0, Cell const& C, uint16_t const& grainSize, size_t& index) {
		C_iter begin = C0 + C.ICHILD;
	  C_iter end   = C0 + C.ICHILD + C.NCHILD;
	  cells.insert(cells.end(),begin,end);
	  index+=C.NCHILD;
		for(C_iter cc = begin; cc<end; ++cc) 
			if(index < grainSize) packSubtree(cells,C0,*cc,grainSize,index);				
	}
#else
	inline void packSubtree(Cells& cells, C_iter const& C0, Cell const& root, uint16_t const& grainSize, size_t& index) {		
		std::queue<Cell> cell_queue;
		C_iter begin = C0 + root.ICHILD;  
		int size = root.NCHILD;
		for (int i = 0; i < size; ++i) 
			cell_queue.push(*(begin + i));	
		while(cell_queue.size() > 0) {
			Cell C = cell_queue.front();		
			begin = C0 + C.ICHILD;  	
			size = C.NCHILD;		
			C.ICELL = index;
			cells.push_back(C);
			cell_queue.pop();
			index++;
			if(index < grainSize) 
				for (int i = 0; i < size; ++i) {										  			 
					cell_queue.push(*(begin + i));	
					cell_queue.back().IPARENT = index - 1;	
				}			
		}	
	}
#endif
} // end exafmm

#endif