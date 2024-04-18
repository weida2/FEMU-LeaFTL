#ifndef __DFTL_CACHE_H
#define __DFTL_CACHE_H
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "../nvme.h"
#include "ftl.h"


/*slice L2P map in page size*/
#define ENTRY_PER_PAGE (1024) // 4KB / 4B = 1024
#define INVALID_PPA     (~(0ULL))
#define INVALID_LPN     (~(0ULL))
#define UNMAPPED_PPA    (~(0ULL))

enum {
    NAND_READ =  0,
    NAND_WRITE = 1,
    NAND_ERASE = 2,

    NAND_READ_LATENCY = 40000,
    NAND_PROG_LATENCY = 200000,
    NAND_ERASE_LATENCY = 2000000,
};

typedef struct Node {
    uint64_t page_idx;
    bool update_in_flash;
    uint64_t *l2p_entries;   // pair<int, int>[0~1024] l2p
    struct Node *pre;
    struct Node *next;
    
    struct Node *hash_next;  // ??链式哈希值, 用于链式哈希的同一个%hash=key的下一个指针结点
} Node;

typedef struct LRUCache {
    uint64_t capacity;
    uint64_t count;
    
    Node *head;
    Node *tail;

    Node **hashTable; // ??链式哈希？

    uint8_t *flash_bmp;
} LRUCache;


// DFTL_MAPPING_Table
typedef struct G_map_entry
{
    uint64_t page_idx;
    uint64_t ppa;
    // plr     segs; 


} G_map_entry;

typedef struct d_l2p_tbl {
    LRUCache *CMT;    

    LRUCache *nand_cache; // GTD ?

    G_map_entry *GTD; 
    
    int num_gtd_entries; 
} DFTLTable;





Node *createNode(uint64_t page_idx);
LRUCache* createLRUCache(uint64_t capacity,uint32_t bmpCnt);
void deleteLRUCache(LRUCache* cache);
void addToFront(LRUCache *cache, Node *node);
void moveToFront(LRUCache *cache, Node *node);
void MoveToFrontByKey(LRUCache *cache, uint64_t key);
void removeByKey(LRUCache *cache, uint64_t key);

uint64_t hash(uint64_t key, uint64_t capacity);

uint64_t addNodeToCache(LRUCache *cache, LRUCache *nand_cache, Node *node, bool *evict_happen, uint64_t *evict_key ,uint64_t* lat);
uint64_t addNodeToNandCache(LRUCache *nand_cache, Node *node);

uint64_t lru_get(LRUCache *cache, uint64_t key);
void lru_put(LRUCache *cache,LRUCache *nand_cache,uint64_t key, uint64_t value ,uint64_t*lat);
uint64_t dftl_evict(LRUCache *cache, LRUCache* nand_cache,uint64_t* lat);


void dftl_node_map_set_bit(uint8_t *map, uint8_t bit_idx);
uint8_t dftl_node_map_check_bit(uint8_t *map, uint8_t bit_idx);
void dftl_node_map_clear_bit(uint8_t *map, uint8_t bit_idx);

uint64_t move_node_from_nand_to_cache(LRUCache *cache, LRUCache *nand_cache, uint64_t key,uint64_t* lat);


// DFTL_Table_method

DFTLTable* initialize_dftl_table(uint32_t lmaCnt);

uint64_t dftl_get(DFTLTable *table, uint64_t lpa ,uint64_t* lat);
void dftl_put(DFTLTable* table, uint64_t lpa, uint64_t ppa , uint64_t* lat);




#endif /* __DFTL_CACHE_H */