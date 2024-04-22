#include "dftl.h"
#include "./ftl.h"



/**
 * @brief DFTL_method
 * 
 * 创建CMT_table_Node
 */
uint64_t hash(uint64_t key, uint64_t capacity)
{
    return key % capacity;
}
void dftl_node_map_set_bit(uint8_t *map, uint8_t bit_idx) {
    *map |= (1 << bit_idx);
}

uint8_t dftl_node_map_check_bit(uint8_t *map, uint8_t bit_idx) {
    return *map & (1 << bit_idx);
}

void dftl_node_map_clear_bit(uint8_t *map, uint8_t bit_idx) {
    *map &= ~(1 << bit_idx);    
}


Node *createNode(uint64_t page_idx, struct ssd* ssd) {
    Node *new_node = (Node *)malloc(sizeof(Node));
    
    struct ppa ppa;
    ppa = get_new_page(ssd);
    uint64_t vpgidx = ppa2vpgidx(ssd, &ppa);
    ssd_advance_write_pointer(ssd);

    //ssd->rmap[vpgidx] = page_idx; // rmap.g == page.type ? 0: translate page, 1: data page 
                                    // rmap.vpn, rmap.lpn


    new_node->ppa = vpgidx;        // ppn
    new_node->page_idx = page_idx; // vpn

    new_node->next = NULL;
    new_node->pre = NULL;
    
    new_node->l2p_entries = (uint64_t *)malloc(sizeof(uint64_t) * ENTRY_PER_PAGE);
    memset(new_node->l2p_entries, 0xff, sizeof(uint64_t) * ENTRY_PER_PAGE); // 不能将每个元素设置成0xFFFFFFFF
    
    new_node->update_in_flash = true;
    new_node->hash_next = NULL;
    
    return new_node;
}


// 创建LRU 采用哈希链式结构
LRUCache* createLRUCache(uint64_t capacity, uint32_t bmpCnt) {
    LRUCache *cache = (LRUCache *)malloc(sizeof(LRUCache));
    cache->capacity = capacity;
    cache->count = 0;
    cache->head = NULL;
    cache->tail = NULL;

    cache->hashTable = (Node **)malloc(sizeof(Node *) * capacity);
    for (int i = 0; i < capacity; i++) {
        cache->hashTable[i] = NULL;
    }

    // 初始化位图 用来判断该页在 闪存上 或者 cache 上
    cache->flash_bmp = (uint8_t *)malloc(sizeof(uint8_t) * bmpCnt);
    memset(cache->flash_bmp, 0, sizeof(uint8_t) * bmpCnt);

    assert(cache->flash_bmp[0] == (uint8_t)0);

    return cache;
}


void deleteLRUCache(LRUCache* cache) {
    for (int i = 0; i < cache->capacity; i++) {
        Node *node = cache->hashTable[i];
        while (node != NULL) {
            Node *node_next = node->hash_next;
            free(node);
            node = node_next;
        }
    }
    free(cache->hashTable);
    free(cache->flash_bmp);
    free(cache);
}

// 头插法 新插入的直接addtoFront
void addToFront(LRUCache *cache, Node *node) {
    node->next = cache->head;
    node->pre = NULL;

    if (cache->head != NULL) {
        cache->head->pre = node;
    }
    cache->head = node;
    if (cache->tail == NULL) {
        cache->tail = node;
    }
}

// 已经存在的 moveToFront
void moveToFront(LRUCache *cache, Node *node) {
    if (node == cache->head) {
        return ;
    }
    if (node == cache->tail) {
        cache->tail = node->pre;
    }

    node->pre->next = node->next;
    if (node->next != NULL) {
        node->next->pre = node->pre;
    }
    addToFront(cache, node);
}

// tmp no use
void MoveToFrontByKey(LRUCache *cache, uint64_t key) {
    uint64_t page_idx = key / ENTRY_PER_PAGE;
    uint64_t index = hash(page_idx, cache->capacity);

    Node *cur = cache->hashTable[index];
    while (cur != NULL) {
        if (cur->page_idx == page_idx) {
            moveToFront(cache, cur);
            break;
        }
        cur = cur->hash_next;
    }

    return ;
}

// TODO: tmp don't know to use
void removeByKey(LRUCache *cache, uint64_t key) {
    uint64_t page_idx = key / ENTRY_PER_PAGE;
    uint64_t index = hash(page_idx, cache->capacity);

    Node *pre = NULL;
    Node *cur = cache->hashTable[index];

    // first: update hash[key]_list
    while (cur != NULL) {
        if (cur->page_idx == page_idx) {
            break;
        }
        pre = cur;
        cur = cur->hash_next;
    }
    if (cur == NULL) {
        // node with key no exist in cache
        return ;
    }
    if (pre == NULL) {
        // has_head node == page_idx
        cache->hashTable[index] = cache->hashTable[index]->hash_next;
    }else {
        pre->hash_next = cur->hash_next;
    }

    // second: update LRU node list (two-dir)
    if (cache->head == cur) {
        cache->head = cur->next;
    }
    if (cache->tail == cur) {
        cache->tail = cur->pre;
    }
    if (cur->pre != NULL) {
        cur->pre->next = cur->next;
    }
    if (cur->next != NULL) {
        cur->next->pre = cur->pre;
    }

    free(cur->l2p_entries);
    free(cur);

    cache->count--;
}

// 新的Node-> page_idx 不管新的旧的
uint64_t addNodeToCache(DFTLTable *d_maptbl, LRUCache *cache, LRUCache *nand_cache, Node *node, bool *evict_happen, uint64_t *evict_key ,uint64_t* lat) {
    if (cache->count >= cache->capacity) {
        *evict_happen = true;
        *evict_key = dftl_evict(d_maptbl, cache, nand_cache, lat);
    }

    uint64_t index = hash(node->page_idx, cache->capacity);
    // 头插法到hastbl[i]中
    Node *cur = cache->hashTable[index];
    node->hash_next = cur;
    cache->hashTable[index] = node;

    // 将新插入的node 放到LRU_list 的top
    addToFront(cache, node);
    cache->count++;

    // update_map
    uint64_t bmp_idx = node->page_idx / 8;
    uint8_t  bit_idx  = node->page_idx % 8;
    dftl_node_map_set_bit(&(cache->flash_bmp[bmp_idx]), bit_idx);
    
    return node->page_idx;
}

// TODO: full gc
// same node.page_idx remove
uint64_t addNodeToNandCache(LRUCache *nand_cache, Node *node) {
    uint64_t index = hash(node->page_idx, nand_cache->capacity); // index = node->page_idx;

    Node *cur = nand_cache->hashTable[index];
    node->hash_next = cur;
    nand_cache->hashTable[index] = node;

    addToFront(nand_cache, node);
    nand_cache->count++;

    return node->page_idx;
}

uint64_t lru_get(DFTLTable *d_maptbl, LRUCache *cache, uint64_t key, struct ssd* ssd) {
    uint64_t page_idx = key / ENTRY_PER_PAGE;
    uint64_t index = hash(page_idx, cache->capacity);

    Node *cur = cache->hashTable[index];
    while (cur != NULL) {
        if (cur->page_idx == page_idx) {
            moveToFront(cache, cur);
            return cur->l2p_entries[key % ENTRY_PER_PAGE];  // ppa = l2p[oft]
        }
        cur = cur->hash_next;
    }

    // lpn还未写入 读缺失
    return UNMAPPED_PPA;
}

// (key, value) -> (lpa, ppa)
void lru_put(DFTLTable *d_maptbl, LRUCache *cache,LRUCache *nand_cache,uint64_t key, uint64_t value ,uint64_t*lat, struct ssd* ssd) {
    uint64_t page_idx = key / ENTRY_PER_PAGE;
    uint64_t hash_index = hash(page_idx, cache->capacity);

    femu_log("[lru_put] lpa: %lu-> ppa: %lu, page_idx: %lu, cache_size: %lu\n", key, value, page_idx, cache->capacity);

    Node *cur = cache->hashTable[hash_index];
    while (cur != NULL) {
        if (cur->page_idx == page_idx) {
            // this node([page_idx][l2p_0~1024]) upadte, should be upadte in flash
            cur->l2p_entries[key % ENTRY_PER_PAGE] = value;
            cur->update_in_flash = true;
            moveToFront(cache, cur);
            return ;
        }
        cur = cur->hash_next;
    }
    // this page_idx has not been created (New)
    Node *new_node = createNode(page_idx, ssd);
    new_node->l2p_entries[key % ENTRY_PER_PAGE] = value;
    bool evict_happen = false;
    uint64_t evict_key = UNMAPPED_PPA;

    addNodeToCache(d_maptbl, cache, nand_cache, new_node, &evict_happen, &evict_key, lat);

    return ;
}

// pop()  use lru strategy: least recently use node in cache
uint64_t dftl_evict(DFTLTable *d_maptbl, LRUCache *cache, LRUCache* nand_cache,uint64_t* lat) {
    if (cache->count >= cache->capacity) {
        Node *lastNode = cache->tail;
        if (lastNode) {
            d_maptbl->counter.evit_cnt++;

            uint64_t has_index = hash(lastNode->page_idx, cache->capacity);
            Node *cur = cache->hashTable[has_index];
            Node *pre = NULL;

            while (cur != NULL && cur != lastNode) {
                pre = cur;
                cur = cur->hash_next;
            }
            if (cur == lastNode) {
                // hash_list one-dir
                if (pre) {
                    pre->hash_next = cur->hash_next;
                }else {
                    cache->hashTable[has_index] = cur->hash_next;
                }
                // Lru_node_list two-dir
                if (lastNode->pre) {
                    lastNode->pre->next = lastNode->next;
                }
                if (lastNode->next) {
                    lastNode->next->pre = lastNode->pre;
                }
                cache->tail = lastNode->pre;

                uint64_t page_idx = lastNode->page_idx;      
                
                if (lastNode->update_in_flash) {
                    (*lat) += NAND_PROG_LATENCY;           // don't use lat += ssd_advance_status(ssd, &ppa, &swr); 模拟机制 should be used?
                    
                    // updae to GTD
                    d_maptbl->GTD[page_idx].page_idx    = page_idx; 
                    d_maptbl->GTD[page_idx].ppa         = lastNode->ppa;
                    d_maptbl->GTD[page_idx].segs.no_use = 0; //todo
    
                    addNodeToNandCache(nand_cache, lastNode);  // movetoflahs should add lat


                    // upate bitmap
                    uint64_t bmp_idx = lastNode->page_idx / 8;
                    uint8_t  bit_idx  = lastNode->page_idx % 8;
                    dftl_node_map_clear_bit(&(cache->flash_bmp[bmp_idx]), bit_idx);
                }

                cache->count--;
                
                return page_idx;
            }
        }
    }
    femu_log("EVICT FAILED\n");
    return UNMAPPED_PPA;
}


// 从 flash移出node应该保留 ？
uint64_t move_node_from_nand_to_cache(DFTLTable *d_maptbl, LRUCache *cache, LRUCache *nand_cache, uint64_t key, uint64_t* lat) {
    uint64_t page_idx = key / ENTRY_PER_PAGE;
    uint64_t has_index = hash(page_idx, nand_cache->capacity);
    
    Node *cur = nand_cache->hashTable[has_index];
    Node *pre = NULL;

    while (cur != NULL) {
        if (cur->page_idx == page_idx) {

            break;
        }
        pre = cur;
        cur = cur->hash_next;
    }
    if (cur == NULL) {
        return UNMAPPED_PPA;
    }

    // should not update
    // update has_linklist one-dir
    if (pre == NULL) {
        nand_cache->hashTable[has_index] = nand_cache->hashTable[has_index]->hash_next;
    } else {
        pre->hash_next = cur->hash_next;
    }

    // upate LRU_linklist two-dir
    if (nand_cache->head == cur) {
        nand_cache->head = cur->next;
    }
    if (nand_cache->tail == cur) {
        nand_cache->tail = cur->pre;
    }

    if (cur->pre != NULL) {
        cur->pre->next = cur->next;
    }
    if (cur->next != NULL) {
        cur->next->pre = cur->pre;
    }

    uint64_t evict_key = UNMAPPED_PPA;
    bool evict_happen = false;
    (*lat) += NAND_READ_LATENCY;   // don't use lat += ssd_advance_status(ssd, &ppa, &swr); 模拟机制 should not be used?
    cur->update_in_flash = false; 

    addNodeToCache(d_maptbl, cache, nand_cache, cur, &evict_happen, &evict_key, lat);

    // uint64_t bmp_idx = cur->page_idx / 8;
    // uint8_t  bit_idx = cur->page_idx % 8;
    // dftl_node_map_clear_bit(&(cache->flash_bmp[bmp_idx]), bit_idx);

    //nand_cache->count--;

    return cur->page_idx; // 
}

DFTLTable* dftl_table_init(uint32_t tt_pages) {
    DFTLTable *table = (DFTLTable *)malloc(sizeof(DFTLTable));
    assert(table != NULL);
    uint64_t total_vpn = (tt_pages + ENTRY_PER_PAGE - 1) / ENTRY_PER_PAGE;
    uint64_t bmp_cnt = (total_vpn + 7) / 8;  // tt_pages / 1024 / 8
    
    table->CMT = createLRUCache(total_vpn * CMT_ratio, bmp_cnt);

    table->GTD = (G_map_entry *)malloc(sizeof(G_map_entry) * total_vpn);
    for (int i = 0; i < total_vpn; i++) {
        table->GTD[i].page_idx = i;
        table->GTD[i].ppa      = UNMAPPED_PPA;
        table->GTD[i].segs.no_use = 0;
    }

    table->nand_cache = createLRUCache(total_vpn, bmp_cnt);

    table->counter.group_write_cnt   = 0;
    table->counter.group_read_cnt    = 0;
    table->counter.group_read_miss   = 0;
    table->counter.group_cmt_hit     = 0;
    table->counter.group_double_read = 0;
    table->counter.evit_cnt          = 0;

    return table;
}

uint64_t dftl_get(DFTLTable *table, uint64_t lpa ,uint64_t* lat, struct ssd* ssd) {
    uint64_t page_idx = lpa / ENTRY_PER_PAGE;
    uint64_t bmp_idx = page_idx / 8;
    uint8_t  bit_idx = page_idx % 8;


    if (!dftl_node_map_check_bit(&(table->CMT->flash_bmp[bmp_idx]), bit_idx)) {
        uint64_t evit_page_idx = move_node_from_nand_to_cache(table, table->CMT, table->nand_cache, lpa, lat);
        if (evit_page_idx == UNMAPPED_PPA) {
            return evit_page_idx;
        }
        table->counter.group_double_read++;
    } else {
        table->counter.group_cmt_hit++;
    }

    u_int64_t ppa = lru_get(table, table->CMT, lpa, ssd);
    femu_log("[dftl_get] lpa: %lu, page_idx: %lu, ret_ppa:%ld\n", lpa, page_idx, ppa == UNMAPPED_PPA ? -1 : ppa);

    return ppa;

}
void dftl_put(DFTLTable* table, uint64_t lpa, uint64_t ppa , uint64_t* lat, struct ssd* ssd) {
    // 写入的时候不用检查CMT是否有该页，直接写
    // uint64_t page_idx = lpa / ENTRY_PER_PAGE;
    // uint64_t bmp_idx = page_idx / 8;
    // uint8_t  bit_idx = page_idx % 8;
    // if (dftl_node_map_check_bit(&table->CMT->flash_bmp[bmp_idx], bit_idx)) {
    //     move_node_from_nand_to_cache(table, table->nand_cache, lpa, lat);
    // }
    lru_put(table, table->CMT, table->nand_cache, lpa, ppa, lat, ssd);
    return ;
}

void dftl_static(DFTLTable *d_maptbl) {
    femu_log("[DFTL]write_cnt:%d, read_cnt: %d, read_miss: %d, cmt_hit: %d, double_read: %d, evit_cnt: %d\n",
    d_maptbl->counter.group_write_cnt, d_maptbl->counter.group_read_cnt,
    d_maptbl->counter.group_read_miss, d_maptbl->counter.group_cmt_hit,
    d_maptbl->counter.group_double_read, d_maptbl->counter.evit_cnt);
    femu_log("[DFTL]cache_cnt: %lu/%lu, nand_cnt: %lu/%lu\n", d_maptbl->CMT->count, d_maptbl->CMT->capacity, 
    d_maptbl->nand_cache->count, d_maptbl->nand_cache->capacity);
}