#include"dftl.h"

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


Node *createNode(uint64_t page_idx) {
    Node *new_node = (Node *)malloc(sizeof(Node));
    
    new_node->page_idx = page_idx;
    new_node->next = NULL;
    new_node->pre = NULL;
    
    new_node->l2p_entries = (uint64_t *)malloc(sizeof(uint64_t) * ENTRY_PER_PAGE);
    memset(new_node->l2p_entries, 0xff, sizeof(uint64_t) * ENTRY_PER_PAGE); // 不能将每个元素设置成0xFFFFFFFF
    
    new_node->update_in_flash = true;
    new_node->hash_next = NULL;
    
    return new_node;
}


// 创建LRU 采用哈希链式结构
LRUCache* createLRUCache(uint64_t capacity,uint32_t bmpCnt) {
    LRUCache *cache = (LRUCache *)malloc(sizeof(LRUCache));
    cache->capacity = capacity;
    cache->count = 0;
    cache->head = NULL;
    cache->tail = NULL;

    cache->hashTable = (Node **)malloc(sizeof(Node *) * capacity);
    for (int i = 0; i < capacity; i++) {
        cache->hashTable[i] = NULL;
    }

    // 初始化位图 用来判断该页是否在闪存上?
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
uint64_t addNodeToCache(LRUCache *cache, LRUCache *nand_cache, Node *node, bool *evict_happen, uint64_t *evict_key ,uint64_t* lat) {
    if (cache->count >= cache->capacity) {
        *evict_happen = true;
        *evict_key = dftl_evict(cache, nand_cache, *lat);
    }

    uint64_t index = hash(node->page_idx, cache->capacity);
    // 头插法到hastbl[i]中
    Node *cur = cache->hashTable[index];
    node->hash_next = cur;
    cache->hashTable[index] = node;

    // 将新插入的node 放到LRU_list 的top
    addToFront(cache, node);
    cache->count++;
    
    return node->page_idx;
}

// TODO: 满了GC ？
uint64_t addNodeToNandCache(LRUCache *nand_cache, Node *node) {
    uint64_t index = hash(node->page_idx, nand_cache->capacity); // index = node->page_idx;

    Node *cur = nand_cache->hashTable[index];
    node->hash_next = cur;
    nand_cache->hashTable[index] = node;

    addToFront(nand_cache, node);
    nand_cache->count++;

    return node->page_idx;
}

uint64_t lru_get(LRUCache *cache, uint64_t key) {
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
void lru_put(LRUCache *cache,LRUCache *nand_cache,uint64_t key, uint64_t value ,uint64_t*lat) {
    uint64_t page_idx = key / ENTRY_PER_PAGE;
    uint64_t index = hash(page_idx, cache->capacity);

    Node *cur = cache->hashTable[index];
    while (cur != NULL) {
        if (cur->page_idx == page_idx) {
            cur->l2p_entries[key % ENTRY_PER_PAGE] = value;
            // this node([page_idx][l2p_0~1024]) upadte, should be upadte in flash
            cur->update_in_flash = true;
            moveToFront(cache, cur);
            return ;
        }
        cur = cur->hash_next;
    }
    // this page_idx has not been created (New)
    Node *new_node = createNode(page_idx);
    new_node->l2p_entries[key % ENTRY_PER_PAGE] = value;
    bool evict_happen = false;
    bool evict_key = UNMAPPED_PPA;

    addNodeToCache(cache, nand_cache, new_node, evict_happen, evict_key, lat);

    return ;
}

// pop()  use lru strategy: least recently use node in cache
uint64_t dftl_evict(LRUCache *cache, LRUCache* nand_cache,uint64_t* lat) {
    if (cache->count >= cache->capacity) {
        Node *lastNode = cache->tail;
        if (lastNode) {
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
                    cache->hashTable[has_index] = cache->hashTable[has_index]->hash_next;
                }
                // Lru_node_list two-dir
                if (lastNode->pre) {
                    lastNode->pre->next = lastNode->next;
                }
                if (lastNode->next) {
                    lastNode->next->pre = lastNode->pre;
                }
                cache->tail = lastNode->pre;

                uint64_t result = lastNode->page_idx;
                addNodeToNandCache(nand_cache, lastNode);  // movetoflahs should add lat
                if (lastNode->update_in_flash) {
                    (*lat) += NAND_PROG_LATENCY;
                }

                // upate bitmap
                uint64_t bmp_idx = lastNode->page_idx / 8;
                uint8_t  bit_idx  = lastNode->page_idx % 8;
                dftl_node_map_set_bit(&(cache->flash_bmp[bmp_idx]), bit_idx);
                cache->count--;
                
                return result;
            }
        }
    }
    femu_log("EVICT FAILED\n");
    return UNMAPPED_PPA;
}


// 从 flash移出node应该保留 ？
uint64_t move_node_from_nand_to_cache(LRUCache *cache, LRUCache *nand_cache, uint64_t key,uint64_t* lat) {
    uint64_t page_idx = key / ENTRY_PER_PAGE;
    uint64_t has_index = has(page_idx, nand_cache->capacity);
    
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
        nand_cache->tail = cur.pre;
    }

    if (cur->pre != NULL) {
        cur->pre->next = cur->next;
    }
    if (cur->next != NULL) {
        cur->next->pre = cur->pre;
    }

    uint64_t evict_key = UNMAPPED_PPA;
    bool evict_happen = false;
    (*lat) += NAND_READ_LATENCY;
    addNodeToCache(cache, nand_cache, cur, &evict_happen, &evict_key, lat);
    cur->update_in_flash = false; // should be before ??

    uint64_t bmp_idx = cur->page_idx / 8;
    uint8_t  bit_idx = cur->page_idx % 8;
    dftl_node_map_clear_bit(&(cache->flash_bmp[bmp_idx]), bit_idx);

    nand_cache->count--;

    return evict_key; // ?? no cur->page_idx ?
}

DFTLTable* initialize_dftl_table(uint32_t lmaCnt) {

}

uint64_t dftl_get(DFTLTable *table, uint64_t lpa ,uint64_t* lat) {

}
void dftl_put(DFTLTable* table, uint64_t lpa, uint64_t ppa , uint64_t* lat) {
    
}