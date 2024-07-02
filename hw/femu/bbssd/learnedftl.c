#include "learnedftl.h"
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
    for (int i = 0; i < ENTRY_PER_PAGE; i++) {
        new_node->l2p_entries[i] = UNMAPPED_PPA;
    }
    
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

    //femu_log("[lru_put] lpa: %lu-> ppa: %lu, page_idx: %lu, cache_size: %lu\n", key, value, page_idx, cache->capacity);

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

/**
 * @brief dftl_evict
 * 
 * pop()  use lru strategy: least recently use node in cache
 * Learned_FTL for learn: 可以有多种更新的学习策略在CMT_evict的时候，(好时更，用占用大小来做阈值评判标准，) 
 * 更新的时候如果已经有该范围内的索引，则涉及到索引的更新，(就地更新，延迟更新，增量更新)
 * 如果没有该范围的索引则训练拟合，涉及索引模型的选择，(线性模型（线性回归（多种算法，最优算法，贪心线性算法），线性不回归）和非线性模型(神经网络模型NN(GPU更利于这种矩阵数学并行运算), CDF, B+树..)和混合模型)
 * 主要区别:拟合预测的精度，和拟合模型的大小和训练耗时，(预测错误有开销，要想精度高则训练的模型大小更大以及耗时更大和处理器需求)
 * 
 */


uint64_t dftl_evict(DFTLTable *d_maptbl, LRUCache *cache, LRUCache* nand_cache,uint64_t* lat) {
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
                    d_maptbl->counter.evit_cnt++;
                    
                    // updae GTD, split point
                    Point points[ENTRY_PER_PAGE];
                    int num_points = 0;

                    for (uint64_t i = 0; i < ENTRY_PER_PAGE; i++) {
                        uint64_t lpn, ppn;
                        lpn = i; // lpn % ENTRY_PER_PAGE
                        ppn = lastNode->l2p_entries[i];

                        if (ppn != UNMAPPED_PPA) {
                            Point pt;
                            pt.x = lpn, pt.y = ppn;
                            points[num_points++] = pt;
                        }
                    }

                    if (d_maptbl->swt_io) {
                        d_maptbl->GTD[page_idx].page_idx = page_idx;
                        d_maptbl->GTD[page_idx].ppa      = lastNode->ppa;

                        dftl_GTD_update(&d_maptbl->GTD[page_idx], points, num_points);
                    }


                    addNodeToNandCache(nand_cache, lastNode);  // movetoflahs should add lat
                    (*lat) += NAND_PROG_LATENCY;           // don't use lat += ssd_advance_status(ssd, &ppa, &swr); 模拟机制 should be used?

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


// 从 flash移出node应该保留
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
        // 这里应该读缺失时也要更新延迟，缓存映射表的缺点，(可以用位图标记?)
        (*lat) += NAND_READ_LATENCY;   // don't use lat += ssd_advance_status(ssd, &ppa, &swr); 模拟机制 should not be used?
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

    return cur->page_idx; 
}

DFTLTable* dftl_table_init(uint32_t tt_pages) {
    DFTLTable *table = (DFTLTable *)malloc(sizeof(DFTLTable));
    assert(table != NULL);

    uint64_t total_vpn = (tt_pages + ENTRY_PER_PAGE - 1) / ENTRY_PER_PAGE;
    uint64_t bmp_cnt = (total_vpn + 7) / 8;  // tt_pages / 1024 / 8
    
    table->CMT = createLRUCache(total_vpn * CMT_ratio, bmp_cnt);

    table->GTD = (G_map_entry *)malloc(sizeof(G_map_entry) * total_vpn);
    for (uint64_t i = 0; i < total_vpn; i++) {
        gtd_entry_init(&table->GTD[i], i);
    }

    table->nand_cache = createLRUCache(total_vpn, bmp_cnt);

    table->counter.group_write_cnt   = 0;
    table->counter.group_read_cnt    = 0;
    table->counter.group_read_miss   = 0;
    table->counter.group_cmt_hit     = 0;
    table->counter.group_double_read = 0;
    table->counter.evit_cnt          = 0;
    table->counter.GTD_hit           = 0;

    table->swt_io = 0;  // 0:DFTL, 1:DFTL_LEA

    return table;
}

uint64_t dftl_get(DFTLTable *table, uint64_t lpa ,uint64_t* lat, struct ssd* ssd) {
    uint64_t page_idx = lpa / ENTRY_PER_PAGE;
    uint64_t bmp_idx = page_idx / 8;
    uint8_t  bit_idx = page_idx % 8;

    if (!dftl_node_map_check_bit(&(table->CMT->flash_bmp[bmp_idx]), bit_idx)) {
        if (table->swt_io) {
            // check GTD bit_map
            uint64_t ret_ppa = gtd_lookup(&table->GTD[page_idx], lpa);
            // check flash
            if (!ret_ppa) {
                uint64_t evit_page_idx = move_node_from_nand_to_cache(table, table->CMT, table->nand_cache, lpa, lat);
                if (evit_page_idx == UNMAPPED_PPA) {
                    //这里应该也是得加上一次额外的多余的闪存读取，读缺失惩罚比较大
                    return evit_page_idx;
                }
                table->counter.group_double_read++;
            } else {
                table->counter.GTD_hit++;
                return ret_ppa;
            }
        } else {
            uint64_t evit_page_idx = move_node_from_nand_to_cache(table, table->CMT, table->nand_cache, lpa, lat);
            if (evit_page_idx == UNMAPPED_PPA) {
                return evit_page_idx;
            }
            table->counter.group_double_read++; 
        }
    } else {
        table->counter.group_cmt_hit++;
    }

    u_int64_t ppa = lru_get(table, table->CMT, lpa, ssd);
    // femu_log("[dftl_get] lpa: %lu, page_idx: %lu, ret_ppa:%ld\n", lpa, page_idx, ppa == UNMAPPED_PPA ? -1 : ppa);
    
    // 发现一个问题
    // 如果此时ppa == UNMAP说明 读的该lpa所在页号page_idx是存在，但是对应的lpa不存在(访问不存在的lpa)，粗粒度读缺失，造成多余的闪存访问延迟

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



// PLR_method
int dftl_Segment_is_valid(Segment *seg, uint32_t x) {
    if (!(x >= seg->x1 && x <= seg->x2))
        return 0;
    if (seg->accurate) {
        int k = round(1.0/seg->k);
        if (k != 0) {
            if (((x - seg->x1) % k) != 0)
                return 0;
        } else {
            // TODO
        }
    }
    else {
        if (seg->filter.length > 0 && (x - seg->x1) <= seg->filter.length) 
            return seg->filter.bitmap[x - seg->x1];
    }
    return 1;
}



uint32_t dftl_Segment_gety(Segment *seg, bool check, uint32_t x) {
    int predict = 0;
    if (!check || dftl_Segment_is_valid(seg, x)) {
            predict = (int)(x*seg->k + seg->b);
            return predict;
    }
    return predict;
}

// 判断段属性: 精确性
void dftl_Segment_check_properties(Segment *seg, Point *points, int num_points) {
    seg->accurate = true, seg->consecutive = true;
    
    for (int i = 0; i < num_points; i++) {
        uint32_t ppa;
        ppa = dftl_Segment_gety(seg, false, points[i].x);
        if (ppa != points[i].y) {
            seg->accurate = false;
        }
        if (!seg->accurate) return ;
    }
}


void dftl_Segment_init(Segment *seg, double k, double b, int x1, int x2, Point *points, int num_points) {
    seg->k = k;
    seg->b = b;
    seg->x1 = x1;
    seg->x2 = x2;
    seg->accurate = true;
      
    seg->filter.length = 0;
    seg->filter.bitmap = NULL;


    if (points != NULL) {
        dftl_Segment_check_properties(seg, points, num_points);

        seg->filter.length = seg->x2 - seg->x1 + 1;
        seg->filter.bitmap = (unsigned char *)malloc(seg->filter.length * sizeof(unsigned char));
        memset(seg->filter.bitmap, 0, seg->filter.length * sizeof(unsigned char));

        if (seg->filter.bitmap != NULL) {
            memset(seg->filter.bitmap, 0, seg->filter.length * sizeof(unsigned char));
            for (int i = 0; i < num_points; i++) {
                seg->filter.bitmap[points[i].x - seg->x1] = 1;
            } 
        }
    }
}

void dftl_Segment_merge(Segment *new_seg, Segment *old_seg, int *samelevel) {
    if (0 && DeBUG) femu_log("new_seg.x1: %u, x2: %u, old_seg.x1: %u, old_seg.x2: %u\n", 
                                    new_seg->x1, new_seg->x2, old_seg->x1, old_seg->x2);

    int start = L_MIN(new_seg->x1, old_seg->x1);
    int end   = L_MAX(new_seg->x2, old_seg->x2);

    if (start == new_seg->x1 && end == new_seg->x2) {
        *samelevel = -1;
        return ;
    }

    if (new_seg->x1 > old_seg->x1 && new_seg->x1 <= old_seg->x2) {
        old_seg->x2 = new_seg->x1 - 1;
        *samelevel = 1;
        
        int k = 0;
        int i = 0;
        for (i = 0; i < old_seg->x2 - old_seg->x1 + 1; i++) {
            old_seg->filter.bitmap[i] = old_seg->filter.bitmap[k++]; 
        }
        while (k < old_seg->filter.length) old_seg->filter.bitmap[k++] = 0;
        return ;
    }

    if (old_seg->x1 > new_seg->x1 && old_seg->x1 <= new_seg->x2) {
        *samelevel = 1;
        
        int k = new_seg->x2 - old_seg->x1 + 1;
        int i = 0;
        old_seg->x1 = new_seg->x2 + 1;
        for (i = 0; i < old_seg->x2 - old_seg->x1 + 1; i++) {
            if (k == old_seg->filter.length) break;
            old_seg->filter.bitmap[i] = old_seg->filter.bitmap[k++];
        }
        while (i < old_seg->filter.length) old_seg->filter.bitmap[i++] = 0;
        return ;
    }
    *samelevel = 0;   // 存在仍重叠 
    return ;

}



/**
 * @brief SimpleSegment_method
 * 
 */

int dftl_get_y(SimpleSegment *simpleseg, int x) {
    int predict;
    predict = round(x * simpleseg->k + simpleseg->b);
    return predict;
}


InsecPoint dftl_inter_section(SimpleSegment *s1, SimpleSegment *s2) {
    InsecPoint insec_pt;
    insec_pt.x = 0, insec_pt.y = 0;
    // fix bug 如果两条直线平行或重叠是不会有交点的
    if (s1->k != s2->k) {
        insec_pt.x = (double) ((s2->b - s1->b) / (s1->k - s2->k));
        insec_pt.y = (double) ((s1->k * s2->b - s2->k * s1->b) / (s1->k - s2->k));
    }
    return insec_pt;
}

bool dftl_is_above(Point *pt, SimpleSegment *s) {
    return pt->y > (int)(s->k * pt->x + s->b);
}

bool dftl_is_below(Point *pt, SimpleSegment *s) {
    return pt->y < (int)(s->k * pt->x + s->b);
}

Point dftl_get_upper_bound(Point *pt, double gamma) {
    Point p;
    p.x = pt->x, p.y = pt->y + gamma;
    return p;
}

Point dftl_get_lower_bound(Point *pt, double gamma) {
    Point p;
    p.x = pt->x, p.y = pt->y - gamma;
    return p;
}

SimpleSegment dftl_frompoints(Point p1, Point p2) {
    SimpleSegment simplesegment;
    if (p2.x != p1.x) {
        simplesegment.k = (double)((double)(int32_t)(p2.y - p1.y) / (p2.x - p1.x));
    }
    if (p2.x != p1.x) {
        simplesegment.b = (double)((double)(p1.y * p2.x) - (double)(p1.x * p2.y)) / (p2.x - p1.x);
    }
    simplesegment.x1 = p1.x;
    simplesegment.x2 = p2.x;
    return simplesegment;
}

SimpleSegment dftl_frompoints_insec(InsecPoint p1, Point p2) {
    SimpleSegment simplesegment;
    if (p2.x != p1.x) {
        simplesegment.k = (double)((double)((int32_t)p2.y - p1.y) / (p2.x - p1.x));
    }
    if (p2.x != p1.x) {
        simplesegment.b = (double)((double)(p1.y * p2.x) - (double)(p1.x * p2.y)) / (p2.x - p1.x);
    }
    simplesegment.x1 = p1.x;
    simplesegment.x2 = p2.x;
    return simplesegment;
}


/**
 * @brief PLR_method
 * 
 */
void dftl_plr_init(PLR* plr, double gamma) {
    plr->gamma = gamma;
    plr->max_length = ENTRY_PER_PAGE;

    dftl_plr__init(plr);
}


// Tem_struct, use for build segments, and then destry after insert
void dftl_plr__init(PLR* plr) {

    // plr_destroy(plr);

    plr->segments = g_malloc0(sizeof(Segment) * (plr->max_length));
    for (int i = 0; i < plr->max_length; i++) {
        dftl_Segment_init(&plr->segments[i], 0, 0, 0, 0, NULL, 0);
    }
    plr->num_segments = 0;

    plr->s0.x = 0, plr->s0.y = 0;   // init point
    plr->s1.x = 0, plr->s1.y = 0;
    plr->rho_upper.k = 0, plr->rho_upper.b = 0, plr->rho_upper.x1 = 0, plr->rho_upper.x2 = 0;
    plr->rho_lower.k = 0, plr->rho_lower.b = 0, plr->rho_lower.x1 = 0, plr->rho_lower.x2 = 0;
    plr->sint.x = 0, plr->sint.y = 0;
    plr->state = PLR_CONSTANTS_FIRST;
    
    plr->points = NULL;
    plr->num_points = 0;

}


void dftl_plr_add_segment(PLR *plr, Segment *seg) {
    plr->segments[plr->num_segments % plr->max_length] = *seg;
    plr->num_segments++;  
}

void dftl_plr_destroy(PLR* plr) {
    if (plr->segments != NULL) {
        free(plr->segments);
        plr->segments = NULL;
        plr->num_segments = 0;
    }
    if (plr->points != NULL) {
        free(plr->points);
        plr->points = NULL;
        plr->num_points = 0;        
    }
}

int dftl_build_segment(PLR* plr, Segment *seg) {
        if (plr->state == PLR_CONSTANTS_FIRST) { 
            seg = NULL;
            return 0;
        }
        // 建立单点段
        else if (plr->state == PLR_CONSTANTS_SECOND) {
            dftl_Segment_init(seg, 1, plr->s0.y - plr->s0.x, plr->s0.x, plr->s0.x, plr->points, plr->num_points);
        }

        // 建立多点段
        else if (plr->state == PLR_CONSTANTS_READY) {
            double avg_slope = ((plr->rho_lower.k + plr->rho_upper.k) / 2.0);
            double intercept = 0;
            if (plr->sint.x == 0 && plr->sint.y == 0) {
                // intercept = -avg_slope * plr->s0.x + plr->s0.y;    // 不精确
                avg_slope = plr->rho_lower.k;
                intercept = plr->rho_lower.b;
            } else {
                intercept = -avg_slope * plr->sint.x + plr->sint.y;
                
            }
            dftl_Segment_init(seg, avg_slope, intercept, plr->s0.x, plr->s1.x, plr->points, plr->num_points);
        }
        return 1;
}

bool dftl_should_stop(PLR* plr, Point *point) {
    if (plr->s1.x == 0) {
        if (point->x > plr->s0.x + plr->max_length || point->x <= plr->s0.x) return true;  // 段的横向长度
    }else if (point->x > plr->s1.x + plr->max_length || point->x <= plr->s1.x) return true;

    return false;
}


// up and lower 有界误差范围的线性回归学习
int dftl_process_point(PLR* plr, Point* point, Segment *seg) {

    int ret = 0;
    if (plr->state == PLR_CONSTANTS_FIRST) {
        plr->s0 = *point;
        plr->state = PLR_CONSTANTS_SECOND;
    } else if (plr->state == PLR_CONSTANTS_SECOND) {
        if (dftl_should_stop(plr, point)) {
            ret = dftl_build_segment(plr, seg);    //这里建

            plr->s0 = *point;
            plr->s1.x = 0, plr->s1.y = 0;
            plr->rho_upper.k = 0, plr->rho_upper.b = 0, plr->rho_upper.x1 = 0, plr->rho_upper.x2 = 0;
            plr->rho_lower.k = 0, plr->rho_lower.b = 0, plr->rho_lower.x1 = 0, plr->rho_lower.x2 = 0;
            plr->sint.x = 0, plr->sint.y = 0;
            plr->state = PLR_CONSTANTS_SECOND;

            free(plr->points);
            plr->points = NULL;
            plr->num_points = 0;
        }else{
            plr->s1 = *point;
            plr->rho_lower = dftl_frompoints(dftl_get_upper_bound(&plr->s0, plr->gamma), 
                                                    dftl_get_lower_bound(&plr->s1, plr->gamma));
            plr->rho_upper = dftl_frompoints(dftl_get_lower_bound(&plr->s0, plr->gamma), 
                                                    dftl_get_upper_bound(&plr->s1, plr->gamma));            
            plr->sint = dftl_inter_section(&plr->rho_upper, &plr->rho_lower);

            plr->state = PLR_CONSTANTS_READY;
        }
    } else if (plr->state == PLR_CONSTANTS_READY) {
        // 如果这个点既不在下界上面也不在上界下面而且x坐标距离已经超出了s0开始的值() 则建段
        if (dftl_is_below(point, &plr->rho_lower) || dftl_is_above(point, &plr->rho_upper) || dftl_should_stop(plr, point)) {
            ret = dftl_build_segment(plr, seg);     

            plr->s0 = *point;
            plr->s1.x = 0, plr->s1.y = 0;
            plr->rho_upper.k = 0, plr->rho_upper.b = 0, plr->rho_upper.x1 = 0, plr->rho_upper.x2 = 0;
            plr->rho_lower.k = 0, plr->rho_lower.b = 0, plr->rho_lower.x1 = 0, plr->rho_lower.x2 = 0;
            plr->sint.x = 0, plr->sint.y = 0;
            plr->state = PLR_CONSTANTS_SECOND;

            free(plr->points);
            plr->points = NULL;
            plr->num_points = 0;
        }else {
            plr->s1 = *point;

            Point s_upper = dftl_get_upper_bound(point, plr->gamma);
            Point s_lower = dftl_get_lower_bound(point, plr->gamma);

            if (dftl_is_below(&s_upper, &plr->rho_upper)) 
                plr->rho_upper = dftl_frompoints_insec(plr->sint, s_upper);
            if (dftl_is_above(&s_lower, &plr->rho_lower)) 
                plr->rho_lower = dftl_frompoints_insec(plr->sint, s_lower);
        }
    }

    if (plr->num_points == 0) {
            plr->points = (Point *)malloc(sizeof(Point));
    } else {
        plr->points = (Point *)realloc(plr->points, (plr->num_points + 1) * sizeof(Point));
    }
    plr->points[plr->num_points] = *point;
    plr->num_points++;

    return ret;
}

void dftl_plr_learn(PLR* plr, Point* points, int num_points) {
    // femu_log("[plr_learn]: in\n");

    for (int i = 0; i < num_points; i++) {
        Segment seg;
        dftl_Segment_init(&seg, 0, 0, 0, 0, NULL, 0);
        int ret = dftl_process_point(plr, &points[i], &seg);
        if (ret) {
            dftl_plr_add_segment(plr, &seg);
        }
    }

    Segment final_seg;
    dftl_Segment_init(&final_seg, 0, 0, 0, 0, NULL, 0);
    int ret = dftl_build_segment(plr, &final_seg);
   //femu_log("[process_pt]: 段学习完成, 总共点数:%d\n", plr->num_points);
    if (ret) {
           // femu_log("[plr_learn]: 学习成功后添加段(最后的段)\n");
            dftl_plr_add_segment(plr, &final_seg);
    }

    return ;
}

void dftl_GTD_update(G_map_entry *gtd, Point* points, int num_points) {
    
    dftl_plr__init(&gtd->plr);
    dftl_plr_learn(&gtd->plr, points, num_points);
    
    femu_log("[gtd_update]: plr学习完成, 索引数量%u:\n", gtd->plr.num_segments);
    if (1 && DeBUG) {
        for (int i = 0; i < gtd->plr.num_segments; i++) {
            Segment seg = gtd->plr.segments[i];
            femu_log("[gtd_update][plr_Segment %d]: k: %.2f, b: %.2f, x1: %u, x2: %u, %s, %s\n", i,
                seg.k, seg.b, seg.x1, seg.x2, seg.accurate ? "精确":"不精确", seg.consecutive ? "连续":"不连续");
        }
    }

    gtd_add_segments(gtd, gtd->plr.segments, gtd->plr.num_segments); 
    
    if (1 && DeBUG) {
        femu_log("[gtd_update]: 此次segs索引段添加到gtd完成, 此时的gtd[%lu]信息如下, num_segs: %u\n", gtd->page_idx, gtd->num_segments);
        dftl_print_gtd(gtd);
    }

    dftl_plr_destroy(&gtd->plr);    
}

int32_t dftl_binary_search(G_map_entry *gtd, Segment *seg) {
    int32_t left = 0;
    int32_t right = gtd->num_segments - 1;
    int32_t mid = 0;
    int32_t target = seg->x1;
    while (left < right) {
        mid = (left + right) >> 1;
        //femu_log("[begin] left:%u, right:%u, mid:%u\n",left, right, mid);
        if (gtd->segments[mid].x1 >= target) {
            right = mid;
        } else {
            left = mid + 1;
        }
       // fprintf("[end] left:%u, right:%u, mid:%u\n",left, right, mid);
    }
    return left;
}


void gtd_add_segment(G_map_entry *gtd, Segment *seg, int *index) {
    if (gtd->num_segments == 0) {
        gtd->segments = (Segment *)malloc(1 * sizeof(Segment));
        gtd->segments[0] = *seg;  // C语言结构浅拷贝,有可能出bug 
        gtd->num_segments++; 
        *index = 0;
    } else {
        gtd->segments = (Segment *)realloc(gtd->segments, (gtd->num_segments + 1) * sizeof(Segment));
        // 二分查找插入位置,永远是右半区的最左边界
        int32_t pos = dftl_binary_search(gtd, seg);
        if (pos == gtd->num_segments - 1 && seg->x1 > gtd->segments[pos].x1) {
            pos = gtd->num_segments;
            *index = pos;
            gtd->num_segments++;
            gtd->segments[pos] = *seg;
        } else {
            *index = pos;
            gtd->num_segments++;
            for (int i = gtd->num_segments - 1; i >= pos + 1; i--) {
                gtd->segments[i] = gtd->segments[i - 1];
            }
            gtd->segments[pos] = *seg;          
        }
    }
    
    // update Bitmap
    uint32_t k = 0;
    uint32_t i = seg->x1;
    while (i <= seg->x2 && k < seg->filter.length) {
        gtd->bitmap_filter.bitmap[i++] = seg->filter.bitmap[k++];
    }  
}

void dftl_Segs_add_segment(Segs *segs, Segment *seg, int seg_id) {
    if (segs->num_segments == 0) {
        segs->segments = (Segment *)malloc(1 * sizeof(Segment));
        segs->segments[0] = *seg;  
        segs->segment_id[0] = seg_id; 
        segs->num_segments++; 
    } else {
        segs->segments = (Segment *)realloc(segs->segments, (segs->num_segments + 1) * sizeof(Segment));
        segs->segments[segs->num_segments] = *seg;
        segs->segment_id[segs->num_segments] = seg_id;
        segs->num_segments++;
    }
}

void gtd_del_segment(G_map_entry *gtd, int pos) {
    if (pos < 0 || pos >= gtd->num_segments) {
        femu_log("删除位置不对\n");
        return ;
    }
    for (int i = pos; i < gtd->num_segments - 1; i++) {
        gtd->segments[i] = gtd->segments[i + 1];
    }
    gtd->num_segments--;
    gtd->segments = (Segment *)realloc(gtd->segments, gtd->num_segments * sizeof(Segment));

    if (gtd->num_segments == 0) {
        free(gtd->segments);
        gtd->segments = NULL;
    }
}

void gtd_add_segments(G_map_entry *gtd, Segment *segments, int num_segments) {

    for (int i = 0; i < num_segments; i++) {
        int index = 0;

        if (gtd->num_segments == 0) {
            gtd_add_segment(gtd, &segments[i], &index);
            continue;
        }
        Segs overlap_segs;
        overlap_segs.num_segments = 0;
        overlap_segs.segments = NULL;
        memset(overlap_segs.segment_id, 0, Write_Buffer_Entries * sizeof(int));

        gtd_add_segment(gtd, &segments[i], &index);
        if (index != 0) {
            if (gtd->segments[index].x1 <= gtd->segments[index - 1].x2) {
                dftl_Segs_add_segment(&overlap_segs, &gtd->segments[index - 1], index - 1);
            }
        }
        for (int j = index + 1; j < gtd->num_segments; j++) {
            if (gtd->segments[j].x1 > segments[i].x2) {
                break;
            }
            dftl_Segs_add_segment(&overlap_segs, &gtd->segments[j], j);            
        }
        

        uint32_t indices_to_delete[ENTRY_PER_PAGE];
        int      indect_pointer = 0;

        for (int j = 0; j < overlap_segs.num_segments; j++) {
            int same_level = 0;
            dftl_Segment_merge(&segments[i], &overlap_segs.segments[j], &same_level);

            if (same_level == -1) {
                indices_to_delete[indect_pointer++] = overlap_segs.segment_id[j]; 
            } else if (same_level == 0) {
                //TODO               
            } else if (same_level == 1) {
                // 合并后不重叠
                int pos = overlap_segs.segment_id[j];
                gtd->segments[pos] = overlap_segs.segments[j];
            }

        }
    
        if (overlap_segs.segments != NULL) {
            free(overlap_segs.segments);
            overlap_segs.num_segments = 0;
            overlap_segs.segments = NULL;
        }

        for (int j = indect_pointer - 1; j >= 0; j--) {
            gtd_del_segment(gtd, indices_to_delete[j]);
        }

    }
    return ;
}

void gtd_entry_init(G_map_entry *gtd, uint64_t page_idx) {
    gtd->page_idx = page_idx;
    gtd->ppa      = UNMAPPED_PPA;
    
    // 精确段 误差0
    dftl_plr_init(&gtd->plr, PLR_Error_bound);

    gtd->bitmap_filter.length = ENTRY_PER_PAGE + 1;
    gtd->bitmap_filter.bitmap = (unsigned char *)malloc(gtd->bitmap_filter.length * sizeof(unsigned char));
    memset(gtd->bitmap_filter.bitmap, 0, gtd->bitmap_filter.length * sizeof(unsigned char));

    gtd->num_segments = 0;
    gtd->segments = NULL;
}

uint64_t gtd_lookup(G_map_entry *gtd, uint64_t lpn) {
    uint32_t page_idx = lpn / ENTRY_PER_PAGE;
    uint8_t find_lpa = lpn % ENTRY_PER_PAGE;   
    uint64_t ppa = 0;

    // check bit_map first;
    if (gtd->bitmap_filter.bitmap[find_lpa]) {
        Segment find_seg;
        dftl_Segment_init(&find_seg, 0, 0, find_lpa, find_lpa, NULL, 0);
        int pos = 0;
        // 第一个大于等于tar的位置
        pos = dftl_binary_search(gtd, &find_seg);
        if (pos == 0 || (pos < gtd->num_segments && gtd->segments[pos].x1 == find_lpa))
            pos = pos;
        else 
            pos -= 1;
        find_seg = gtd->segments[pos];

        // need verify
        ppa = dftl_Segment_gety(&find_seg, false, find_lpa);
    }

    if (ppa == 0) {
        femu_log("gtd[%lu]上未找到对于 lpn: %lu的写入信息\n", gtd->page_idx, lpn);
    }

    return ppa;    
}


void dftl_print_gtd(G_map_entry *gtd) {
    int valid_bits = 0;
    for (int i = 0; i < gtd->bitmap_filter.length; i++) {
        valid_bits += gtd->bitmap_filter.bitmap[i];
    }
    femu_log("[gtd %lu]: valid_bits: %d\n", gtd->page_idx, valid_bits);

    for (int i = 0; i < gtd->num_segments; i++) {
        Segment seg = gtd->segments[i];
        femu_log("[gtd_Segment %d]: k: %.2f, b: %.2f, x1: %u, x2: %u, %s, %s\n", i,
                seg.k, seg.b, seg.x1, seg.x2, seg.accurate ? "精确":"不精确", seg.consecutive ? "连续":"不连续");
    }
}

void dftl_static(DFTLTable *d_maptbl) {
    femu_log("[DFTL]write_cnt:%d, read_cnt: %d, read_miss: %d, cmt_hit: %d, gtd_hit: %d, double_read: %d, evit_cnt: %d\n",
    d_maptbl->counter.group_write_cnt, d_maptbl->counter.group_read_cnt,
    d_maptbl->counter.group_read_miss, 
    d_maptbl->counter.group_cmt_hit, d_maptbl->counter.GTD_hit,
    d_maptbl->counter.group_double_read, d_maptbl->counter.evit_cnt);
    femu_log("[DFTL]cache_cnt: %lu/%lu, nand_cnt: %lu/%lu\n", d_maptbl->CMT->count, d_maptbl->CMT->capacity, 
    d_maptbl->nand_cache->count, d_maptbl->nand_cache->capacity);

    uint64_t num_segs = 0;
    for (int i = 0; i < d_maptbl->nand_cache->capacity; i++) {
        num_segs += d_maptbl->GTD[i].num_segments;
    }
    femu_log("[DFTL]GTD_seg_cnt: %lu\n", num_segs);
}


void o_static(struct ssd* ssd) {
    femu_log("[Normal_FTL]write_cnt:%d, map_cnt: %d, read_cnt: %d, read_miss: %d",
                ssd->write_cnt, ssd->map_cnt, ssd->read_cnt, ssd->read_miss);
    return ;
}