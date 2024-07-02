#include "leaftl.h"
/**
 * @brief CRB_method
 * 
 * 还未用上 
 */

void crb_init(CRB *crb) {
    crb->num_segments = 0;
}

void crb_insert_segment(CRB *crb, Segment* seg) {
    uint8_t start_lpa = seg->x1;
    // Check if the CRB is already full
    if (crb->num_segments == MAX_CRB_SEGMENTS - 1) {
        femu_log("CRB is full. Cannot insert segment.\n");
        return;
    }

    // Check if the start_lpa already exists in the CRB
    for (unsigned char i = 0; i < crb->num_segments; i++) {
        if (crb->segments[i].x1 == start_lpa) {
            femu_log("Segment with start LPA %u already exists. Updating start LPA.\n", start_lpa);
            // 这一步应该是更新成往右相邻的点的LPA, 但是Segment的信息只有s,l,k,b并不知道下个点是啥有点迷
            crb->segments[i].x1++;  
            return;
        }
    }

    // Insert the new segment
    crb->segments[crb->num_segments].x1 = start_lpa;
    crb->segments[crb->num_segments].x2 = seg->x2;
    crb->num_segments++;

    // Sort the segments by their starting LPA using bubble sort
    // 冒泡排序
    for (unsigned char i = 0; i < crb->num_segments - 1; i++) {
        for (unsigned char j = 0; j < crb->num_segments - i - 1; j++) {
            if (crb->segments[j].x1 > crb->segments[j + 1].x1) {
                // Swap the segments
                Segment temp = crb->segments[j];
                crb->segments[j] = crb->segments[j + 1];
                crb->segments[j + 1] = temp;
            }
        }
    }
}

// 找到属于哪个近似段或者是该近似段的起始LPA
unsigned int crb_find_segment(CRB *crb, uint8_t lpa) {
    // Binary search to find the segment containing the given LPA
    unsigned int left = 0;
    unsigned int right = crb->num_segments - 1;

    while (left <= right) {
        unsigned int mid = left + (right - left) / 2;

        if (crb->segments[mid].x1 <= lpa && lpa < crb->segments[mid].x2) {
            // Found the segment containing the LPA
            // return crb->segments[mid].x1;
            return mid;
        }
        else if (crb->segments[mid].x1 > lpa) {
            // LPA is in the left half
            right = mid - 1;
        }
        else {
            // LPA is in the right half
            left = mid + 1;
        }
    }

    // LPA is not found in any segment
    return -1;
}

/**
 * @brief Segment_method
 * 
 */


// 可能出现bug
int Segment_is_valid(Segment *seg, uint32_t x) {
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
    // 
    // 添加位图判断 判断 lpn 是否在段上存在
    else {
        if (seg->filter.length > 0 && (x - seg->x1) <= seg->filter.length) 
            return seg->filter.bitmap[x - seg->x1];
    }
    return 1;
}



uint32_t Segment_gety(Segment *seg, bool check, uint32_t x) {
    int predict = 0;
    if (!check || Segment_is_valid(seg, x)) {
            predict = (int)(x*seg->k + seg->b);
            return predict;
    }
    return predict;
}

uint32_t Segment_gety_upper(Segment *seg, bool check, uint32_t x) {
    int predict = 0;
    if (!check || Segment_is_valid(seg, x)) {
            int pos = seg->filter.bitmap_upper[x - seg->x1];
            if (pos == 1) predict = (int)(x*seg->k + seg->b);
            else predict = round(x*seg->k + seg->b);
            return predict;
    }
    return predict;
}

int isConsecutive(Point *points, int num_points) {
    int i;
    int diff = points[1].x - points[0].x;  // 计算第一个相邻点的 x 坐标差异
    for (i = 2; i < num_points; i++) {
        if (points[i].x - points[i-1].x != diff) {
            return 0;  // 存在不连续的 x 坐标差异
        }
    }

    return 1;  // 所有的 x 坐标差异都相等，即连续
}


// 判断段属性: 精确性，连续性
void Segment_check_properties(Segment *seg, Point *points, int num_points) {
    seg->accurate = true, seg->consecutive = true;
    
    for (int i = 0; i < num_points; i++) {
        uint32_t ppa;
        ppa = Segment_gety(seg, false, points[i].x);
        if (ppa != points[i].y) {
            seg->accurate = false;
        }
        if (!seg->accurate) return ;
    }
    if (!isConsecutive(points, num_points)) {
            seg->consecutive = false;
    }
}


void Segment_init(Segment *seg, double k, double b, int x1, int x2, Point *points, int num_points) {
    seg->k = k;
    seg->b = b;
    seg->x1 = x1;
    seg->x2 = x2;
    seg->accurate = true;
      
    seg->filter.length = 0;
    seg->filter.bitmap = NULL;


    if (points != NULL) {
        Segment_check_properties(seg, points, num_points);

        if (1 || !seg->consecutive) {
            // 初始化位图
            if (seg->filter.bitmap != NULL) {
                free(seg->filter.bitmap);
                seg->filter.length = 0;
                seg->filter.bitmap = NULL;
            }
            seg->filter.length = seg->x2 - seg->x1 + 1;
            seg->filter.bitmap = (unsigned char *)malloc(seg->filter.length * sizeof(unsigned char));
            seg->filter.bitmap_upper = (unsigned char *)malloc(seg->filter.length * sizeof(unsigned char));
            memset(seg->filter.bitmap, 0, seg->filter.length * sizeof(unsigned char));
            memset(seg->filter.bitmap_upper, 0, seg->filter.length * sizeof(unsigned char));

            if (seg->filter.bitmap != NULL) {
                memset(seg->filter.bitmap, 0, seg->filter.length * sizeof(unsigned char));
                for (int i = 0; i < num_points; i++) {
                    seg->filter.bitmap[points[i].x - seg->x1] = 1;
                } 
            }
            if (seg->filter.bitmap_upper != NULL) {
                memset(seg->filter.bitmap_upper, 0, seg->filter.length * sizeof(unsigned char));
                for (int i = 0; i < num_points; i++) {
                    uint32_t ppa;
                    ppa = Segment_gety(seg, false, points[i].x);
                    if (ppa >= points[i].y) {
                        seg->filter.bitmap_upper[points[i].x - seg->x1] = 1;
                    } else {
                        seg->filter.bitmap_upper[points[i].x - seg->x1] = 2;
                    } 
                }                
            }
        }
    }
}

bool Segment_overlaps(Segment *seg1, Segment *seg2) {
    return L_MIN(seg1->x2, seg2->x2) - L_MAX(seg1->x1, seg2->x1) >= 0;
}


// may has a bug // fix
void Segment_merge(Segment *new_seg, Segment *old_seg, int *samelevel) {
    // if (!old_seg || !old_seg->k || !old_seg->x1 || !old_seg->filter.bitmap) {
    //     *samelevel = 2;
    //     return;
    // }

    //else if (new_seg->consecutive)  no use new_seg->consecutive
    if (0 && DeBUG) femu_log("new_seg.x1: %u, x2: %u, old_seg.x1: %u, old_seg.x2: %u\n", 
                                    new_seg->x1, new_seg->x2, old_seg->x1, old_seg->x2);

    int start = L_MIN(new_seg->x1, old_seg->x1);
    int end   = L_MAX(new_seg->x2, old_seg->x2);

    Bitmap new_bm, old_bm;
    new_bm.length = end - start + 1;
    old_bm.length = end - start + 1;
    new_bm.bitmap = (unsigned char *)malloc(new_bm.length * sizeof(unsigned char));
    old_bm.bitmap = (unsigned char *)malloc(old_bm.length * sizeof(unsigned char));
    memset(new_bm.bitmap, 0, new_bm.length * sizeof(unsigned char)); 
    memset(old_bm.bitmap, 0, old_bm.length * sizeof(unsigned char));

    if (0 && DeBUG) {
        femu_log("new_seg.bitmap:");
        for (int i = 0; i < new_seg->filter.length; i++) {
            printf("%u", new_seg->filter.bitmap[i]);
        }
        printf("\n");
        femu_log("old_seg.bitmap:");
        for (int i = 0; i < old_seg->filter.length; i++) {
            printf("%u", old_seg->filter.bitmap[i]);
        }
        printf("\nnew_bm.length: %d, new_bm.bitmap:", new_bm.length);
        for (int i = 0; i < new_bm.length; i++) {
            printf("%u", new_bm.bitmap[i]);
        }    
        printf("\n");
    }

    int j = 0;
    for (int i = new_seg->x1 - start; i < new_seg->x2 - start + 1; i++) {
        new_bm.bitmap[i] = new_seg->filter.bitmap[j];
        j++;
    }
    j = 0;
    for (int i = old_seg->x1 - start; i < old_seg->x2 - start + 1; i++) {
        old_bm.bitmap[i] = old_seg->filter.bitmap[j];
        j++;
    }
    j = 0;
 
    for (int i = 0; i < end - start + 1; i++) {
        old_bm.bitmap[i] = old_bm.bitmap[i] & (~new_bm.bitmap[i]);
    }
    
    // Debug
    if (0 && DeBUG) {
        printf("new_bm:");
        for (int i = 0; i < end - start + 1; i++) {
            printf(" %u", new_bm.bitmap[i]);
        }
        printf("\n");
        printf("old_bm:");
        for (int i = 0; i < end - start + 1; i++) {
            printf(" %u", old_bm.bitmap[i]);
        }
        printf("\n");
    }

    int first_valid = -1;
    for (int i = 0; i < end - start + 1; i++) {
        if (old_bm.bitmap[i]) {
            first_valid = i;
            break;
        }
    }
    if (first_valid == -1) {
        *samelevel = -1;
        return ;
    }
    int last_valid =  -1;
    for (int i = end - start; i >= 0; i--) {
        if (old_bm.bitmap[i]) {
            last_valid = i;
            break;
        }
    }
    old_seg->x1 = first_valid + start;
    old_seg->x2 = last_valid  + start;  // 在这里对old_seg处理,只是在overlap_segs中的old_seg做处理,
                                        // 并没有对原来的logplr中的segs中的old_seg做处理


    for (int i = first_valid; i < last_valid + 1; i++) {
        old_seg->filter.bitmap[j] =  old_bm.bitmap[i];
        j++;
    }

    if (new_seg->x1 > old_seg->x2 || new_seg->x2 < old_seg->x1) {
        *samelevel = 1;  // 存在不重叠
        return ;
    }     

    *samelevel = 0;   // 存在仍重叠 
    return ;

}



void SimpleSegment_init(SimpleSegment *simpseg, double k, double b, int x1, int x2) {
    simpseg->k = k;
    simpseg->b = b;
    simpseg->x1 = x1;
    simpseg->x2 = x2;    
}


int get_y(SimpleSegment *simpleseg, int x) {
    int predict;
    predict = round(x * simpleseg->k + simpleseg->b);
    return predict;
}


InsecPoint inter_section(SimpleSegment *s1, SimpleSegment *s2) {
    InsecPoint insec_pt;
    insec_pt.x = 0, insec_pt.y = 0;
    // fix bug 如果两条直线平行或重叠是不会有交点的
    if (s1->k != s2->k) {
        insec_pt.x = (double) ((s2->b - s1->b) / (s1->k - s2->k));
        insec_pt.y = (double) ((s1->k * s2->b - s2->k * s1->b) / (s1->k - s2->k));
    }
    return insec_pt;
}

bool is_above(Point *pt, SimpleSegment *s) {
    return pt->y > (int)(s->k * pt->x + s->b);
}

bool is_below(Point *pt, SimpleSegment *s) {
    return pt->y < (int)(s->k * pt->x + s->b);
}

Point get_upper_bound(Point *pt, double gamma) {
    Point p;
    p.x = pt->x, p.y = pt->y + gamma;
    return p;
}

Point get_lower_bound(Point *pt, double gamma) {
    Point p;
    p.x = pt->x, p.y = pt->y - gamma;
    return p;
}

SimpleSegment frompoints(Point p1, Point p2) {
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

SimpleSegment frompoints_insec(InsecPoint p1, Point p2) {
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

// typedef struct {
//     const char* FIRST;
//     const char* SECOND;
//     const char* READY;
// } PLR_CONSTANTS;



/**
 * @brief PLR_method
 * 
 */
void plr_init(PLR* plr, double gamma) {
    plr->gamma = gamma;
    plr->max_length = Write_Buffer_Entries;

    plr__init(plr);
}


// Tem_struct, use for build segments, and then destry after insert LSM struc
void plr__init(PLR* plr) {

    // plr_destroy(plr);

    plr->segments = g_malloc0(sizeof(Segment) * (plr->max_length));
    for (int i = 0; i < plr->max_length; i++) {
        Segment_init(&plr->segments[i], 0, 0, 0, 0, NULL, 0);
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


void plr_add_segment(PLR *plr, Segment *seg) {
    // if (plr->num_segments == 0) {
    //     plr->segments = (Segment *)malloc(1 * sizeof(Segment));
    //     plr->segments[0] = *seg;  // C语言结构浅拷贝,有可能出bug
    //     plr->num_segments++; 
    // } else {
    //     plr->segments = (Segment *)realloc(plr->segments, (plr->num_segments + 1) * sizeof(Segment));
    //     plr->segments[plr->num_segments] = *seg;
    //     plr->num_segments++;
    // }
    plr->segments[plr->num_segments % plr->max_length] = *seg;
    plr->num_segments++;  
}

void plr_destroy(PLR* plr) {
   // femu_log("[plr_destroy]: 此次学习完成, 删除临时plr->segments\n");
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
   // free(plr);
}

int build_segment(PLR* plr, Segment *seg) {
        if (plr->state == PLR_CONSTANTS_FIRST) { 
            seg = NULL;
            return 0;
        }
        // 建立单点段
        else if (plr->state == PLR_CONSTANTS_SECOND) {
            Segment_init(seg, 1, plr->s0.y - plr->s0.x, plr->s0.x, plr->s0.x, plr->points, plr->num_points);
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
            if (avg_slope < 0) {
                femu_log("seg.k: %.2f, s0.x: %u, s0.y: %u, s1.x: %u, s1.y: %u\n", avg_slope, 
                    plr->s0.x, plr->s0.y, plr->s1.x, plr->s1.y);
            }
            Segment_init(seg, avg_slope, intercept, plr->s0.x, plr->s1.x, plr->points, plr->num_points);
        }
        return 1;
}

bool should_stop(PLR* plr, Point *point) {
    if (plr->s1.x == 0) {
        if (point->x > plr->s0.x + plr->max_length || point->x <= plr->s0.x) return true;  // 段的横向长度
    }else if (point->x > plr->s1.x + plr->max_length || point->x <= plr->s1.x) return true;

    return false;
}


// up and lower 有界误差范围的线性回归学习
int process_point(PLR* plr, Point* point, Segment *seg) {

   // femu_log("[process_point]: in\n");

    int ret = 0;
    if (plr->state == PLR_CONSTANTS_FIRST) {
        plr->s0 = *point;
        plr->state = PLR_CONSTANTS_SECOND;
    } else if (plr->state == PLR_CONSTANTS_SECOND) {
        if (should_stop(plr, point)) {
            ret = build_segment(plr, seg);    //这里建

            plr->s0 = *point;
            plr->s1.x = 0, plr->s1.y = 0;
            plr->rho_upper.k = 0, plr->rho_upper.b = 0, plr->rho_upper.x1 = 0, plr->rho_upper.x2 = 0;
            plr->rho_lower.k = 0, plr->rho_lower.b = 0, plr->rho_lower.x1 = 0, plr->rho_lower.x2 = 0;
            plr->sint.x = 0, plr->sint.y = 0;
            plr->state = PLR_CONSTANTS_SECOND;

            //femu_log("[process_pt]: 段学习完成, 总共点数:%d\n", plr->num_points);
            free(plr->points);
            plr->points = NULL;
            plr->num_points = 0;
        }else{
            plr->s1 = *point;
            // rho_lower 和 rho_upper 和 sint  与gamma的斜率相关
            plr->rho_lower = frompoints(get_upper_bound(&plr->s0, plr->gamma), 
                                                    get_lower_bound(&plr->s1, plr->gamma));
            plr->rho_upper = frompoints(get_lower_bound(&plr->s0, plr->gamma), 
                                                    get_upper_bound(&plr->s1, plr->gamma));            
            plr->sint = inter_section(&plr->rho_upper, &plr->rho_lower);

            plr->state = PLR_CONSTANTS_READY;
        }
    } else if (plr->state == PLR_CONSTANTS_READY) {
        // 如果这个点既不在下界上面也不在上界下面而且x坐标距离已经超出了s0开始的值(256) 则建段
        if (is_below(point, &plr->rho_lower) || is_above(point, &plr->rho_upper) || should_stop(plr, point)) {
            ret = build_segment(plr, seg);     

            plr->s0 = *point;
            plr->s1.x = 0, plr->s1.y = 0;
            plr->rho_upper.k = 0, plr->rho_upper.b = 0, plr->rho_upper.x1 = 0, plr->rho_upper.x2 = 0;
            plr->rho_lower.k = 0, plr->rho_lower.b = 0, plr->rho_lower.x1 = 0, plr->rho_lower.x2 = 0;
            plr->sint.x = 0, plr->sint.y = 0;
            plr->state = PLR_CONSTANTS_SECOND;

            //femu_log("[process_pt]: 段学习完成, 总共点数:%d\n", plr->num_points);
            free(plr->points);
            plr->points = NULL;
            plr->num_points = 0;
        }else {
            // QS
            plr->s1 = *point;

            Point s_upper = get_upper_bound(point, plr->gamma);
            Point s_lower = get_lower_bound(point, plr->gamma);

            if (is_below(&s_upper, &plr->rho_upper)) 
                plr->rho_upper = frompoints_insec(plr->sint, s_upper);
            if (is_above(&s_lower, &plr->rho_lower)) 
                plr->rho_lower = frompoints_insec(plr->sint, s_lower);
        }
    }

    if (plr->num_points == 0) {
            plr->points = (Point *)malloc(sizeof(Point));
    } else {
        plr->points = (Point *)realloc(plr->points, (plr->num_points + 1) * sizeof(Point));
    }
    plr->points[plr->num_points] = *point;
    plr->num_points++;

  //  femu_log("[process_point]: out\n");
    return ret;
}

void plr_learn(PLR* plr, Point* points, int num_points) {
    // femu_log("[plr_learn]: in\n");

    for (int i = 0; i < num_points; i++) {
        Segment seg;
        Segment_init(&seg, 0, 0, 0, 0, NULL, 0);
        int ret = process_point(plr, &points[i], &seg);
        if (ret) {
            //femu_log("[plr_learn]: 学习成功后添加段\n");
            plr_add_segment(plr, &seg);
        }
    }

    Segment final_seg;
    Segment_init(&final_seg, 0, 0, 0, 0, NULL, 0);
    int ret = build_segment(plr, &final_seg);
   //femu_log("[process_pt]: 段学习完成, 总共点数:%d\n", plr->num_points);
    if (ret) {
           // femu_log("[plr_learn]: 学习成功后添加段(最后的段)\n");
            plr_add_segment(plr, &final_seg);
    }

    // femu_log("[plr_learn]: out\n");
   // *num_segments = plr->num_segments;
    return ;
}

// TODO
void plr_sorted(PLR* plr) {

}


/**
 * @brief LogPLR_method
 * 
 * 有关日志结构索引段的方法函数
 * 
 * @param logplr 输入日志结构学习索引段
 * @param other 
 * 
 * @return void
 * 
 * @details 
 * LogPLR_init();
 * binary_search();
 * LogPLR_add_segment();
 * LogPLR_add_segment
 * LogPLR_del_segment
 * Segs_add_segment();   // Segs's method
 */



void LogPLR_init(LogPLR *logplr, int level) {
    logplr->level = level;
    logplr->segments = NULL;
    logplr->num_segments = 0;
    logplr->max_length = 256;
}


// fix: 寻找大于等于target的右半区的左边界点
int32_t binary_search(LogPLR *logplr, Segment *seg) {
    int32_t left = 0;
    int32_t right = logplr->num_segments - 1;
    int32_t mid = 0;
    int32_t target = seg->x1;
    while (left < right) {
        mid = (left + right) >> 1;
        //femu_log("[begin] left:%u, right:%u, mid:%u\n",left, right, mid);
        if (logplr->segments[mid].x1 >= target) {
            right = mid;
        } else {
            left = mid + 1;
        }
       // fprintf("[end] left:%u, right:%u, mid:%u\n",left, right, mid);
    }
    return left;
}

// L_i添加段
// 如果释放掉PLR之后，LogPLR中的段不会不见 C语言结构的浅拷贝
void LogPLR_add_segment(LogPLR *logplr, Segment *seg, int *index) {
    if (logplr->num_segments == 0) {
        logplr->segments = (Segment *)malloc(1 * sizeof(Segment));
        logplr->segments[0] = *seg;  // C语言结构浅拷贝,有可能出bug 
        logplr->num_segments++; 
        *index = 0;
    } else {
        logplr->segments = (Segment *)realloc(logplr->segments, (logplr->num_segments + 1) * sizeof(Segment));
        // 二分查找插入位置,永远是右半区的最左边界
        int32_t pos = binary_search(logplr, seg);
        if (pos == logplr->num_segments - 1 && seg->x1 > logplr->segments[pos].x1) {
            pos = logplr->num_segments;
            *index = pos;
            logplr->num_segments++;
            logplr->segments[pos] = *seg;
        } else {
            *index = pos;
            logplr->num_segments++;
            for (int i = logplr->num_segments - 1; i >= pos + 1; i--) {
                logplr->segments[i] = logplr->segments[i - 1];
            }
            logplr->segments[pos] = *seg;          
        }   
    }   
}

void LogPLR_del_segment(LogPLR *logplr, int pos) {
    if (pos < 0 || pos >= logplr->num_segments) {
        femu_log("删除位置不对\n");
        return ;
    }
    for (int i = pos; i < logplr->num_segments - 1; i++) {
        logplr->segments[i] = logplr->segments[i + 1];
    }
    logplr->num_segments--;
    logplr->segments = (Segment *)realloc(logplr->segments, logplr->num_segments * sizeof(Segment));

    if (logplr->num_segments == 0) {
        free(logplr->segments);
        logplr->segments = NULL;
    }
    //femu_log("删除成功: %d\n", pos);
}


void Segs_add_segment(Segs *segs, Segment *seg, int seg_id) {
    if (segs->num_segments == 0) {
        segs->segments = (Segment *)malloc(1 * sizeof(Segment));
        segs->segments[0] = *seg;  // C语言结构浅拷贝,有可能出bug
        segs->segment_id[0] = seg_id; 
        segs->num_segments++; 
    } else {
        segs->segments = (Segment *)realloc(segs->segments, (segs->num_segments + 1) * sizeof(Segment));
        segs->segments[segs->num_segments] = *seg;
        segs->segment_id[segs->num_segments] = seg_id;
        segs->num_segments++;
    }
}




/**
 * @brief Group_method
 * 
 */


void Group_lookup(void) {
    
}

void Group_seg_merge(void) {
    
}

void Group_gc(void) {

}

void Group_init(Group *group, double gamma, int group_id) {
    plr_init(&group->plr, gamma);
    group->L = NULL;
    group->num_levels = 0;
    group->max_levels = 256;

    group->group_id = group_id;
}

void Group_add_LogPLR(Group *group) {
    if (group->num_levels == 0) {
        group->L = (LogPLR *)malloc(1 * sizeof(LogPLR));
        LogPLR_init(&group->L[group->num_levels], group->num_levels);
        group->num_levels++;
    } else {
        group->L = (LogPLR *)realloc(group->L, (group->num_levels + 1) * sizeof(LogPLR));
        LogPLR_init(&group->L[group->num_levels], group->num_levels);
        group->num_levels++;        
    }
}

void Group_del_segments(Group *group, int level, int pos) {
    LogPLR_del_segment(&group->L[level], pos);
    return ;
}

void Group_add_segments(Group *group, int level, Segment *segments, int num_segments, bool recursive) {

    //femu_log("[Group_add_segments]: in\n");
    if (group->num_levels >= group->max_levels) {
        // TODO
        // Compact合并
        //Group_compact(group);
        femu_log("[Group_add_segments]:group[%d]->num_levels: %d, overflow!!!\n", group->group_id, group->num_levels);
    }

    while (group->num_levels <= level) {
        Group_add_LogPLR(group);
    }

   // LogPLR *Li = &group->L[level];
    Segs_Simple confict_segs;
    confict_segs.num_segments = 0;
    confict_segs.segments = NULL;

    // femu_log("[Group_add_segments]: num_segments:%d, %s递归添加\n", num_segments, recursive ? "是" : "不是");
    for (int i = 0; i < num_segments; i++) {
        //femu_log("segment:[%d].x1 = %u, add to level: L%d\n", i,segments[i].x1, level);
        int index = 0;

        // debug
        if (group->L[level].num_segments >= group->L[level].max_length) {
            // TODO 替换策略 继续添加统计段
            femu_log("[Group_add_segments]:g[%d]->L[%d].num_segs: %d, error 超出最大值! \n", group->group_id, level, group->L[level].num_segments);

            if (group->L[level].num_segments >= 1024) {
                Print_Group(group);
                break;
            }
            //break;
        }

        if (group->L[level].num_segments == 0) {
            LogPLR_add_segment(&group->L[level], &segments[i], &index);
            continue;
        }
        Segs overlap_segs;
        overlap_segs.num_segments = 0;
        overlap_segs.segments = NULL;
        memset(overlap_segs.segment_id, 0, Write_Buffer_Entries * sizeof(int));

        LogPLR_add_segment(&group->L[level], &segments[i], &index);
        if (index != 0) {
            // 在前面的最多只有一个重叠
            // Segs_add_segment 处理重复段函数
            if (group->L[level].segments[index].x1 <= group->L[level].segments[index - 1].x2) {
                Segs_add_segment(&overlap_segs, &group->L[level].segments[index - 1], index - 1);
            }
        }

        for (int j = index + 1; j < group->L[level].num_segments; j++) {
            // 添加后续的重叠区间段
            if (group->L[level].segments[j].x1 > segments[i].x2) {
                break;
            }
            Segs_add_segment(&overlap_segs, &group->L[level].segments[j], j);            
        }
        
        // debug
        if (0 && DeBUG) {
            femu_log("group[%d]->L[%d] 目前重叠段数为: %d\n", group->group_id, level, overlap_segs.num_segments);
        }

        uint8_t indices_to_delete[Write_Buffer_Entries];
        int     indect_pointer = 0;

        for (int j = 0; j < overlap_segs.num_segments; j++) {
            int same_level = 2;
            // 段合并 same_level 判断合并操作后的旧段状态 2:原本就不重叠, 1: 处理后同一层不重叠, 0:处理后同一层仍重叠, -1:处理后旧段不存在，被新段完全包含 
            Segment_merge(&segments[i], &overlap_segs.segments[j], &same_level);
            // femu_log("same_level: %d\n", same_level);
            if (same_level == -1) {
                indices_to_delete[indect_pointer++] = overlap_segs.segment_id[j]; 
            } else if (same_level == 0) {
                
                if (confict_segs.num_segments == 0) {
                    confict_segs.segments = (Segment *)malloc(1 * sizeof(Segment));
                    confict_segs.segments[0] = overlap_segs.segments[j];
                    confict_segs.num_segments++;
                } else {
                    confict_segs.segments = (Segment *)realloc(confict_segs.segments, (confict_segs.num_segments + 1) * sizeof(Segment));
                    confict_segs.segments[confict_segs.num_segments] = overlap_segs.segments[j];
                    confict_segs.num_segments++;
                }

                indices_to_delete[indect_pointer++] = overlap_segs.segment_id[j];                   
            } else if (same_level == 1) {
                // 合并后不重叠此时应该处理原group->L[level].segments[j]的段
                // 将overlap_segs.segments[j]段的值同步到group->L[level].segments[j]的段的值中
                int pos = overlap_segs.segment_id[j];
                group->L[level].segments[pos] = overlap_segs.segments[j];

            } else if (same_level == 2) {
                // 基本不会进入此个条件
            }

        }
    
        if (overlap_segs.segments != NULL) {
            free(overlap_segs.segments);
            overlap_segs.num_segments = 0;
            overlap_segs.segments = NULL;
        }

        // 新段插入时overlaps重叠矛盾时立即处理, 状态为 仍存在重叠0 或者 旧段被新段完全覆盖-1, 的这些段需要删除
        // 从尾往后删除
        for (int j = indect_pointer - 1; j >= 0; j--) {
            Group_del_segments(group, level, indices_to_delete[j]);
        }

    }

    //  如果添加段时出现区间重叠，处理后准太则将旧段往下层推，如果在下层还是重叠，则为该段创建新的Level
    if (recursive) {
        if (confict_segs.num_segments != 0) {
            //femu_log("将 %d 个重叠区间添加到 L[%d]\n", confict_segs.num_segments, level + 1);
            Group_add_segments(group, level + 1, confict_segs.segments, confict_segs.num_segments, false);
        }
    } else {
        // 如果递归添加后仍重叠 则直接为这些重叠段建立新的L[level]索引水平L_i
        if (confict_segs.num_segments != 0) {
            Group_add_LogPLR(group);
            for (int i = 0; i < confict_segs.num_segments; i++) {
                int no_use = 0;
                LogPLR_add_segment(&group->L[group->num_levels - 1], &confict_segs.segments[i], &no_use); // 表示用不上id                
            }       
        }
    }
    if (confict_segs.segments != NULL) {
        free(confict_segs.segments);
        confict_segs.num_segments = 0;
        confict_segs.segments = NULL;
    }
    // femu_log("[Group_add_segments]: out\n");
    return ;
}

void Group_update(Group *group, Point* points, int num_points) {
    // Points 应该是已经排序好的 lpa : ppa
    u_int64_t stime = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
    plr__init(&group->plr);
    plr_learn(&group->plr, points, num_points);
    plr_sorted(&group->plr);
    u_int64_t etime = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
    
   // femu_log("[Group_update]: plr学习完成, num_segments:%d, learn_time: %lu ns\n", group->plr.num_segments, etime - stime);
    
    // print
    if (0 && DeBUG) {
        for (int i = 0; i < group->plr.num_segments; i++) {
            Segment seg = group->plr.segments[i];
            femu_log("[Group_update][plr_Segment %d]: k: %.2f, b: %.2f, x1: %u, x2: %u, %s, %s\n", i,
                seg.k, seg.b, seg.x1, seg.x2, seg.accurate ? "精确":"不精确", seg.consecutive ? "连续":"不连续");
        }
    }
    Group_add_segments(group, 0, group->plr.segments, group->plr.num_segments, true);  

    
    if (0 && DeBUG) {
        femu_log("[Group_update]: 此次segs索引段添加到group完成, 此时的group[%d]信息如下:\n", group->group_id);
        Print_Group(group);
    }
    plr_destroy(&group->plr);
}

void Print_Group(Group *group) {
    for (int i = 0; i < group->num_levels; i++) {
        for (int j = 0; j < group->L[i].num_segments; j++) {
            Segment seg = group->L[i].segments[j];
            femu_log("[Group_update] group[%d]->L[%d][Segment %d]: k: %.2f, b: %.2f, x1: %u, x2: %u, %s, %s\n", group->group_id, 
            i, j, seg.k, seg.b, seg.x1, seg.x2, seg.accurate ? "精确":"不精确", seg.consecutive ? "连续":"不连续");            
        }
    }    
}

void FrameGroup_init(FrameGroup *framegroup, double gamma) {
    framegroup->gamma = gamma;
    framegroup->num_segments = 0;
    framegroup->frame_length = 256;  // 1B maybe不够 frame_LBA_length
    framegroup->max_size = Mapping_TABLE_SIZE;
    
    framegroup->counter.group_write_cnt = 0;
    framegroup->counter.group_read_cnt = 0;
    framegroup->counter.group_read_acc_hit = 0;
    framegroup->counter.group_read_noacc_hit = 0;
    framegroup->counter.group_read_miss = 0;
    framegroup->counter.group_read_noacc_miss = 0;
    framegroup->counter.group_double_read = 0;

    framegroup->cnt_groups = 0;

    // framegroup->num_groups = framegroup->frame_length;
    // framegroup->num_groups = 16 * 1024 * 1024; // 2^32 / 2^8 = 2^24 3B
    // 16 GB / 4KB = 4M = 4 * 1024 * 1024 = 4194304
    // 64 GB = 16 * 1024 * 1024 => ngroups = 16 * 4 * 1024
    framegroup->num_groups = 32 * 16 * 1024;  // 512G

    framegroup->groups = (Group *)malloc(framegroup->num_groups * sizeof(Group));
    for (int i = 0; i < framegroup->num_groups; i++) {
        Group_init(&framegroup->groups[i], framegroup->gamma, i);
    }

    return ;
}

void FrameGroup_update(FrameGroup *framegroup, Point* points, int num_points) {
    Split_Points spt;
    spt.num_split = 0;
    spt.s_points = NULL;
   // spt.st = (bool *)malloc(framegroup->num_groups * sizeof(bool));
    spt.pos_split = (int *)malloc((framegroup->num_groups + 2) * sizeof(int));
    for (int i = 0; i < framegroup->num_groups + 2; i++) {
        spt.pos_split[i] = -1;
    }

    //femu_log("[FrameGroup_update]: num_points:%d\n", num_points);
    for (int i = 0; i < num_points; i++) {
        uint32_t split_id = points[i].x / framegroup->frame_length;
        uint8_t split_lpa = points[i].x % framegroup->frame_length;   // 好像不用
        uint32_t split_ppa = points[i].y;

        if (0 && DeBUG) {
            femu_log("[FrameGroup_update]: [LOOP %d] num_split:%d, point.x: %u, split_LPA: %u, split_id: %u\n", 
                        i, spt.num_split, points[i].x, split_lpa, split_id);
        }

        if (split_id >= framegroup->num_groups) {
            femu_log("[FrameGroup_update]: lba范围溢出\n");
            continue;
        }
        
        if (spt.pos_split[split_id] != -1) {
            int pos = spt.pos_split[split_id];
            spt.s_points[pos].points[spt.s_points[pos].num_points].x = split_lpa;
            spt.s_points[pos].points[spt.s_points[pos].num_points].y = split_ppa;
            spt.s_points[pos].num_points++;
        } else {
            if (spt.s_points == NULL) {
                spt.s_points = (Points *)malloc((spt.num_split + 1) * sizeof(Points));
                spt.s_points[spt.num_split].num_points = 0;
                spt.s_points[spt.num_split].split_id = split_id;
                spt.s_points[spt.num_split].points[0].x = split_lpa;
                spt.s_points[spt.num_split].points[0].y = split_ppa;
                spt.s_points[spt.num_split].num_points++;
                
                spt.pos_split[split_id] = spt.num_split;
                spt.num_split++;                
            } else {
                spt.s_points = (Points *)realloc(spt.s_points, (spt.num_split + 1) * sizeof(Points));
                
                spt.s_points[spt.num_split].num_points = 0;
                spt.s_points[spt.num_split].split_id = split_id;
                spt.s_points[spt.num_split].points[0].x = split_lpa;
                spt.s_points[spt.num_split].points[0].y = split_ppa;
                spt.s_points[spt.num_split].num_points++;
                
                spt.pos_split[split_id] = spt.num_split;
                spt.num_split++;
            }
        }
    }

    for (int i = 0; i < spt.num_split; i++) {
        Point *s_points = spt.s_points[i].points;
        int num_s_points = spt.s_points[i].num_points;
        int group_id = spt.s_points[i].split_id;
        //femu_log("-----------[FrameGroup_update]: spt.num_split:%d------------\n", i);
        Group_update(&framegroup->groups[group_id], s_points, num_s_points);
    }
    free(spt.s_points);
    free(spt.pos_split);
   // femu_log("[FrameGroup_update]: end !!!!!!!\n");
}

void FrameGroup_static(FrameGroup *framegroup) {
    uint32_t num = 0, cnt = 0;
    for (int i = 0; i < framegroup->num_groups; i++) {
        for (int j = 0; j < framegroup->groups[i].num_levels; j++) {
            num += framegroup->groups[i].L[j].num_segments;   
            cnt ++;
        }
    }
    uint32_t o_cnt = 0;
    for (int i = 0; i < framegroup->o_maptbl_maxLBA + 10; i++) {
        if (framegroup->o_maptbl[i]) o_cnt++;
    }
    framegroup->num_segments = num;
    framegroup->cnt_groups   = cnt;
    femu_log("[FrameGroup_static][lea_write]: 误差范围: %.1f, w_cnt: %d, framegroup->cnt_groups: %d, 总共的学习索引段数量: %d, o_maptbl: %lu, o_maptbl_maxLBA: %lu, o_maptbl_minLBA: %lu\n", 
                                    framegroup->gamma, framegroup->counter.group_write_cnt,
                                    framegroup->cnt_groups, framegroup->num_segments, 
                                    o_cnt, framegroup->o_maptbl_maxLBA, framegroup->o_maptbl_minLBA);
    femu_log("[FrameGroup_static][lea_read] : r_cnt: %d, lac_hit: %d, nac_hit: %d, nac_miss: %d, double_read: %d, omap_hit:%lu, miss: %d\n", 
                                    framegroup->counter.group_read_cnt, 
                                    framegroup->counter.group_read_acc_hit, framegroup->counter.group_read_noacc_hit,
                                    framegroup->counter.group_read_noacc_miss, framegroup->counter.group_double_read,
                                    framegroup->o_maptbl_hit, framegroup->counter.group_read_miss);    
}

uint64_t FrameGroup_lookup(FrameGroup *framegroup, uint64_t lpn, Segment *seg)
{
    uint32_t group_id = lpn / framegroup->frame_length;
    uint8_t find_lpa = lpn % framegroup->frame_length;   
    uint64_t ppa = 0;

    if (group_id >= framegroup->num_groups) {
        femu_log("[FrameGroup_lookup]: lba范围溢出\n");
        return ppa;
    }
    Group group = framegroup->groups[group_id];
    for (int i = 0; i < group.num_levels; i++) {
        LogPLR logplr = group.L[i];
        Segment find_seg;
        Segment_init(&find_seg, 0, 0, find_lpa, find_lpa, NULL, 0);
        int pos = 0;
        // 第一个大于等于tar的位置
        pos = binary_search(&logplr, &find_seg);
        if (pos == 0 || (pos < logplr.num_segments && logplr.segments[pos].x1 == find_lpa))
            pos = pos;
        else 
            pos -= 1;
        find_seg = logplr.segments[pos];

        //ppa = Segment_gety(&find_seg, true, find_lpa);
        // 上下取整预测
        ppa = Segment_gety_upper(&find_seg, true, find_lpa);

        if (ppa == 0) {
           // femu_log("group[%d]->L[%d] 未找到 lpa: %u\n", group.group_id, i, find_lpa);
            continue;
        } else {
            // TODO
            // 如果是不精确的还得先检查CRB 
            // Check bit_map
            *seg = find_seg;
            break;
        }
    }
    if (ppa == 0) {
       // femu_log("group[%d]上未找到对于 lpn: %lu的写入信息\n", group.group_id, lpn);
    }

    return ppa;
}


/**
 * @brief ssd_method
 * 
*/

static void *ftl_thread(void *arg);

static inline bool should_gc(struct ssd *ssd)
{
    return (ssd->lm.free_line_cnt <= ssd->sp.gc_thres_lines);
}

static inline bool should_gc_high(struct ssd *ssd)
{
    return (ssd->lm.free_line_cnt <= ssd->sp.gc_thres_lines_high);
}

static inline struct ppa get_maptbl_ent(struct ssd *ssd, uint64_t lpn)
{
    return ssd->maptbl[lpn];
}

static inline void set_maptbl_ent(struct ssd *ssd, uint64_t lpn, struct ppa *ppa)
{
    assert(lpn < ssd->sp.tt_pgs);
    ssd->maptbl[lpn] = *ppa;
}

static uint64_t ppa2pgidx(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    uint64_t pgidx;

    pgidx = ppa->g.ch  * spp->pgs_per_ch  + \
            ppa->g.lun * spp->pgs_per_lun + \
            ppa->g.pl  * spp->pgs_per_pl  + \
            ppa->g.blk * spp->pgs_per_blk + \
            ppa->g.pg;

    assert(pgidx < spp->tt_pgs);

    return pgidx;
}

static inline struct ppa pgidx2ppa(struct ssd *ssd, uint64_t pgidx)
{
    struct ssdparams *spp = &ssd->sp;
    struct ppa ppa;
    ppa.ppa = 0;
    ppa.g.ch = pgidx / spp->pgs_per_ch;
    ppa.g.lun = (pgidx - ppa.g.ch * spp->pgs_per_ch) / spp->pgs_per_lun;
    ppa.g.pl = (pgidx - ppa.g.ch * spp->pgs_per_ch - ppa.g.lun * spp->pgs_per_lun) / spp->pgs_per_pl;
    ppa.g.blk = (pgidx - ppa.g.ch * spp->pgs_per_ch - ppa.g.lun * spp->pgs_per_lun - 
                ppa.g.pl  * spp->pgs_per_pl) / spp->pgs_per_blk;
    ppa.g.pg = (pgidx - ppa.g.ch * spp->pgs_per_ch - ppa.g.lun * spp->pgs_per_lun - 
                ppa.g.pl  * spp->pgs_per_pl - ppa.g.blk * spp->pgs_per_blk);

    assert(pgidx < spp->tt_pgs);

    return ppa;
}

uint64_t ppa2vpgidx(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    uint64_t vpgidx;

    // pgidx = ppa->g.ch  * spp->pgs_per_ch  + 
    //         ppa->g.lun * spp->pgs_per_lun + 
    //         ppa->g.pl  * spp->pgs_per_pl  + 
    //         ppa->g.blk * spp->pgs_per_blk + 
    //         ppa->g.pg;

    
    vpgidx = ppa->g.blk * spp->pgs_per_line  + \
             ppa->g.pg  * spp->blks_per_line + \
             ppa->g.pl  * spp->luns_per_line + \
             ppa->g.lun * spp->nchs          + \
             ppa->g.ch;

    assert(vpgidx < spp->tt_pgs);

    return vpgidx;
}

static inline struct ppa vpgidx2ppa(struct ssd *ssd, uint64_t vpgidx)
{
    struct ssdparams *spp = &ssd->sp;
    struct ppa ppa;
    ppa.ppa = 0;
    // ppa.g.ch = pgidx / spp->pgs_per_ch;
    // ppa.g.lun = (pgidx - ppa.g.ch * spp->pgs_per_ch) / spp->pgs_per_lun;
    // ppa.g.pl = (pgidx - ppa.g.ch * spp->pgs_per_ch - ppa.g.lun * spp->pgs_per_lun) / spp->pgs_per_pl;
    // ppa.g.blk = (pgidx - ppa.g.ch * spp->pgs_per_ch - ppa.g.lun * spp->pgs_per_lun - 
    //             ppa.g.pl  * spp->pgs_per_pl) / spp->pgs_per_blk;
    // ppa.g.pg = (pgidx - ppa.g.ch * spp->pgs_per_ch - ppa.g.lun * spp->pgs_per_lun - 
    //             ppa.g.pl  * spp->pgs_per_pl - ppa.g.blk * spp->pgs_per_blk);

    ppa.g.blk = vpgidx / spp->pgs_per_line;
    ppa.g.pg  = (vpgidx - ppa.g.blk * spp->pgs_per_line) / spp->blks_per_line;
    ppa.g.pl  = (vpgidx - ppa.g.blk * spp->pgs_per_line - ppa.g.pg * spp->blks_per_line) / spp->luns_per_line;
    ppa.g.lun = (vpgidx - ppa.g.blk * spp->pgs_per_line - ppa.g.pg * spp->blks_per_line -
                 ppa.g.pl * spp->luns_per_line) / spp->nchs;
    ppa.g.ch  = (vpgidx - ppa.g.blk * spp->pgs_per_line - ppa.g.pg * spp->blks_per_line -
                 ppa.g.pl * spp->luns_per_line - ppa.g.lun * spp->nchs);

    assert(vpgidx < spp->tt_pgs);

    return ppa;
}



static inline uint64_t get_rmap_ent(struct ssd *ssd, struct ppa *ppa)
{
    uint64_t pgidx = ppa2pgidx(ssd, ppa);

    return ssd->rmap[pgidx];
}

/* set rmap[page_no(ppa)] -> lpn */
static inline void set_rmap_ent(struct ssd *ssd, uint64_t lpn, struct ppa *ppa)
{
    uint64_t pgidx = ppa2pgidx(ssd, ppa);

    ssd->rmap[pgidx] = lpn;
}

static inline int victim_line_cmp_pri(pqueue_pri_t next, pqueue_pri_t curr)
{
    return (next > curr);
}

static inline pqueue_pri_t victim_line_get_pri(void *a)
{
    return ((struct line *)a)->vpc;
}

static inline void victim_line_set_pri(void *a, pqueue_pri_t pri)
{
    ((struct line *)a)->vpc = pri;
}

static inline size_t victim_line_get_pos(void *a)
{
    return ((struct line *)a)->pos;
}

static inline void victim_line_set_pos(void *a, size_t pos)
{
    ((struct line *)a)->pos = pos;
}

static void ssd_init_lines(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;
    struct line_mgmt *lm = &ssd->lm;
    struct line *line;

    lm->tt_lines = spp->blks_per_pl;
    assert(lm->tt_lines == spp->tt_lines);
    lm->lines = g_malloc0(sizeof(struct line) * lm->tt_lines);

    QTAILQ_INIT(&lm->free_line_list);
    lm->victim_line_pq = pqueue_init(spp->tt_lines, victim_line_cmp_pri,
            victim_line_get_pri, victim_line_set_pri,
            victim_line_get_pos, victim_line_set_pos);
    QTAILQ_INIT(&lm->full_line_list);

    lm->free_line_cnt = 0;
    for (int i = 0; i < lm->tt_lines; i++) {
        line = &lm->lines[i];
        line->id = i;
        line->ipc = 0;
        line->vpc = 0;
        line->pos = 0;
        /* initialize all the lines as free lines */
        QTAILQ_INSERT_TAIL(&lm->free_line_list, line, entry);
        lm->free_line_cnt++;
    }

    assert(lm->free_line_cnt == lm->tt_lines);
    lm->victim_line_cnt = 0;
    lm->full_line_cnt = 0;
}

static void ssd_init_write_pointer(struct ssd *ssd)
{
    struct write_pointer *wpp = &ssd->wp;
    struct line_mgmt *lm = &ssd->lm;
    struct line *curline = NULL;

    curline = QTAILQ_FIRST(&lm->free_line_list);
    QTAILQ_REMOVE(&lm->free_line_list, curline, entry);
    lm->free_line_cnt--;

    /* wpp->curline is always our next-to-write super-block */
    wpp->curline = curline;
    wpp->ch = 0;
    wpp->lun = 0;
    wpp->pg = 0;
    wpp->blk = 0;
    wpp->pl = 0;
}

static inline void check_addr(int a, int max)
{
    assert(a >= 0 && a < max);
}

static struct line *get_next_free_line(struct ssd *ssd)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *curline = NULL;

    curline = QTAILQ_FIRST(&lm->free_line_list);
    if (!curline) {
        FrameGroup_static(&ssd->l_maptbl);
        femu_err("No free lines left in [%s] !!!!\n", ssd->ssdname);
        return NULL;
    }

    QTAILQ_REMOVE(&lm->free_line_list, curline, entry);
    lm->free_line_cnt--;
    return curline;
}

void ssd_advance_write_pointer(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;
    struct write_pointer *wpp = &ssd->wp;
    struct line_mgmt *lm = &ssd->lm;

    check_addr(wpp->ch, spp->nchs);
    wpp->ch++;
    if (wpp->ch == spp->nchs) {
        wpp->ch = 0;
        check_addr(wpp->lun, spp->luns_per_ch);
        wpp->lun++;
        /* in this case, we should go to next lun */
        if (wpp->lun == spp->luns_per_ch) {
            wpp->lun = 0;
            /* go to next page in the block */
            check_addr(wpp->pg, spp->pgs_per_blk);
            wpp->pg++;
            if (wpp->pg == spp->pgs_per_blk) {
                wpp->pg = 0;
                /* move current line to {victim,full} line list */
                // 超级块写满了
                if (wpp->curline->vpc == spp->pgs_per_line) {
                    /* all pgs are still valid, move to full line list */
                    assert(wpp->curline->ipc == 0);
                    QTAILQ_INSERT_TAIL(&lm->full_line_list, wpp->curline, entry);
                    lm->full_line_cnt++;
                } else {
                    assert(wpp->curline->vpc >= 0 && wpp->curline->vpc < spp->pgs_per_line);
                    /* there must be some invalid pages in this line */
                    assert(wpp->curline->ipc > 0);
                    // 将超级块放到gc victim优先队列中
                    pqueue_insert(lm->victim_line_pq, wpp->curline);
                    lm->victim_line_cnt++;
                }
                /* current line is used up, pick another empty line */
                check_addr(wpp->blk, spp->blks_per_pl);
                wpp->curline = NULL;
                wpp->curline = get_next_free_line(ssd);
                if (!wpp->curline) {
                    /* TODO */
                    abort();
                }
                wpp->blk = wpp->curline->id;
                check_addr(wpp->blk, spp->blks_per_pl);
                /* make sure we are starting from page 0 in the super block */
                assert(wpp->pg == 0);
                assert(wpp->lun == 0);
                assert(wpp->ch == 0);
                /* TODO: assume # of pl_per_lun is 1, fix later */
                assert(wpp->pl == 0);
            }
        }
    }
}

struct ppa get_new_page(struct ssd *ssd)
{
    struct write_pointer *wpp = &ssd->wp;
    struct ppa ppa;
    ppa.ppa = 0;
    ppa.g.ch = wpp->ch;
    ppa.g.lun = wpp->lun;
    ppa.g.pg = wpp->pg;
    ppa.g.blk = wpp->blk;
    ppa.g.pl = wpp->pl;
    assert(ppa.g.pl == 0);

    return ppa;
}

static void check_params(struct ssdparams *spp)
{
    /*
     * we are using a general write pointer increment method now, no need to
     * force luns_per_ch and nchs to be power of 2
     */

    //assert(is_power_of_2(spp->luns_per_ch));
    //assert(is_power_of_2(spp->nchs));
}

static void ssd_init_params(struct ssdparams *spp)
{
    spp->secsz = 512;
    spp->secs_per_pg = 8;
    spp->pgs_per_blk = 2048; // *2 *  = 512 GB
    spp->blks_per_pl = 1024; /*2 = 64GB */
    spp->pls_per_lun = 1;
    spp->luns_per_ch = 8;
    spp->nchs = 8;

    spp->pg_rd_lat = NAND_READ_LATENCY;
    spp->pg_wr_lat = NAND_PROG_LATENCY;
    spp->blk_er_lat = NAND_ERASE_LATENCY;
    spp->ch_xfer_lat = 0;

    /* calculated values */
    spp->secs_per_blk = spp->secs_per_pg * spp->pgs_per_blk;
    spp->secs_per_pl = spp->secs_per_blk * spp->blks_per_pl;
    spp->secs_per_lun = spp->secs_per_pl * spp->pls_per_lun;
    spp->secs_per_ch = spp->secs_per_lun * spp->luns_per_ch;
    spp->tt_secs = spp->secs_per_ch * spp->nchs;

    spp->pgs_per_pl = spp->pgs_per_blk * spp->blks_per_pl;
    spp->pgs_per_lun = spp->pgs_per_pl * spp->pls_per_lun;
    spp->pgs_per_ch = spp->pgs_per_lun * spp->luns_per_ch;
    spp->tt_pgs = spp->pgs_per_ch * spp->nchs;

    spp->blks_per_lun = spp->blks_per_pl * spp->pls_per_lun;
    spp->blks_per_ch = spp->blks_per_lun * spp->luns_per_ch;
    spp->tt_blks = spp->blks_per_ch * spp->nchs;

    spp->pls_per_ch =  spp->pls_per_lun * spp->luns_per_ch;
    spp->tt_pls = spp->pls_per_ch * spp->nchs;

    spp->tt_luns = spp->luns_per_ch * spp->nchs;  // 16GB  8 * 8 = 64

    /* line is special, put it at the end */
    spp->blks_per_line = spp->tt_luns * 1; /* TODO: to fix under multiplanes */ // no多通道 plane = 1
    spp->pgs_per_line = spp->blks_per_line * spp->pgs_per_blk;
    spp->secs_per_line = spp->pgs_per_line * spp->secs_per_pg;
    spp->tt_lines = spp->blks_per_lun; /* TODO: to fix under multiplanes */

    // super_block == line
    // 32 GB
    spp->luns_per_line = spp->nchs * spp->luns_per_ch;          // 64

    spp->pls_per_line = spp->luns_per_line * spp->pls_per_lun;  // 64 
    spp->blks_per_line = spp->pls_per_line;                     // 64 (blks_per_line = pages_per_superpagd) 

    spp->pgs_per_line = spp->blks_per_line * spp->pgs_per_blk;  // 64 * 512 = 32768 


    spp->gc_thres_pcent = 0.75;
    spp->gc_thres_lines = (int)((1 - spp->gc_thres_pcent) * spp->tt_lines);
    spp->gc_thres_pcent_high = 0.95;
    spp->gc_thres_lines_high = (int)((1 - spp->gc_thres_pcent_high) * spp->tt_lines);
    spp->enable_gc_delay = true;


    check_params(spp);
}

static void ssd_init_nand_page(struct nand_page *pg, struct ssdparams *spp)
{
    pg->nsecs = spp->secs_per_pg;
    pg->sec = g_malloc0(sizeof(nand_sec_status_t) * pg->nsecs);
    for (int i = 0; i < pg->nsecs; i++) {
        pg->sec[i] = SEC_FREE;
    }
    pg->status = PG_FREE;
}

static void ssd_init_nand_blk(struct nand_block *blk, struct ssdparams *spp)
{
    blk->npgs = spp->pgs_per_blk;
    blk->pg = g_malloc0(sizeof(struct nand_page) * blk->npgs);
    for (int i = 0; i < blk->npgs; i++) {
        ssd_init_nand_page(&blk->pg[i], spp);
    }
    blk->ipc = 0;
    blk->vpc = 0;
    blk->erase_cnt = 0;
    blk->wp = 0;
}

static void ssd_init_nand_plane(struct nand_plane *pl, struct ssdparams *spp)
{
    pl->nblks = spp->blks_per_pl;
    pl->blk = g_malloc0(sizeof(struct nand_block) * pl->nblks);
    for (int i = 0; i < pl->nblks; i++) {
        ssd_init_nand_blk(&pl->blk[i], spp);
    }
}

static void ssd_init_nand_lun(struct nand_lun *lun, struct ssdparams *spp)
{
    lun->npls = spp->pls_per_lun;
    lun->pl = g_malloc0(sizeof(struct nand_plane) * lun->npls);
    for (int i = 0; i < lun->npls; i++) {
        ssd_init_nand_plane(&lun->pl[i], spp);
    }
    lun->next_lun_avail_time = 0;
    lun->busy = false;
}

static void ssd_init_ch(struct ssd_channel *ch, struct ssdparams *spp)
{
    ch->nluns = spp->luns_per_ch;
    ch->lun = g_malloc0(sizeof(struct nand_lun) * ch->nluns);
    for (int i = 0; i < ch->nluns; i++) {
        ssd_init_nand_lun(&ch->lun[i], spp);
    }
    ch->next_ch_avail_time = 0;
    ch->busy = 0;
}

static void ssd_init_maptbl(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;

    ssd->maptbl = g_malloc0(sizeof(struct ppa) * spp->tt_pgs);
    for (int i = 0; i < spp->tt_pgs; i++) {
        ssd->maptbl[i].ppa = UNMAPPED_PPA;
    }

}

static void ssd_init_rmap(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;

    ssd->rmap = g_malloc0(sizeof(uint64_t) * spp->tt_pgs);
    for (int i = 0; i < spp->tt_pgs; i++) {
        ssd->rmap[i] = INVALID_LPN;
    }
}

static void ssd_init_o_maptbl(struct ssd *ssd) {
    struct ssdparams *spp = &ssd->sp;
    ssd->l_maptbl.o_maptbl = g_malloc0(sizeof(uint64_t) * spp->tt_pgs);
    for (int i = 0; i < spp->tt_pgs; i++) {
        ssd->l_maptbl.o_maptbl[i] = 0;
    }
    ssd->l_maptbl.o_maptbl_cnt = 0;
    ssd->l_maptbl.o_maptbl_hit = 0;

    ssd->l_maptbl.o_maptbl_maxLBA = 0;
    ssd->l_maptbl.o_maptbl_minLBA = 16777216; // 64GB
    
}

void ssd_init(FemuCtrl *n)
{
    struct ssd *ssd = n->ssd;
    struct ssdparams *spp = &ssd->sp;

    assert(ssd);


    ssd_init_params(spp);
    
    

    /* initialize ssd internal layout architecture */
    ssd->ch = g_malloc0(sizeof(struct ssd_channel) * spp->nchs);
    for (int i = 0; i < spp->nchs; i++) {
        ssd_init_ch(&ssd->ch[i], spp);
    }

    /* initialize maptbl */
    ssd_init_maptbl(ssd);


    // 初始化日志结构索引映射表    
    FrameGroup_init(&ssd->l_maptbl, error_bound);
    ssd_init_o_maptbl(ssd);
    ssd->num_write_entries = 0;
    ssd->flush = false;
    for (int i = 0; i < WB_Entries + 2; i++) {
        ssd->WB[i].LPA = 0;
        ssd->WB[i].PPA = 0;
    }
    ssd->hit_wb = 0;



    /* initialize rmap */
    ssd_init_rmap(ssd);

    /* initialize all the lines */
    ssd_init_lines(ssd);

    /* initialize write pointer, this is how we allocate new pages for writes */
    ssd_init_write_pointer(ssd);

    qemu_thread_create(&ssd->ftl_thread, "FEMU-FTL-Thread", ftl_thread, n,
                       QEMU_THREAD_JOINABLE);
}

static inline bool valid_ppa(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    int ch = ppa->g.ch;
    int lun = ppa->g.lun;
    int pl = ppa->g.pl;
    int blk = ppa->g.blk;
    int pg = ppa->g.pg;
    int sec = ppa->g.sec;

    if (ch >= 0 && ch < spp->nchs && lun >= 0 && lun < spp->luns_per_ch && pl >=
        0 && pl < spp->pls_per_lun && blk >= 0 && blk < spp->blks_per_pl && pg
        >= 0 && pg < spp->pgs_per_blk && sec >= 0 && sec < spp->secs_per_pg)
        return true;

    return false;
}

static inline bool valid_lpn(struct ssd *ssd, uint64_t lpn)
{
    return (lpn < ssd->sp.tt_pgs);
}

static inline bool mapped_ppa(struct ppa *ppa)
{
    return !(ppa->ppa == UNMAPPED_PPA);
}

static inline struct ssd_channel *get_ch(struct ssd *ssd, struct ppa *ppa)
{
    return &(ssd->ch[ppa->g.ch]);
}

static inline struct nand_lun *get_lun(struct ssd *ssd, struct ppa *ppa)
{
    struct ssd_channel *ch = get_ch(ssd, ppa);
    return &(ch->lun[ppa->g.lun]);
}

static inline struct nand_plane *get_pl(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_lun *lun = get_lun(ssd, ppa);
    return &(lun->pl[ppa->g.pl]);
}

static inline struct nand_block *get_blk(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_plane *pl = get_pl(ssd, ppa);
    return &(pl->blk[ppa->g.blk]);
}

static inline struct line *get_line(struct ssd *ssd, struct ppa *ppa)
{
    return &(ssd->lm.lines[ppa->g.blk]);
}

static inline struct nand_page *get_pg(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_block *blk = get_blk(ssd, ppa);
    return &(blk->pg[ppa->g.pg]);
}

static uint64_t ssd_advance_status(struct ssd *ssd, struct ppa *ppa, struct
        nand_cmd *ncmd)
{
    int c = ncmd->cmd;
    uint64_t cmd_stime = (ncmd->stime == 0) ? \
        qemu_clock_get_ns(QEMU_CLOCK_REALTIME) : ncmd->stime;
    uint64_t nand_stime;
    struct ssdparams *spp = &ssd->sp;
    struct nand_lun *lun = get_lun(ssd, ppa);
    uint64_t lat = 0;

    switch (c) {
    case NAND_READ:
        /* read: perform NAND cmd first */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->pg_rd_lat;
        lat = lun->next_lun_avail_time - cmd_stime;
#if 0
        lun->next_lun_avail_time = nand_stime + spp->pg_rd_lat;

        /* read: then data transfer through channel */
        chnl_stime = (ch->next_ch_avail_time < lun->next_lun_avail_time) ? \
            lun->next_lun_avail_time : ch->next_ch_avail_time;
        ch->next_ch_avail_time = chnl_stime + spp->ch_xfer_lat;

        lat = ch->next_ch_avail_time - cmd_stime;
#endif
        break;

    case NAND_WRITE:
        /* write: transfer data through channel first */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        if (ncmd->type == USER_IO) {
            lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;
        } else {
            lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;
        }
        lat = lun->next_lun_avail_time - cmd_stime;

#if 0
        chnl_stime = (ch->next_ch_avail_time < cmd_stime) ? cmd_stime : \
                     ch->next_ch_avail_time;
        ch->next_ch_avail_time = chnl_stime + spp->ch_xfer_lat;

        /* write: then do NAND program */
        nand_stime = (lun->next_lun_avail_time < ch->next_ch_avail_time) ? \
            ch->next_ch_avail_time : lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;

        lat = lun->next_lun_avail_time - cmd_stime;
#endif
        break;

    case NAND_ERASE:
        /* erase: only need to advance NAND status */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->blk_er_lat;

        lat = lun->next_lun_avail_time - cmd_stime;
        break;

    default:
        femu_err("Unsupported NAND command: 0x%x\n", c);
    }

    return lat;
}

/* update SSD status about one page from PG_VALID -> PG_INVALID */
static void mark_page_invalid(struct ssd *ssd, struct ppa *ppa)
{
    struct line_mgmt *lm = &ssd->lm;
    struct ssdparams *spp = &ssd->sp;
    struct nand_block *blk = NULL;
    struct nand_page *pg = NULL;
    bool was_full_line = false;
    struct line *line;

    /* update corresponding page status */
    pg = get_pg(ssd, ppa);
    assert(pg->status == PG_VALID);
    pg->status = PG_INVALID;

    /* update corresponding block status */
    blk = get_blk(ssd, ppa);
    assert(blk->ipc >= 0 && blk->ipc < spp->pgs_per_blk);
    blk->ipc++;
    assert(blk->vpc > 0 && blk->vpc <= spp->pgs_per_blk);
    blk->vpc--;

    /* update corresponding line status */
    line = get_line(ssd, ppa);
    assert(line->ipc >= 0 && line->ipc < spp->pgs_per_line);
    if (line->vpc == spp->pgs_per_line) {
        assert(line->ipc == 0);
        was_full_line = true;
    }
    line->ipc++;
    assert(line->vpc > 0 && line->vpc <= spp->pgs_per_line);
    /* Adjust the position of the victime line in the pq under over-writes */
    if (line->pos) {
        /* Note that line->vpc will be updated by this call */
        pqueue_change_priority(lm->victim_line_pq, line->vpc - 1, line);
    } else {
        line->vpc--;
    }

    if (was_full_line) {
        /* move line: "full" -> "victim" */
        QTAILQ_REMOVE(&lm->full_line_list, line, entry);
        lm->full_line_cnt--;
        pqueue_insert(lm->victim_line_pq, line);
        lm->victim_line_cnt++;
    }
}

static void mark_page_valid(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_block *blk = NULL;
    struct nand_page *pg = NULL;
    struct line *line;

    /* update page status */
    pg = get_pg(ssd, ppa);
    assert(pg->status == PG_FREE);
    pg->status = PG_VALID;

    /* update corresponding block status */
    blk = get_blk(ssd, ppa);
    assert(blk->vpc >= 0 && blk->vpc < ssd->sp.pgs_per_blk);
    blk->vpc++;

    /* update corresponding line status */
    line = get_line(ssd, ppa);
    assert(line->vpc >= 0 && line->vpc < ssd->sp.pgs_per_line);
    line->vpc++;
}

static void mark_block_free(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    struct nand_block *blk = get_blk(ssd, ppa);
    struct nand_page *pg = NULL;

    for (int i = 0; i < spp->pgs_per_blk; i++) {
        /* reset page status */
        pg = &blk->pg[i];
        assert(pg->nsecs == spp->secs_per_pg);
        pg->status = PG_FREE;
    }

    /* reset block status */
    assert(blk->npgs == spp->pgs_per_blk);
    blk->ipc = 0;
    blk->vpc = 0;
    blk->erase_cnt++;
}

static void gc_read_page(struct ssd *ssd, struct ppa *ppa)
{
    /* advance ssd status, we don't care about how long it takes */
    if (ssd->sp.enable_gc_delay) {
        struct nand_cmd gcr;
        gcr.type = GC_IO;
        gcr.cmd = NAND_READ;
        gcr.stime = 0;
        ssd_advance_status(ssd, ppa, &gcr);
    }
}

/* move valid page data (already in DRAM) from victim line to a new page */
static uint64_t gc_write_page(struct ssd *ssd, struct ppa *old_ppa)
{
    struct ppa new_ppa;
    struct nand_lun *new_lun;
    uint64_t lpn = get_rmap_ent(ssd, old_ppa);

    assert(valid_lpn(ssd, lpn));
    new_ppa = get_new_page(ssd);
    /* update maptbl */
    set_maptbl_ent(ssd, lpn, &new_ppa);
    /* update rmap */
    set_rmap_ent(ssd, lpn, &new_ppa);

    mark_page_valid(ssd, &new_ppa);

    /* need to advance the write pointer here */
    ssd_advance_write_pointer(ssd);

    if (ssd->sp.enable_gc_delay) {
        struct nand_cmd gcw;
        gcw.type = GC_IO;
        gcw.cmd = NAND_WRITE;
        gcw.stime = 0;
        ssd_advance_status(ssd, &new_ppa, &gcw);
    }

    /* advance per-ch gc_endtime as well */
#if 0
    new_ch = get_ch(ssd, &new_ppa);
    new_ch->gc_endtime = new_ch->next_ch_avail_time;
#endif

    new_lun = get_lun(ssd, &new_ppa);
    new_lun->gc_endtime = new_lun->next_lun_avail_time;

    return 0;
}

static struct line *select_victim_line(struct ssd *ssd, bool force)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *victim_line = NULL;

    victim_line = pqueue_peek(lm->victim_line_pq);
    if (!victim_line) {
        return NULL;
    }

    if (!force && victim_line->ipc < ssd->sp.pgs_per_line / 8) {
        return NULL;
    }

    pqueue_pop(lm->victim_line_pq);
    victim_line->pos = 0;
    lm->victim_line_cnt--;

    /* victim_line is a danggling node now */
    return victim_line;
}

/* here ppa identifies the block we want to clean */
static void clean_one_block(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    struct nand_page *pg_iter = NULL;
    int cnt = 0;

    for (int pg = 0; pg < spp->pgs_per_blk; pg++) {
        ppa->g.pg = pg;
        pg_iter = get_pg(ssd, ppa);
        /* there shouldn't be any free page in victim blocks */
        assert(pg_iter->status != PG_FREE);
        if (pg_iter->status == PG_VALID) {
            gc_read_page(ssd, ppa);
            /* delay the maptbl update until "write" happens */
            gc_write_page(ssd, ppa);
            cnt++;
        }
    }

    assert(get_blk(ssd, ppa)->vpc == cnt);
}

static void mark_line_free(struct ssd *ssd, struct ppa *ppa)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *line = get_line(ssd, ppa);
    line->ipc = 0;
    line->vpc = 0;
    /* move this line to free line list */
    QTAILQ_INSERT_TAIL(&lm->free_line_list, line, entry);
    lm->free_line_cnt++;
}

static int do_gc(struct ssd *ssd, bool force)
{
    struct line *victim_line = NULL;
    struct ssdparams *spp = &ssd->sp;
    struct nand_lun *lunp;
    struct ppa ppa;
    int ch, lun;

    victim_line = select_victim_line(ssd, force);
    if (!victim_line) {
        return -1;
    }

    ppa.g.blk = victim_line->id;
    femu_debug("GC-ing line:%d,ipc=%d,victim=%d,full=%d,free=%d\n", ppa.g.blk,
              victim_line->ipc, ssd->lm.victim_line_cnt, ssd->lm.full_line_cnt,
              ssd->lm.free_line_cnt);

    /* copy back valid data */
    // superblok level clean block
    for (ch = 0; ch < spp->nchs; ch++) {
        for (lun = 0; lun < spp->luns_per_ch; lun++) {
            ppa.g.ch = ch;
            ppa.g.lun = lun;
            ppa.g.pl = 0;
            lunp = get_lun(ssd, &ppa);
            clean_one_block(ssd, &ppa);
            mark_block_free(ssd, &ppa);

            // 擦除时间
            if (spp->enable_gc_delay) {
                struct nand_cmd gce;
                gce.type = GC_IO;
                gce.cmd = NAND_ERASE;
                gce.stime = 0;
                ssd_advance_status(ssd, &ppa, &gce);
            }

            lunp->gc_endtime = lunp->next_lun_avail_time;
        }
    }

    /* update line status */
    mark_line_free(ssd, &ppa);

    return 0;
}

static uint64_t ssd_read(struct ssd *ssd, NvmeRequest *req)
{
    struct ssdparams *spp = &ssd->sp;
    uint64_t lba = req->slba;
    int nsecs = req->nlb;
    struct ppa ppa;
    uint64_t start_lpn = lba / spp->secs_per_pg;
    uint64_t end_lpn = (lba + nsecs - 1) / spp->secs_per_pg;
    uint64_t lpn;
    uint64_t sublat, maxlat = 0;
    uint64_t ret_ppa;
    uint64_t mem_lat = 0;

    if (end_lpn >= spp->tt_pgs) {
        return maxlat;
        femu_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
    }

    if (0 && DeBUG) {
        for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
            femu_log("[ssd_lea_read] READ_LPN: %lu\n", lpn);
        }
    }

    // 统计 验证
    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        if (ssd->l_maptbl.o_maptbl[lpn] == 1) {
            ssd->l_maptbl.o_maptbl_hit++;
        }
    }


    /* leaFTL IO read path */
    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        Segment seg;
        
        
        for (int i = 0; i < ssd->num_write_entries; i++) {
            if (lpn == ssd->WB[ssd->num_write_entries].LPA) {
                ssd->hit_wb++;
                maxlat += mem_lat;
                goto hit_wb;
            }
        }

       // uint64_t stime = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
        ret_ppa = FrameGroup_lookup(&ssd->l_maptbl, lpn, &seg);
       // uint64_t etime = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
        ssd->l_maptbl.counter.group_read_cnt++;

        if (!ret_ppa) {
            //printf("%s,lpn(%" PRId64 ") not mapped to valid ppa\n", ssd->ssdname, lpn);
            //printf("Invalid ppa,ch:%d,lun:%d,blk:%d,pl:%d,pg:%d,sec:%d\n",
            //ppa.g.ch, ppa.g.lun, ppa.g.blk, ppa.g.pl, ppa.g.pg, ppa.g.sec);
            ssd->l_maptbl.counter.group_read_miss++;
            continue;
        }
        else if (!valid_lpn(ssd, ret_ppa)) {
            femu_log("Invalid ppA: %lu\n", ret_ppa);
            continue;
        }
        else {
            //  femu_log("[lookup_cost]LA:%lu -> PA:%lu, et-st:%lu ns\n", lpn, ret_ppa, etime - stime);
            if (!seg.accurate) {
                //femu_log("seg不精确,寻找ppa: %lu 对应的oob\n", ret_ppa);
                uint64_t ed = (uint64_t)ssd->l_maptbl.gamma;
                uint64_t y;
                // uint64_t t_ppa;
                // struct ppa tmp_ppa = vpgidx2ppa(ssd, ret_ppa);
                // t_ppa = ppa2pgidx(ssd, &tmp_ppa);

                for (y = 0; y <= ed; y++) {
                    uint64_t i = ret_ppa - y, j = ret_ppa + y;  
                    if (ssd->rmap[i] == lpn || ssd->rmap[j] == lpn) {
                        if (y == 0) {
                            //femu_log("预测精确,直接读 ppa: %lu -> oob.lpa: %lu,\n", ret_ppa, lpn);
                            ppa = vpgidx2ppa(ssd, ret_ppa);
                            struct nand_cmd srd;
                            srd.type = USER_IO;
                            srd.cmd = NAND_READ;
                            srd.stime = req->stime;
                            sublat = ssd_advance_status(ssd, &ppa, &srd);
                            maxlat = (sublat > maxlat) ? sublat : maxlat;
                        } else {
                            if (0 && DeBUG) {
                                femu_log("预测错误，需要双读,第一次读 ppa: %lu -> oob.lpa: %lu,\n\t\t第二次读ppa: %lu -> oob.lpa: %lu\n ", ret_ppa, ssd->rmap[ret_ppa],
                                        ssd->rmap[i] == lpn ? i : j, lpn);
                            }
                            ppa = vpgidx2ppa(ssd, ret_ppa);
                            struct nand_cmd srd;
                            srd.type = USER_IO;
                            srd.cmd = NAND_READ;
                            srd.stime = req->stime;
                            sublat = ssd_advance_status(ssd, &ppa, &srd);
                            maxlat = (sublat > maxlat) ? sublat : maxlat;
                            
                            ret_ppa = ssd->rmap[i] == lpn ? i : j;
                            ppa = vpgidx2ppa(ssd, ret_ppa);
                            sublat = ssd_advance_status(ssd, &ppa, &srd);
                            maxlat = (sublat > maxlat) ? sublat : maxlat;
                            ssd->l_maptbl.counter.group_double_read++;

                        }
                        ssd->l_maptbl.counter.group_read_noacc_hit++;
                        break;
                    } 
                }
                if (y > ed) {
                   // femu_log("查找失败,寻找 ppa: %lu 的oob空间未找到对应的lpn: %lu\n", ret_ppa, lpn);
                    ssd->l_maptbl.counter.group_read_noacc_miss++;
                    // TODO 未作处理
                }
            } else {
                ssd->l_maptbl.counter.group_read_acc_hit++;
               // femu_log("seg精确,直接读取 lpn: %lu -> flash_ppa: %lu \n", lpn, ret_ppa);
                ppa = vpgidx2ppa(ssd, ret_ppa);
                struct nand_cmd srd;
                srd.type = USER_IO;
                srd.cmd = NAND_READ;
                srd.stime = req->stime;
                sublat = ssd_advance_status(ssd, &ppa, &srd);
                maxlat = (sublat > maxlat) ? sublat : maxlat;
            }
        }

    }

hit_wb:
    return maxlat;    
}

static uint64_t ssd_write(struct ssd *ssd, NvmeRequest *req)
{
    uint64_t lba = req->slba;
    struct ssdparams *spp = &ssd->sp;
    int len = req->nlb;
    uint64_t start_lpn = lba / spp->secs_per_pg;
    uint64_t end_lpn = (lba + len - 1) / spp->secs_per_pg;
   // struct ppa ppa;
    uint64_t lpn;
    uint64_t curlat = 0, maxlat = 0;

    if (end_lpn >= spp->tt_pgs) {
        return maxlat;
        femu_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
    }

    // 统计
    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        if (ssd->l_maptbl.o_maptbl[lpn] == 0) {
            ssd->l_maptbl.o_maptbl[lpn] = 1;
            ssd->l_maptbl.o_maptbl_cnt++;
            if (ssd->l_maptbl.o_maptbl_maxLBA < lpn) ssd->l_maptbl.o_maptbl_maxLBA = lpn;
            if (ssd->l_maptbl.o_maptbl_minLBA > lpn) ssd->l_maptbl.o_maptbl_minLBA = lpn;
        }
    }

    // 模拟写入并分配ppa然后计算闪存时延,
    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {

        ssd->WB[ssd->num_write_entries].LPA = lpn;
        ssd->WB[ssd->num_write_entries].PPA = 0;

        ssd->num_write_entries++;
        //缓冲区写满后 排序LPA然后分配PBA 之后学习索引
        if (ssd->num_write_entries == WB_Entries) {

            uint64_t stime = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);

            int n = ssd->num_write_entries;
            // sort lpn
            for (int i = 0; i < n - 1; i++) {
                for (int j = 0; j < n - i - 1; j++) {
                    if (ssd->WB[j].LPA > ssd->WB[j + 1].LPA) {
                        struct Write_Buffer temp;
                        temp = ssd->WB[j];
                        ssd->WB[j] = ssd->WB[j + 1];
                        ssd->WB[j + 1] = temp;
                    }
                }
            }
            // remove duplicate data
            int cur_idx = 0, uni_idx = 0;
            while (cur_idx < n) {
                ssd->WB[uni_idx] = ssd->WB[cur_idx];
                while (cur_idx < n - 1 && ssd->WB[cur_idx].LPA == ssd->WB[cur_idx + 1].LPA) {
                    cur_idx++;
                }
                uni_idx++;
                cur_idx++; 
            }
            n = uni_idx;

            struct ppa PPA;
            struct ppa TPPA[1030];
            for (int i = 0; i < n; i++) {
                PPA = get_new_page(ssd);
                TPPA[i] = PPA;
                uint64_t vpgidx = ppa2vpgidx(ssd, &PPA);
                ssd->WB[i].PPA = vpgidx;
                // For GC
                // set_rmap_ent(ssd, ssd->WB[i].LPA, &PPA);
                
                ssd->rmap[vpgidx] = ssd->WB[i].LPA;

                ssd_advance_write_pointer(ssd);

                ssd->l_maptbl.counter.group_write_cnt++;
            }

            Point *points = (Point *)(ssd->WB);
            
            if (0 && DeBUG) {
                femu_log("[ssd_lea_write]:排序后LPA-PPA信息\n");
                for (int i = 0; i < n; i++) {
                    femu_log("lpa:%u -> ppa:%u\n", points[i].x, points[i].y);
                }
            }

            FrameGroup_update(&ssd->l_maptbl, points, n);

            uint64_t etime = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
            femu_log("[training_cost]stime:%lu, etime:%lu, et-st:%lu ns\n", stime, etime, etime - stime);
           // FrameGroup_static(&ssd->l_maptbl);
            for (int i = 0; i < WB_Entries + 2; i++) {
                ssd->WB[i].LPA = 0;
                ssd->WB[i].PPA = 0;
            }
            ssd->num_write_entries = 0;

            // simulator lat            
            for (int i = 0; i < n; i++) {
                struct nand_cmd swr;
                swr.type = USER_IO;
                swr.cmd = NAND_WRITE;
                swr.stime = req->stime;
                /* get latency statistics */
                curlat = ssd_advance_status(ssd, &TPPA[i], &swr);
                maxlat = (curlat > maxlat) ? curlat : maxlat;
            }
           // maxlat = (n * NAND_PROG_LATENCY * 10000) / (spp->nchs * spp->luns_per_ch);
        }


        // 不满也刷 最后一次刷
        if (0 && ssd->flush) {
            // TODO
        }

    }     
    return maxlat;
}


static void *ftl_thread(void *arg)
{
    FemuCtrl *n = (FemuCtrl *)arg;
    struct ssd *ssd = n->ssd;
    NvmeRequest *req = NULL;
    uint64_t lat = 0;
    int rc;
    int i;

    while (!*(ssd->dataplane_started_ptr)) {
        usleep(100000);
    }

    /* FIXME: not safe, to handle ->to_ftl and ->to_poller gracefully */
    ssd->to_ftl = n->to_ftl;
    ssd->to_poller = n->to_poller;

    while (1) {
        for (i = 1; i <= n->num_poller; i++) {
            if (!ssd->to_ftl[i] || !femu_ring_count(ssd->to_ftl[i]))
                continue;

            rc = femu_ring_dequeue(ssd->to_ftl[i], (void *)&req, 1);
            if (rc != 1) {
                printf("FEMU: FTL to_ftl dequeue failed\n");
            }

            assert(req);
            
            switch (req->cmd.opcode) {
            case NVME_CMD_WRITE:
                lat = ssd_write(ssd, req);
                break;
            case NVME_CMD_READ:
                lat = ssd_read(ssd, req);
                break;
            case NVME_CMD_DSM:
                lat = 0;
                break;
            default:
                //femu_err("FTL received unkown request type, ERROR\n");
                ;
            }

            req->reqlat = lat;
            req->expire_time += lat;

            rc = femu_ring_enqueue(ssd->to_poller[i], (void *)&req, 1);
            if (rc != 1) {
                femu_err("FTL to_poller enqueue failed\n");
            }

            /* clean one line if needed (in the background) */
            if (should_gc(ssd)) {
                do_gc(ssd, false);
            }
        }
    }

    return NULL;
}

