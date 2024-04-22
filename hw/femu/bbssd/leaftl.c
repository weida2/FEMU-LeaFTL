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
                ApproximateSegment temp = crb->segments[j];
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
    if (!old_seg || !old_seg->k || !old_seg->x1 || !old_seg->filter.bitmap) {
        *samelevel = 2;
        return;
    }

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



/**
 * @brief SimpleSegment_method
 * 
 */
void SimpleSegment_init(SimpleSegment *simpseg, double k, double b, int x1, int x2) {
    simpseg->k = k;
    simpseg->b = b;
    simpseg->x1 = x1;
    simpseg->x2 = x2;    
}

// int round(double x) {
//     return (int)(x + 0.5);
// }

int get_y(SimpleSegment *simpleseg, int x) {
    int predict;
    predict = round(x * simpleseg->k + simpleseg->b);
    return predict;
}

Point intersection(SimpleSegment *s1, SimpleSegment *s2) {
    Point p;
    p.x = (int) ((s2->b - s1->b) / (s1->k - s2->k));
    p.y = (int) ((s1->k * s2->b - s2->k * s1->b) / (s1->k - s2->k));
    return p;
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
    //double gamma = 0;

    // femu_log("[Group_update]: in\n");
    plr__init(&group->plr);
    plr_learn(&group->plr, points, num_points);
    plr_sorted(&group->plr);
    
   // femu_log("[Group_update]: plr学习完成, num_segments:%d\n", group->plr.num_segments);
    
    // print
    if (0 && DeBUG) {
        for (int i = 0; i < group->plr.num_segments; i++) {
            Segment seg = group->plr.segments[i];
            femu_log("[Group_update][plr_Segment %d]: k: %.2f, b: %.2f, x1: %u, x2: %u, %s, %s\n", i,
                seg.k, seg.b, seg.x1, seg.x2, seg.accurate ? "精确":"不精确", seg.consecutive ? "连续":"不连续");
        }
    }
    Group_add_segments(group, 0, group->plr.segments, group->plr.num_segments, true);  // should be true

    
    if (0 && DeBUG) {
        femu_log("[Group_update]: 此次segs索引段添加到group完成, 此时的group[%d]信息如下:\n", group->group_id);
        for (int i = 0; i < group->num_levels; i++) {
            for (int j = 0; j < group->L[i].num_segments; j++) {
                Segment seg = group->L[i].segments[j];
                femu_log("[Group_update] group[%d]->L[%d][Segment %d]: k: %.2f, b: %.2f, x1: %u, x2: %u, %s, %s\n", group->group_id, 
                i, j, seg.k, seg.b, seg.x1, seg.x2, seg.accurate ? "精确":"不精确", seg.consecutive ? "连续":"不连续");            
            }
        }
    }
    plr_destroy(&group->plr);
    // femu_log("[Group_update]: out\n");
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
    framegroup->num_groups = 16 * 4 * 1024;  // 不太对上面那样算 是4G的条目 下面是按照8MB的映射表大小计算 

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