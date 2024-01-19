#ifndef __FEMU_LeaFTL_H
#define __FEMU_LeaFTL_H

#include "../nvme.h"
#include <math.h>
#include <stdlib.h>

#define MB *1024*1024
#define KB *1024

#define MAX_CRB_SEGMENTS 256
#define ApproximateSegment Segment
#define Write_Buffer_size    4*1024*1024
#define Write_Buffer_Entries 1024
#define Mapping_TABLE_SIZE   8*1024*1024

#define L_MIN(a, b) ((a) < (b) ? (a) : (b))
#define L_MAX(a, b) ((a) > (b) ? (a) : (b))


typedef struct LRUCache {
    int size;
    int maxsize;
}LRUCache;   // todo 

typedef struct Point {
    uint32_t x;        // lpn
    uint32_t y;        // ppn
} Point;


// no use
typedef struct Points {
    Point points[256];
    int num_points;
    int split_id;
} Points;

typedef struct Split_Points{
    Points *s_points;
    int num_split;
    //bool *st;
    int *pos_split;
} Split_Points;

// 交接点
typedef struct InsecPoint {
    float x;        // lpn
    float y;        // ppn
} InsecPoint;


struct Avg_slope {
    unsigned int k : 15;         // 整数类型，占15位
    unsigned int indicator : 1;  // 指示位，占1位
};

typedef union {
    struct Avg_slope avg_slope;
    float value;
} Avg_slope_union;



typedef struct Bitmap {
    int length;
    unsigned char *bitmap;
} Bitmap;


typedef struct Segment {
    float  k;      // 2 Bytes  use Avg_slope ?
    double b;      // 4 Bytes
    uint8_t   x1;       // start_index 1 Bytes(256)  ?
    uint8_t   x2;       // end_inedex  1 Bytes 
    bool  accurate;       // True;  Segment Type 1:accurate, 2:approximate 
    bool  consecutive;      // 判断横坐标连续间隔是否相等
    Bitmap   filter;      // 位图判断

}Segment;

typedef struct SimpleSegment {
    float  k;       // 2B
    double b;       // 4B
    int   x1;       // start_index 1B
    int   x2;       // end_inedex  1B
}SimpleSegment;

typedef struct Segs {
    Segment *segments;
    int segment_id[20];
    int num_segments;
} Segs;

enum {
    PLR_CONSTANTS_FIRST = 1,
    PLR_CONSTANTS_SECOND = 2,
    PLR_CONSTANTS_READY = 3
};

typedef struct PLR{
    double gamma;
    int max_length;  // 256
    // vector<Segment> segments;
    Segment *segments;
    uint8_t num_segments;

    Point s0;
    Point s1;
    SimpleSegment rho_upper;
    SimpleSegment rho_lower;
    InsecPoint sint;
    int state;

    Point  *points;
    int num_points;  

} PLR;




/**
 * @brief 日志结构学习索引段 LogPLR 层级结构体
 * 
 */
typedef struct LogPLR{
    double gamma;
    int max_length;  // 256
    // vector<Segment> segments;
    uint8_t  level;

    Segment *segments;
    uint8_t num_segments;

} LogPLR;

typedef struct Group    {
    double gamma;
    PLR plr;    // 专门用来学习段,之后插日志结构LogPLR, (用来学习索引段的临时数据结构，临时工)

    LogPLR *L;  
    int num_levels;

    int group_id;

    // method
} Group;

/**
 * @brief 学习索引段Frame包装结构 mapping table
 * 
 */

typedef struct Counter {
    int group_write_miss;
    int group_write_hit;
}Counter;


typedef struct FrameGroup {

     Group *groups;
     int num_groups;

     double gamma;
     int num_segments;

     int max_size;
     int frame_length;
     struct Counter counter; 

    // method

} FrameGroup;



// typedef struct FlashMetadata
// {
//     int counter;
//     int flash_num_blocks;
//     int flash_npage_per_block;
    
//     // lmapping table
//     double gamma;
//     //self.mapping_table = LogPLR(frame_no=0, gamma=self.gamma)
//     FrameGroup l_maptbl;
    
//     // @method
// };


typedef struct CRB{
    ApproximateSegment segments[MAX_CRB_SEGMENTS];  
    uint16_t num_segments;  // 256
} CRB;



void crb_init(CRB *crb);
void crb_insert_segment(CRB *crb, Segment* seg);
unsigned int crb_find_segment(CRB *crb, uint8_t lpa);

bool Segment_is_valid(Segment *seg, uint32_t x);
uint32_t Segment_gety(Segment *seg, uint32_t x);
int isConsecutive(Point *points, int num_points);
void Segment_check_properties(Segment *seg, Point *points, int num_points);
void Segment_init(Segment *seg, float k, double b, int x1, int x2, Point *points, int num_points);
bool Segment_overlaps(Segment *seg1, Segment *seg2);
void Segment_merge(Segment *new_seg, Segment *old_seg, int *samelevel);
void SimpleSegment_init(SimpleSegment *simpseg, double k, double b, int x1, int x2);
int get_y(SimpleSegment *simpleseg, int x);

Point intersection(SimpleSegment *s1, SimpleSegment *s2);
InsecPoint inter_section(SimpleSegment *s1, SimpleSegment *s2);
bool is_above(Point *pt, SimpleSegment *s);
bool is_below(Point *pt, SimpleSegment *s);
Point get_upper_bound(Point *pt, double gamma);
Point get_lower_bound(Point *pt, double gamma);
SimpleSegment frompoints(Point p1, Point p2);
SimpleSegment frompoints_insec(InsecPoint p1, Point p2);

void plr_init(PLR* plr, double gamma);
void plr__init(PLR* plr);
void plr_destroy(PLR* plr);
void plr_add_segment(PLR *plr, Segment *seg);

int build_segment(PLR* plr, Segment *seg);
bool should_stop(PLR* plr, Point *point);
int process_point(PLR* plr, Point* point, Segment *seg);
void plr_learn(PLR* plr, Point* points, int num_points);
void plr_sorted(PLR* plr);

void LogPLR_init(LogPLR *logplr, int level);
uint8_t binary_search(LogPLR *logplr, Segment *seg);
void LogPLR_add_segment(LogPLR *logplr, Segment *seg, int *index);
void LogPLR_del_segment(LogPLR *logplr, int pos);

void Segs_add_segment(Segs *segs, Segment *seg, int seg_id);

void Group_init(Group *group, double gamma, int group_id);
void Group_lookup(void);
void Group_seg_merge(void);
void Group_gc(void);
void Group_add_LogPLR(Group *group);
void Group_del_segments(Group *group, int level, int pos);
void Group_add_segments(Group *group, int level, Segment *segments, int num_segments, bool recursive);
void Group_update(Group *group, Point* points, int num_points);
LRUCache *createLRUCache(int maxsize);
void FrameGroup_init(FrameGroup *framegroup, double gamma);
void FrameGroup_static(FrameGroup *framegroup);
void FrameGroup_update(FrameGroup *framegroup, Point* points, int num_points);




// /**
//  * @brief CRB_method
//  * 
//  * 还未用上 
//  */

// void crb_init(CRB *crb) {
//     crb->num_segments = 0;
// }

// void crb_insert_segment(CRB *crb, Segment* seg) {
//     uint8_t start_lpa = seg->x1;
//     // Check if the CRB is already full
//     if (crb->num_segments >= MAX_CRB_SEGMENTS) {
//         femu_log("CRB is full. Cannot insert segment.\n");
//         return;
//     }

//     // Check if the start_lpa already exists in the CRB
//     for (unsigned char i = 0; i < crb->num_segments; i++) {
//         if (crb->segments[i].x1 == start_lpa) {
//             femu_log("Segment with start LPA %u already exists. Updating start LPA.\n", start_lpa);
//             // 这一步应该是更新成往右相邻的点的LPA, 但是Segment的信息只有s,l,k,b并不知道下个点是啥有点迷
//             crb->segments[i].x1++;  
//             return;
//         }
//     }

//     // Insert the new segment
//     crb->segments[crb->num_segments].x1 = start_lpa;
//     crb->segments[crb->num_segments].x2 = seg->x2;
//     crb->num_segments++;

//     // Sort the segments by their starting LPA using bubble sort
//     // 冒泡排序
//     for (unsigned char i = 0; i < crb->num_segments - 1; i++) {
//         for (unsigned char j = 0; j < crb->num_segments - i - 1; j++) {
//             if (crb->segments[j].x1 > crb->segments[j + 1].x1) {
//                 // Swap the segments
//                 ApproximateSegment temp = crb->segments[j];
//                 crb->segments[j] = crb->segments[j + 1];
//                 crb->segments[j + 1] = temp;
//             }
//         }
//     }
// }

// // 找到属于哪个近似段或者是该近似段的起始LPA
// unsigned int crb_find_segment(CRB *crb, uint8_t lpa) {
//     // Binary search to find the segment containing the given LPA
//     unsigned int left = 0;
//     unsigned int right = crb->num_segments - 1;

//     while (left <= right) {
//         unsigned int mid = left + (right - left) / 2;

//         if (crb->segments[mid].x1 <= lpa && lpa < crb->segments[mid].x2) {
//             // Found the segment containing the LPA
//             // return crb->segments[mid].x1;
//             return mid;
//         }
//         else if (crb->segments[mid].x1 > lpa) {
//             // LPA is in the left half
//             right = mid - 1;
//         }
//         else {
//             // LPA is in the right half
//             left = mid + 1;
//         }
//     }

//     // LPA is not found in any segment
//     return -1;
// }

// /**
//  * @brief Segment_method
//  * 
//  */


// // 可能出现bug
// bool Segment_is_valid(Segment *seg, uint32_t x) {
//     if (!(x >= seg->x1 && x <= seg->x2))
//         return false;
//     if (seg->accurate) {
//         int k = round(1.0/seg->k)   ; 
//         if (((x - seg->x1) % k) != 0)        // 这里可能很多不通过
//             return false;
//     }
//     // 暂时未添加 判断 x 是否在段上的条件
//     // else {
//     //     if (seg->filter.length > 0 && (x - seg->x1) <= seg->filter.length) 
//     //         return seg->filter.bitmap[x - seg->x1];
//     // }
//     return true;
// }


// uint32_t Segment_gety(Segment *seg, uint32_t x) {
//     int predict = 0;
//     if (Segment_is_valid(seg, x)) {
//             predict = round(x*seg->k + seg->b);
//             return predict;
//     }
//     return -1;
// }

// int isConsecutive(Point *points, int num_points) {
//     int i;
//     int diff = points[1].x - points[0].x;  // 计算第一个相邻点的 x 坐标差异
//     for (i = 2; i < num_points; i++) {
//         if (points[i].x - points[i-1].x != diff) {
//             return 0;  // 存在不连续的 x 坐标差异
//         }
//     }

//     return 1;  // 所有的 x 坐标差异都相等，即连续
// }


// // 判断段属性: 精确性，连续性
// void Segment_check_properties(Segment *seg, Point *points, int num_points) {
//     seg->accurate = true, seg->consecutive = true;
    
//     for (int i = 0; i < num_points; i++) {
//         uint32_t ppa;
//         ppa = Segment_gety(seg, points[i].x);
//         if (ppa != points[i].y) {
//             seg->accurate = false;
//         }
//         if (!isConsecutive(points, num_points)) {
//             seg->consecutive = false;
//         }

//         if (!seg->accurate && !seg->consecutive) return ;
//     }
// }


// void Segment_init(Segment *seg, float k, double b, int x1, int x2, Point *points, int num_points) {
//     seg->k = k;
//     seg->b = b;
//     seg->x1 = x1;
//     seg->x2 = x2;
//     seg->accurate = true;
      
//     seg->filter.length = 0;
//     seg->filter.bitmap = NULL;


//     if (points != NULL) {
//         Segment_check_properties(seg, points, num_points);

//         if (!seg->consecutive) {
//             // 初始化位图
//             if (seg->filter.bitmap != NULL) {
//                 free(seg->filter.bitmap);
//                 seg->filter.length = 0;
//                 seg->filter.bitmap = NULL;
//             }
//             seg->filter.length = seg->x2 - seg->x1 + 1;
//             seg->filter.bitmap = (unsigned char *)malloc(seg->filter.length * sizeof(unsigned char));
//             memset(seg->filter.bitmap, 0, seg->filter.length * sizeof(unsigned char));
            
//             if (seg->filter.bitmap != NULL) {
//                 memset(seg->filter.bitmap, 0, seg->filter.length * sizeof(unsigned char));
            
//                 for (int i = 0; i < num_points; i++) {
//                     seg->filter.bitmap[points[i].x - seg->x1] = 1;
//                 }
//             }
//         }
//     }
// }

// bool Segment_overlaps(Segment *seg1, Segment *seg2) {
//     return L_MIN(seg1->x2, seg2->x2) - L_MAX(seg1->x1, seg2->x1) >= 0;
// }


// // may has a bug
// void Segment_merge(Segment *new_seg, Segment *old_seg, int *samelevel) {
//     //不重叠
//     if (new_seg->x1 >= old_seg->x2 || new_seg->x2 <= old_seg->x1) {
//         *samelevel = 1;  
//         return ;
//     } 
//     //else if (new_seg->consecutive)  no use new_seg->consecutive

//     int start = L_MIN(new_seg->x1, old_seg->x1);
//     int end   = L_MAX(new_seg->x2, new_seg->x2);

//     Bitmap new_bm, old_bm;
//     new_bm.length = end - start + 1;
//     old_bm.length = end - start + 1;
//     new_bm.bitmap = (unsigned char *)malloc(new_bm.length * sizeof(unsigned char));
//     old_bm.bitmap = (unsigned char *)malloc(old_bm.length * sizeof(unsigned char));
//     memset(new_bm.bitmap, 0, new_bm.length * sizeof(unsigned char)); 
//     memset(old_bm.bitmap, 0, old_bm.length * sizeof(unsigned char));

//     int j = 0;
//     for (int i = new_seg->x1 - start; i < new_seg->x2 - start + 1; i++) {
//         new_bm.bitmap[i] = new_seg->filter.bitmap[j];
//         j++;
//     }
//     j = 0;
//     for (int i = old_seg->x1 - start; i < old_seg->x2 - start + 1; i++) {
//         old_bm.bitmap[i] = old_seg->filter.bitmap[j];
//         j++;
//     }
//     j = 0;
 
//     for (int i = 0; i < end - start + 1; i++) {
//         old_bm.bitmap[i] = old_bm.bitmap[i] & (~new_bm.bitmap[i]);
//     }
    
//     // Debug
//     if (0) {
//         femu_log("new_bm:");
//         for (int i = 0; i < end - start + 1; i++) {
//             printf(" %u", new_bm.bitmap[i]);
//         }
//         printf("\n");
//         femu_log("old_bm:");
//         for (int i = 0; i < end - start + 1; i++) {
//             printf(" %u", old_bm.bitmap[i]);
//         }
//         printf("\n");
//     }

//     int first_valid = -1;
//     for (int i = 0; i < end - start + 1; i++) {
//         if (old_bm.bitmap[i]) {
//             first_valid = i;
//             break;
//         }
//     }
//     if (first_valid == -1) {
//         *samelevel = -1;
//         return ;
//     }
//     int last_valid =  -1;
//     for (int i = end - start; i >= 0; i--) {
//         if (old_bm.bitmap[i]) {
//             last_valid = i;
//             break;
//         }
//     }
//     old_seg->x1 = first_valid + start;
//     old_seg->x2 = last_valid  + start;  // 存在bug,在这里对old_seg处理,只是在overlap_segs中的old_seg做处理,
//                                         // 并没有对原来的logplr中的segs中的old_seg做处理


//     for (int i = first_valid; i < last_valid + 1; i++) {
//         old_seg->filter.bitmap[j] =  old_bm.bitmap[i];
//         j++;
//     }

//     if (new_seg->x1 >= old_seg->x2 || new_seg->x2 <= old_seg->x1) {
//         *samelevel = 1;  // 存在不重叠
//         return ;
//     }     

//     *samelevel = 0;   // 存在仍重叠 
//     return ;

// }



// /**
//  * @brief SimpleSegment_method
//  * 
//  */
// void SimpleSegment_init(SimpleSegment *simpseg, double k, double b, int x1, int x2) {
//     simpseg->k = k;
//     simpseg->b = b;
//     simpseg->x1 = x1;
//     simpseg->x2 = x2;    
// }

// // int round(double x) {
// //     return (int)(x + 0.5);
// // }

// int get_y(SimpleSegment *simpleseg, int x) {
//     int predict;
//     predict = round(x * simpleseg->k + simpleseg->b);
//     return predict;
// }

// Point intersection(SimpleSegment *s1, SimpleSegment *s2) {
//     Point p;
//     p.x = (int) ((s2->b - s1->b) / (s1->k - s2->k));
//     p.y = (int) ((s1->k * s2->b - s2->k * s1->b) / (s1->k - s2->k));
//     return p;
// }

// InsecPoint inter_section(SimpleSegment *s1, SimpleSegment *s2) {
//     InsecPoint insec_pt;
//     insec_pt.x = (float) ((s2->b - s1->b) / (s1->k - s2->k));
//     insec_pt.y = (float) ((s1->k * s2->b - s2->k * s1->b) / (s1->k - s2->k));
//     return insec_pt;
// }

// bool is_above(Point *pt, SimpleSegment *s) {
//     return pt->y > (int)(s->k * pt->x + s->b);
// }

// bool is_below(Point *pt, SimpleSegment *s) {
//     return pt->y < (int)(s->k * pt->x + s->b);
// }

// Point get_upper_bound(Point *pt, double gamma) {
//     Point p;
//     p.x = pt->x, p.y = pt->y + gamma;
//     return p;
// }

// Point get_lower_bound(Point *pt, double gamma) {
//     Point p;
//     p.x = pt->x, p.y = pt->y - gamma;
//     return p;
// }

// SimpleSegment frompoints(Point p1, Point p2) {
//     SimpleSegment simplesegment;
//     simplesegment.k = (float)((p2.y - p1.y) / (p2.x - p1.x));
//     simplesegment.b = -simplesegment.k * p1.x + p1.y;
//     simplesegment.x1 = p1.x;
//     simplesegment.x2 = p2.x;
//     return simplesegment;
// }

// SimpleSegment frompoints_insec(InsecPoint p1, Point p2) {
//     SimpleSegment simplesegment;
//     simplesegment.k = (float)((p2.y - p1.y) / (p2.x - p1.x));
//     simplesegment.b = -simplesegment.k * p1.x + p1.y;
//     simplesegment.x1 = p1.x;
//     simplesegment.x2 = p2.x;
//     return simplesegment;
// }

// // typedef struct {
// //     const char* FIRST;
// //     const char* SECOND;
// //     const char* READY;
// // } PLR_CONSTANTS;



// /**
//  * @brief PLR_method
//  * 
//  */
// void plr_init(PLR* plr, double gamma) {
//     plr->gamma = gamma;
//     plr->max_length = 256;

//     plr__init(plr);
// }


// // Tem_struct, use for build segments, and then destry after insert LSM struc
// void plr__init(PLR* plr) {

//     plr_destroy(plr);

//     plr->segments = g_malloc0(sizeof(Segment) * (plr->max_length));

//     plr->s0.x = 0, plr->s0.y = 0;   // init point
//     plr->s1.x = 0, plr->s1.y = 0;
//     plr->rho_upper.x1 = 0, plr->rho_upper.x2 = 0;
//     plr->rho_lower.x1 = 0, plr->rho_lower.x2 = 0;
//     plr->sint.x = 0, plr->sint.y = 0;
//     plr->state = PLR_CONSTANTS_FIRST;
    
//     plr->points = NULL;
//     plr->num_points = 0;

// }


// void plr_add_segment(PLR *plr, Segment *seg) {
//     if (plr->num_segments == 0) {
//         plr->segments = (Segment *)malloc(1 * sizeof(Segment));
//         plr->segments[0] = *seg;  // C语言结构浅拷贝,有可能出bug
//         plr->num_segments++; 
//     } else {
//         plr->segments = (Segment *)realloc(plr->segments, (plr->num_segments + 1) * sizeof(Segment));
//         plr->segments[plr->num_segments] = *seg;
//         plr->num_segments++;
//     }    
// }

// void plr_destroy(PLR* plr) {
//     if (plr->segments != NULL) {
//         free(plr->segments);
//         plr->segments = NULL;
//         plr->num_segments = 0;
//     }
//     if (plr->points != NULL) {
//         free(plr->points);
//         plr->points = NULL;
//         plr->num_points = 0;        
//     }
//    // free(plr);
// }

// int build_segment(PLR* plr, Segment *seg) {
//         if (plr->state == PLR_CONSTANTS_FIRST) { 
//             seg = NULL;
//             return 0;
//         }
//         // 建立单点段
//         else if (plr->state == PLR_CONSTANTS_SECOND) {
//             Segment_init(seg, 1, plr->s0.y - plr->s0.x, plr->s0.x, plr->s0.x, plr->points, plr->num_points);
//         }

//         // 建立多点段
//         else if (plr->state == PLR_CONSTANTS_READY) {
//             float avg_slope = (float16)((plr->rho_lower.k + plr->rho_upper.k) / 2.0);
//             double intercept = -plr->sint.x * avg_slope + plr->sint.y;
//             Segment_init(seg, avg_slope, intercept, plr->s0.x, plr->s1.x, plr->points, plr->num_points);
//         }
//         return 1;
// }

// bool should_stop(PLR* plr, Point *point) {
//     if (plr->s1.x == 0) {
//         if (point->x > plr->s0.x + plr->max_length) return true;
//     }else if (point->x > plr->s1.x + plr->max_length) return true;

//     return false;
// }


// // 处理点的函数,up and lower
// int process_point(PLR* plr, Point* point, Segment *seg) {
//     int ret = 0;
//     if (plr->state == PLR_CONSTANTS_FIRST) {
//         plr->s0 = *point;
//         plr->state = PLR_CONSTANTS_SECOND;
//     } else if (plr->state == PLR_CONSTANTS_SECOND) {
//         if (should_stop(plr, point)) {
//             ret = build_segment(plr, seg);

//             plr->s0 = *point;
//             plr->state = PLR_CONSTANTS_SECOND;

//             free(plr->points);
//             plr->points = NULL;
//             plr->num_points = 0;
//         }else{
//             plr->s1 = *point;
//             // rho_lower 和 rho_upper 和 sint 还不知道作用是啥
//             plr->rho_lower = frompoints(get_upper_bound(&plr->s0, plr->gamma), 
//                                                     get_lower_bound(&plr->s1, plr->gamma));
//             plr->rho_upper = frompoints(get_lower_bound(&plr->s0, plr->gamma), 
//                                                     get_upper_bound(&plr->s1, plr->gamma));            
//             plr->sint = inter_section(&plr->rho_upper, &plr->rho_lower);

//             plr->state = PLR_CONSTANTS_READY;
//         }
//     } else if (plr->state == PLR_CONSTANTS_READY) {
//         if (!is_above(point, &plr->rho_lower) || !is_below(point, &plr->rho_upper) || should_stop(plr, point)) {
//             ret = build_segment(plr, seg);
//             plr->s0 = *point;
//             plr->state = PLR_CONSTANTS_SECOND;

//             free(plr->points);
//             plr->points = NULL;
//             plr->num_points = 0;
//         }else {
//             // QS
//             plr->s1 = *point;

//             Point s_upper = get_upper_bound(point, plr->gamma);
//             Point s_lower = get_lower_bound(point, plr->gamma);

//             if (is_below(&s_upper, &plr->rho_upper)) 
//                 plr->rho_upper = frompoints_insec(plr->sint, s_upper);
//             if (is_above(&s_lower, &plr->rho_lower)) 
//                 plr->rho_lower = frompoints_insec(plr->sint, s_lower);
//         }
//     }

//     if (plr->num_points == 0) {
//             plr->points = malloc(sizeof(Point));
//         } else {
//             plr->points = realloc(plr->points, (plr->num_points + 1) * sizeof(Point));
//         }
//         plr->points[plr->num_points] = *point;
//         plr->num_points++;
    
//     return ret;
// }

// // 学习得到一堆学习索引段 存在PLR中 还没有日志结构
// void plr_learn(PLR* plr, Point* points, int num_points) {
//     // plr_init(plr);
//     //Segment* rejs = malloc(num_points * sizeof(Segment));
//     //int num_rejs = 0;



//     for (int i = 0; i < num_points; i++) {
//         Segment seg;
//         int ret = process_point(plr, &points[i], &seg);
//         if (ret) {
//             femu_log("[plr_learn]: 学习成功后添加段\n");
//             plr_add_segment(plr, &seg);
//         }
//     }

//     Segment final_seg;
//     int ret = build_segment(plr, &final_seg);
//     if (ret) {
//             femu_log("[plr_learn]: 学习成功后添加段(最后的段)\n");
//             plr_add_segment(plr, &final_seg);
//     }

//    // *num_segments = plr->num_segments;
//     return ;
// }

// // TODO
// void plr_sorted(PLR* plr) {

// }


// /**
//  * @brief LogPLR_method
//  * 
//  * 有关日志结构索引段的方法函数
//  * 
//  * @param logplr 输入日志结构学习索引段
//  * @param other 
//  * 
//  * @return void
//  * 
//  * @details 
//  * LogPLR_init();
//  * binary_search();
//  * LogPLR_add_segment();
//  * LogPLR_add_segment
//  * LogPLR_del_segment
//  * Segs_add_segment();   // Segs's method
//  */



// void LogPLR_init(LogPLR *logplr, int level) {
//     logplr->level = level;
//     logplr->segments = NULL;
//     logplr->num_segments = 0;
// }


// uint8_t binary_search(LogPLR *logplr, Segment *seg) {
//     uint8_t left = 0;
//     uint8_t right = logplr->num_segments - 1;
//     uint8_t mid = 0;
//     uint8_t target = seg->x1;
//     while (left <= right) {
//         mid = left + (right - left) / 2;

//         if (logplr->segments[mid].x1 == target) {
//             // 如果找到目标元素，则直接返回其索引
//             return mid;
//         } else if (logplr->segments[mid].x1 < target) {
//             // 如果目标元素大于中间元素，则在右半部分继续查找
//             left = mid + 1;
//         } else {
//             right = mid - 1;
//         }
//     }
//     // 当循环结束时，说明目标元素在数组中不存在
//     // 此时，left表示应该插入的位置
//     return left;
// }

// // L_i添加段
// // 如果释放掉PLR之后，LogPLR中的段会不会也不见？涉及到C语言结构的浅拷贝
// void LogPLR_add_segment(LogPLR *logplr, Segment *seg, int *index) {
//     if (logplr->num_segments == 0) {
//         logplr->segments = (Segment *)malloc(1 * sizeof(Segment));
//         logplr->segments[0] = *seg;  // C语言结构浅拷贝,有可能出bug 
//         logplr->num_segments++; 
//         *index = 0;
//     } else {
//         logplr->segments = (Segment *)realloc(logplr->segments, (logplr->num_segments + 1) * sizeof(Segment));
//         // 二分查找插入位置
//         uint8_t pos = binary_search(logplr, seg);
//         *index = pos;
//         logplr->num_segments++;
//         for (int i = logplr->num_segments - 1; i >= pos + 1; i++) {
//             logplr->segments[i] = logplr->segments[i - 1];
//         }
//         logplr->segments[pos] = *seg;
//     }
// }

// void LogPLR_del_segment(LogPLR *logplr, int pos) {
//     if (pos < 0 || pos >= logplr->num_segments) {
//         femu_log("删除位置不对\n");
//         return ;
//     }
//     for (int i = pos; i < logplr->num_segments - 1; i++) {
//         logplr->segments[i] = logplr->segments[i + 1];
//     }
//     logplr->num_segments--;
//     logplr->segments = (Segment *)realloc(logplr->segments, logplr->num_segments * sizeof(Segment));

//     if (logplr->num_segments == 0) {
//         free(logplr->segments);
//         logplr->segments = NULL;
//     }
// }


// void Segs_add_segment(Segs *segs, Segment *seg, int seg_id) {
//     if (segs->num_segments == 0) {
//         segs->segments = (Segment *)malloc(1 * sizeof(Segment));
//         segs->segments[0] = *seg;  // C语言结构浅拷贝,有可能出bug
//         segs->segment_id[0] = seg_id; 
//         segs->num_segments++; 
//     } else {
//         segs->segments = (Segment *)realloc(segs->segments, (segs->num_segments + 1) * sizeof(Segment));
//         segs->segments[segs->num_segments] = *seg;
//         segs->segment_id[segs->num_segments] = seg_id;
//         segs->num_segments++;
//     }
// }




// /**
//  * @brief Group_method
//  * 
//  */


// void Group_lookup(void) {
    
// }

// void Group_seg_merge(void) {
    
// }

// void Group_gc(void) {

// }

// void Group_init(Group *group, double gamma, int group_id) {
//     plr_init(&group->plr, gamma);
//     group->L = NULL;
//     group->num_levels = 0;

//     group->group_id = group_id;
// }

// void Group_add_LogPLR(Group *group) {
//     if (group->num_levels == 0) {
//         group->L = (LogPLR *)malloc(1 * sizeof(LogPLR));
//         LogPLR_init(&group->L[group->num_levels], group->num_levels);
//         group->num_levels++;
//     } else {
//         group->L = (LogPLR *)realloc(group->L, (group->num_levels + 1) * sizeof(LogPLR));
//         LogPLR_init(&group->L[group->num_levels], group->num_levels);
//         group->num_levels++;        
//     }
// }

// void Group_del_segments(Group *group, int level, int pos) {
//     LogPLR_del_segment(&group->L[level], pos);
//     return ;
// }

// void Group_add_segments(Group *group, int level, Segment *segments, int num_segments, bool recursive) {
//     while (group->num_levels <= level) {
//         Group_add_LogPLR(group);
//     }

//    // LogPLR *Li = &group->L[level];
//     Segs confict_segs;
//     confict_segs.num_segments = 0;
//     confict_segs.segments = NULL;

//     for (int i = 0; i < num_segments; i++) {
//      //   Segment *seg = &segments[i];
//         femu_log("segment:[%d] add to level: %d\n", i, level);
//         int index = 0;
//         if (group->L[level].num_segments == 0) {
//             LogPLR_add_segment(&group->L[level], &segments[i], &index);
//             continue;
//         }
//         Segs overlap_segs;
//         overlap_segs.num_segments = 0;
//         overlap_segs.segments = NULL;

//         LogPLR_add_segment(&group->L[level], &segments[i], &index);
//         if (index != 0) {
//             // 在前面的最多只有一个重叠
//             Segs_add_segment(&overlap_segs, &group->L[level].segments[index - 1], index - 1);
//         }
//         for (int j = index + 1; j < group->L[level].num_segments; j++) {
//             // 添加后续的重叠区间段
//             if (group->L[level].segments[j].x1 >= segments[i].x2) {
//                 break;
//             }
//             Segs_add_segment(&overlap_segs, &group->L[level].segments[j], j);            
//         }

//         uint8_t indices_to_delete[50];
//         int     indect_pointer = 0;
//         for (int j = 0; j < overlap_segs.num_segments; j++) {
//             int same_level = 1;
//             // 段合并 same_level 判断合并操作后的旧段状态 1: 同一层不重叠, 0:同一层重叠, -1:旧段不存在，被新段完全包含 
//             Segment_merge(&segments[i], &overlap_segs.segments[j], &same_level);
//             if (same_level == -1) {
//                 indices_to_delete[indect_pointer++] = overlap_segs.segment_id[j]; 
//             } else if (same_level == 0) {
//                 Segs_add_segment(&confict_segs, &overlap_segs.segments[j], overlap_segs.segment_id[j]);
//                 indices_to_delete[indect_pointer++] = overlap_segs.segment_id[j];                   
//             }

//         }

                
//         // TODO
//         // 插入时overlaps 处理区间重叠矛盾 删除
//         for (int j = 0; j < indect_pointer; j++) {
//             Group_del_segments(group, level, indices_to_delete[j]);
//             // del group->L[level].segment[j];
//         }

//     }

//     // Todo 如果添加段时出现区间重叠，则将旧段往下层推，如果在下层还是重叠，则为该段创建新的Level
//     if (recursive) {
//         if (confict_segs.num_segments != 0) {
//             Group_add_segments(group, level + 1, confict_segs.segments, confict_segs.num_segments, false);
//         }
//     } else {
//         if (confict_segs.num_segments != 0) {
//             Group_add_LogPLR(group);
//             for (int i = 0; i < confict_segs.num_segments; i++) {
//                 int no_use = 0;
//                 LogPLR_add_segment(&group->L[group->num_levels - 1], &confict_segs.segments[i], &no_use); // 表示用不上id                
//             }       
//         }
//     }
//     return ;

// }

// void Group_update(Group *group, int level, Point* points, int num_points) {
//     // Points 应该是已经排序好的 lpa : ppa
//     //double gamma = 0;

//     plr__init(&group->plr);
//     plr_learn(&group->plr, points, num_points);
//     plr_sorted(&group->plr);
    
//     Group_add_segments(group, 0, group->plr.segments, group->plr.num_segments, false);  // should be true
     
//     plr_destroy(&group->plr);

// }

// // no use
// LRUCache *createLRUCache(int maxsize) {
//     LRUCache* cache = g_malloc0(sizeof(LRUCache));
//     cache->maxsize = maxsize;
//     return cache;  
// }


// void FrameGroup_init(FrameGroup *framegroup, double gamma) {
//     framegroup->gamma = gamma;
//     framegroup->num_segments = 0;
//     framegroup->frame_length = 256;  // maybe不够
//     framegroup->max_size = Mapping_TABLE_SIZE;
    
//     framegroup->counter.group_write_hit = 0;
//     framegroup->counter.group_write_miss = 0;
    
//     framegroup->num_groups = framegroup->frame_length;
//     framegroup->groups = (Group *)malloc(framegroup->num_groups * sizeof(Group));
//     for (int i = 0; i < framegroup->num_groups; i++) {
//         Group_init(&framegroup->groups[i], framegroup->gamma, i);
//     }

//     return ;
// }



#endif