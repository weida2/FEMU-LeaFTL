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
#define DeBUG 1

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
    double x;        // lpn
    double y;        // ppn
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
    double k;      // 2 Bytes  use Avg_slope ? for more accurate
    double b;      // 4 Bytes
    uint8_t   x1;       // start_index 1 Bytes(256)  ?
    uint8_t   x2;       // end_inedex  1 Bytes 
    bool  accurate;       // True;  Segment Type 1:accurate, 2:approximate 
    bool  consecutive;      // 判断横坐标连续间隔是否相等
    Bitmap   filter;      // 位图判断

}Segment;

typedef struct SimpleSegment {
    double k;       // 2B 
    double b;       // 4B
    int   x1;       // start_index 1B
    int   x2;       // end_inedex  1B
}SimpleSegment;

typedef struct Segs {
    Segment *segments;
    int segment_id[20];
    int num_segments;
} Segs;

typedef struct Segs_Simple {
    Segment *segments;
    int num_segments;
} Segs_Simple;


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
    int group_read_cnt;
    int group_read_acc_hit;
    int group_reaa_noacc_hit;
    int group_read_noacc_miss;
    int group_read_miss;
}Counter;


typedef struct FrameGroup {

     Group *groups;
     int num_groups;
     int cnt_groups;

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
    uint8_t num_segments;  // 256
} CRB;



void crb_init(CRB *crb);
void crb_insert_segment(CRB *crb, Segment* seg);
unsigned int crb_find_segment(CRB *crb, uint8_t lpa);

bool Segment_is_valid(Segment *seg, uint32_t x);
uint32_t Segment_gety(Segment *seg, bool check, uint32_t x);
int isConsecutive(Point *points, int num_points);
void Segment_check_properties(Segment *seg, Point *points, int num_points);
void Segment_init(Segment *seg, double k, double b, int x1, int x2, Point *points, int num_points);
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
int8_t binary_search(LogPLR *logplr, Segment *seg);
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
uint64_t FrameGroup_lookup(FrameGroup *framegroup, uint64_t lpn, Segment *seg);
uint64_t PPN2VPPN(uint64_t pgidx);


#endif