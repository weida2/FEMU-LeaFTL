#ifndef __FEMU_LeaFTL_H
#define __FEMU_LeaFTL_H

#include <math.h>
#include <stdlib.h>
#include "../nvme.h"
#include "plr.h"

#define MB *1024*1024
#define KB *1024

#define MAX_CRB_SEGMENTS 256
#define Write_Buffer_size    4*1024*1024
#define Write_Buffer_Entries 1024
#define Mapping_TABLE_SIZE   8*1024*1024

#define L_MIN(a, b) ((a) < (b) ? (a) : (b))
#define L_MAX(a, b) ((a) > (b) ? (a) : (b))
#define DeBUG 1

// typedef struct LRUCache {
//     int size;
//     int maxsize;
// }LRUCache;   // todo 



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
    int num_segments;

} LogPLR;

typedef struct Group    {
    double gamma;
    PLR plr;                    // 专门用来学习段,之后插日志结构LogPLR, (用来学习索引段的临时数据结构，临时工)

    LogPLR *L;  
    int num_levels;
    int max_levels;

    int group_id;

    // method
} Group;

/**
 * @brief 学习索引段Frame包装结构 mapping table
 * 
 */

typedef struct Counter {
    int group_write_cnt;
    int group_read_cnt;
    int group_read_acc_hit;
    int group_read_noacc_hit;
    int group_read_noacc_miss;
    int group_double_read;
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
    uint64_t *o_maptbl;
    uint64_t o_maptbl_cnt;
    uint64_t o_maptbl_hit;

    uint64_t o_maptbl_maxLBA; 
    uint64_t o_maptbl_minLBA;


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
    Segment segments[MAX_CRB_SEGMENTS];  
    uint8_t num_segments;  // 256
} CRB;



void crb_init(CRB *crb);
void crb_insert_segment(CRB *crb, Segment* seg);
unsigned int crb_find_segment(CRB *crb, uint8_t lpa);

int Segment_is_valid(Segment *seg, uint32_t x);
uint32_t Segment_gety(Segment *seg, bool check, uint32_t x);
uint32_t Segment_gety_upper(Segment *seg, bool check, uint32_t x);
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
int32_t binary_search(LogPLR *logplr, Segment *seg);
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
// LRUCache *createLRUCache(int maxsize);
void FrameGroup_init(FrameGroup *framegroup, double gamma);
void FrameGroup_static(FrameGroup *framegroup);
void FrameGroup_update(FrameGroup *framegroup, Point* points, int num_points);
uint64_t FrameGroup_lookup(FrameGroup *framegroup, uint64_t lpn, Segment *seg);
uint64_t PPN2VPPN(uint64_t pgidx);


#endif