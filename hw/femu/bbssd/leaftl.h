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



typedef struct CRB{
    Segment segments[MAX_CRB_SEGMENTS];  
    uint8_t num_segments;  // 256
} CRB;



#define  error_bound 0   // 误差范围
// #include "test.hpp"

#define INVALID_PPA     (~(0ULL))
#define INVALID_LPN     (~(0ULL))
#define UNMAPPED_PPA    (~(0ULL))

#define WB_Entries      1024 
// 1k * 4KB = 4MB 假设SSD内写入缓冲区的大小为 4MB 可以缓存1k条lpa：ppa条目

enum {
    NAND_READ =  0,
    NAND_WRITE = 1,
    NAND_ERASE = 2,

    NAND_READ_LATENCY = 40000,
    NAND_PROG_LATENCY = 200000,
    NAND_ERASE_LATENCY = 2000000,
};

enum {
    USER_IO = 0,
    GC_IO = 1,
};

enum {
    NORMAL_IO = 0,
    LEAFTL_IO = 1,
    DFTL_IO   = 2,
    LEAEDFTL_IO = 3,
    
};

enum {
    SEC_FREE = 0,
    SEC_INVALID = 1,
    SEC_VALID = 2,

    PG_FREE = 0,
    PG_INVALID = 1,
    PG_VALID = 2
};

enum {
    FEMU_ENABLE_GC_DELAY = 1,
    FEMU_DISABLE_GC_DELAY = 2,

    FEMU_ENABLE_DELAY_EMU = 3,
    FEMU_DISABLE_DELAY_EMU = 4,

    FEMU_RESET_ACCT = 5,
    FEMU_ENABLE_LOG = 6,
    FEMU_DISABLE_LOG = 7,

    FEMU_ENABLE_LEAFTL_IO  = 8,
    FEMU_ENABLE_DFTL_IO   = 9,
    FEMU_ENABLE_LEAEDFTL_IO = 12,
    FEMU_ENABLE_NOMAL_IO  = 13,
    FEMU_Group_Static     = 10,
    FEMU_DFTL_Static      = 11,


};


#define BLK_BITS    (16)
#define PG_BITS     (16)
#define SEC_BITS    (8)
#define PL_BITS     (8)
#define LUN_BITS    (8)
#define CH_BITS     (7)

/* describe a physical page addr */
struct ppa {
    union {
        struct {
            uint64_t blk : BLK_BITS;
            uint64_t pg  : PG_BITS;
            uint64_t sec : SEC_BITS;
            uint64_t pl  : PL_BITS;
            uint64_t lun : LUN_BITS;
            uint64_t ch  : CH_BITS;
            uint64_t rsv : 1;
        } g;

        uint64_t ppa;
    };
};

typedef int nand_sec_status_t;

struct nand_page {
    nand_sec_status_t *sec;
    int nsecs;
    int status;
};

struct nand_block {
    struct nand_page *pg;
    int npgs;
    int ipc; /* invalid page count */
    int vpc; /* valid page count */
    int erase_cnt;
    int wp; /* current write pointer */
};

struct nand_plane {
    struct nand_block *blk;
    int nblks;
};

struct nand_lun {
    struct nand_plane *pl;
    int npls;
    uint64_t next_lun_avail_time;
    bool busy;
    uint64_t gc_endtime;
};

struct ssd_channel {
    struct nand_lun *lun;
    int nluns;
    uint64_t next_ch_avail_time;
    bool busy;
    uint64_t gc_endtime;
};

struct ssdparams {
    int secsz;        /* sector size in bytes */
    int secs_per_pg;  /* # of sectors per page */
    int pgs_per_blk;  /* # of NAND pages per block */
    int blks_per_pl;  /* # of blocks per plane */
    int pls_per_lun;  /* # of planes per LUN (Die) */
    int luns_per_ch;  /* # of LUNs per channel */
    int nchs;         /* # of channels in the SSD */

    int pg_rd_lat;    /* NAND page read latency in nanoseconds */
    int pg_wr_lat;    /* NAND page program latency in nanoseconds */
    int blk_er_lat;   /* NAND block erase latency in nanoseconds */
    int ch_xfer_lat;  /* channel transfer latency for one page in nanoseconds
                       * this defines the channel bandwith
                       */

    double gc_thres_pcent;
    int gc_thres_lines;
    double gc_thres_pcent_high;
    int gc_thres_lines_high;
    bool enable_gc_delay;

    /* below are all calculated values */
    int secs_per_blk; /* # of sectors per block */
    int secs_per_pl;  /* # of sectors per plane */
    int secs_per_lun; /* # of sectors per LUN */
    int secs_per_ch;  /* # of sectors per channel */
    int tt_secs;      /* # of sectors in the SSD */

    int pgs_per_pl;   /* # of pages per plane */
    int pgs_per_lun;  /* # of pages per LUN (Die) */
    int pgs_per_ch;   /* # of pages per channel */
    int tt_pgs;       /* total # of pages in the SSD */

    int blks_per_lun; /* # of blocks per LUN */
    int blks_per_ch;  /* # of blocks per channel */
    int tt_blks;      /* total # of blocks in the SSD */

    int secs_per_line;
    int pgs_per_line;
    int blks_per_line;
    int tt_lines;

    int pls_per_line;
    int luns_per_line;

    int pls_per_ch;   /* # of planes per channel */
    int tt_pls;       /* total # of planes in the SSD */

    int tt_luns;      /* total # of LUNs in the SSD */
};

typedef struct line {
    int id;  /* line id, the same as corresponding block id */
    int ipc; /* invalid page count in this line */
    int vpc; /* valid page count in this line */
    QTAILQ_ENTRY(line) entry; /* in either {free,victim,full} list */
    /* position in the priority queue for victim lines */
    size_t                  pos;
} line;

/* wp: record next write addr */
struct write_pointer {
    struct line *curline;
    int ch;
    int lun;
    int pg;
    int blk;
    int pl;
};

struct line_mgmt {
    struct line *lines;
    /* free line list, we only need to maintain a list of blk numbers */
    QTAILQ_HEAD(free_line_list, line) free_line_list;
    pqueue_t *victim_line_pq;
    //QTAILQ_HEAD(victim_line_list, line) victim_line_list;
    QTAILQ_HEAD(full_line_list, line) full_line_list;
    int tt_lines;
    int free_line_cnt;
    int victim_line_cnt;
    int full_line_cnt;
};

struct nand_cmd {
    int type;
    int cmd;
    int64_t stime; /* Coperd: request arrival time */
};

typedef struct Write_Buffer {
    uint32_t LPA;
    uint32_t PPA;
} Write_Buffer;

struct ssd {
    char *ssdname;
    struct ssdparams sp;
    struct ssd_channel *ch;
    struct ppa *maptbl; /* page level mapping table */
    uint64_t *rmap;     /* reverse mapptbl, assume it's stored in OOB */ 
    struct write_pointer wp;
    struct line_mgmt lm;


    // leaftl.struct
    struct Write_Buffer WB[WB_Entries + 2];
    int    num_write_entries;
    int    hit_wb;
    FrameGroup l_maptbl;
    bool flush; 



    /* lockless ring for communication with NVMe IO thread */
    struct rte_ring **to_ftl;
    struct rte_ring **to_poller;
    bool *dataplane_started_ptr;
    QemuThread ftl_thread;
};

void ssd_init(FemuCtrl *n);

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


// void Group_compact(Group *group);
// void FrameGroup_do_gc(FrameGroup *framegroup);
// void FrameGroup_compact(FrameGroup *framegrop);

#endif