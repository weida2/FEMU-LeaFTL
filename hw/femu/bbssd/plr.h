#ifndef __FEMU_PLR_H
#define __FEMU_PLR_H

#include <math.h>
#include <stdlib.h>
#include "../nvme.h"
/* For PLR learn
 * 
*/
typedef struct Point {
    uint32_t x;        // lpn
    uint32_t y;        // ppn
} Point;


// no use
typedef struct Points {
    Point points[1024];
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
    unsigned char *bitmap_upper;
} Bitmap;


typedef struct Segment {
    double k;      // 2 Bytes  use Avg_slope ? for more accurate
    double b;      // 4 Bytes
    uint32_t   x1;       // start_index 1 Bytes(256)  ?
    uint32_t   x2;       // end_inedex  1 Bytes 
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
    int segment_id[1024];
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

#endif /* __FEMU_PLR_H */