#pragma once
#include<cmath>
#include<graphics.h>
#include<conio.h>
#include<vector>

#define pmSq(x) ((x)*(x))
#define SQRT_FUZZ (-1.0e-6)
#define TP_TIME_EPSILON 1e-12
#define TC_GET_PROGRESS 0
#define TC_GET_STARTPOINT 1
#define TC_GET_ENDPOINT 2
#define TP_MIN_ARC_ANGLE 1e-3
#define PM_PI      3.14159265358979323846
#define DEFAULT_TC_QUEUE_SIZE 2000
#define TP_POS_EPSILON   1e-12
#define TP_TIME_EPSILON  1e-12
#define TP_VEL_EPSILON   1e-8
#define TP_ACCEL_EPSILON 1e-4
#define ZERO_EMC_POSE(pos) do{\
pos.x=0.0;\
pos.y = 0.0;  }while (0)

typedef struct {
	double x, y;
} PmCartesian;

typedef struct {
	PmCartesian start;
	PmCartesian end;
	double tmag;
	PmCartesian uVec;
} PmCartLine;

typedef struct {
	double cycle_time;
	// 位置
	double target;
	double progress;
	//速度
	double reqvel;			//F命令要求的速度
	double currentvel;		//跟踪当前速度
	double finalvel;		//段末端需要的速度
	double target_vel;
	double maxvel;
	double kink_vel;
	double kink_accel_reduce_prev;
	double kink_accel_reduce;
	double blend_vel;
	//加速度
	double maxaccel;		//按任务计算加速度
	union {
		PmCartLine line;
	} coords;
	int id;
	double tolerance;
	int remove;
	int on_final_decel;
	int blending_next;
	int vel_at_blend_start;
	int is_blending;
} TC_STRUCT;

typedef struct {
	TC_STRUCT* queue;
	int size;
	int _len;
	int start, end;
} TC_QUEUE_STRUCT;

typedef struct {
	TC_QUEUE_STRUCT queue;
	int queueSize;
	double cycletime;
	PmCartesian currentPos;
	PmCartesian goalPos;
	double currentVel;
	int nextId;
	int execId;
	int done;
	int pausing;
} TP_STRUCT;

typedef enum {
	TP_ERR_INVALID = -9,
	TP_ERR_INPUT_TYPE = -8,
	TP_ERR_TOLERANCE = -7,
	TP_ERR_RADIUS_TOO_SMALL = -6,
	TP_ERR_GEOM = -5,
	TP_ERR_RANGE = -4,
	TP_ERR_MISSING_OUTPUT = -3,
	TP_ERR_MISSING_INPUT = -2,
	TP_ERR_FAIL = -1,
	TP_ERR_OK = 0,
	TP_ERR_NO_ACTION,
	TP_ERR_SLOWING,
	TP_ERR_STOPPED,
	TP_ERR_WAITING,
	TP_ERR_ZERO_LENGTH,
	TP_ERR_REVERSE_EMPTY,
	TP_ERR_LAST
} tp_err_t;

typedef struct emcmot_status_t {
	double current_vel;
} emcmot_status_t;

double pmSqrt(double x);
double saturate(double x, double max);
int tcInit(TC_STRUCT* const tc, double cycle_time);
int tcInitKinkProperties(TC_STRUCT* tc);
int tcSetupMotion(TC_STRUCT* const tc, double vel, double ini_maxvel, double acc);
int pmCartCartSub(PmCartesian const* const v1, PmCartesian const* const v2, PmCartesian* const vout);
int pmCartMag(PmCartesian const* const v, double* d);
int pmCartUnitEq(PmCartesian* const v);
double pmLineTarget(PmCartLine* const line, PmCartesian const* const start, PmCartesian const* const end);
int PmCartScalMultEq(PmCartesian* const v, double d);
int pmCartScalMult(PmCartesian const* const v1, double d, PmCartesian* const vout);
int pmCartCartAdd(PmCartesian const* const v1, PmCartesian const* const v2, PmCartesian* const vout);
int pmCartLinePoint(PmCartLine const* const line, double len, PmCartesian* const point);
int tcGetPosReal(TC_STRUCT const* const tc, int of_point, PmCartesian* const pos);
int tcGetEndPoint(TC_STRUCT const* const tc, PmCartesian* const out);
int tcqPut(TC_QUEUE_STRUCT* const tcq, TC_STRUCT const* const tc);
inline int tpAddSegmentToQueue(TP_STRUCT* const tp, TC_STRUCT* const tc, int inc_id);
int tpAddLine(TP_STRUCT* const tp, PmCartesian end, double vel, double ini_maxvel, double acc);
int pmCartCartDot(PmCartesian const* const v1, PmCartesian const* const v2, double* d);
TC_STRUCT* tcqItem(TC_QUEUE_STRUCT const* const tcq, int n);
int tcGetPos(TC_STRUCT const* const tc, PmCartesian* const out);
void tpCalculateTrapezoidalAccel(TP_STRUCT const* const tp, TC_STRUCT* const tc, TC_STRUCT const* const nexttc,
	double* const acc, double* const vel_desired);
int tcUpdateDistFromAccel(TC_STRUCT* const tc, double acc, double vel_desired);
int emcPoseAdd(PmCartesian const* const p1, PmCartesian const* const p2, PmCartesian* const out);
int emcPoseSelfAdd(PmCartesian* const self, PmCartesian const* const p2);
int tpAddCurrentPos(TP_STRUCT* const tp, PmCartesian const* const disp);
int sat_inplace(double* const x, double max);
int tpUpdateCycle(TP_STRUCT* const tp, TC_STRUCT* const tc, TC_STRUCT const* const nexttc);
double tpCalculateTriangleVel(TC_STRUCT const* tc);
int findIntersectionAngel(PmCartesian const* const u1, PmCartesian const* const u2, double* const theta);
int tpComputeBlendVelocity(TC_STRUCT const* tc, TC_STRUCT const* nexttc,
	double target_vel_this, double target_vel_next,
	double* v_blend_this, double* v_blend_next, double* v_blend_net);
void tpUpdateBlend(TP_STRUCT* const tp, TC_STRUCT* const tc, TC_STRUCT* const nexttc);
int tpUpdateMovementStatus(TP_STRUCT* const tp, TC_STRUCT* const tc);
int tpDoParabolicBlending(TP_STRUCT* const tp, TC_STRUCT* const tc, TC_STRUCT* const nexttc);
int tpHandleRegularCycle(TP_STRUCT* const tp, TC_STRUCT* const tc, TC_STRUCT* const nexttc);
int tcqPop(TC_QUEUE_STRUCT* const tcq);
int tcqInit(TC_QUEUE_STRUCT* const tcq);
int tpRunCycle(TP_STRUCT* const tp, long period);
int tcqCreate(TC_QUEUE_STRUCT* const tcq, int _size, TC_STRUCT* const tcSpace);
int tpClear(TP_STRUCT* const tp);
int tpInit(TP_STRUCT* const tp);
int tpCreat(TP_STRUCT* const tp, int _queueSize, int id);
int tpSetPos(TP_STRUCT* const tp, PmCartesian const* const pos);
int emcPoseSelfSub(PmCartesian* const self, PmCartesian const* const p2);
int	emcPoseSub(PmCartesian const* const p1, PmCartesian const* const p2, PmCartesian* const out);
int tcGetEndAccelUnitVector(TC_STRUCT const* const tc, PmCartesian* const out);
int tcGetStartAccelUnitVector(TC_STRUCT const* const tc, PmCartesian* const out);
