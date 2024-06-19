#include "标头.h"

static inline int tcqCheck(TC_QUEUE_STRUCT const* const tcq);
static inline double tpGetRealTargetVel(TP_STRUCT const* const tp, TC_STRUCT const* const tc);
static inline double tpGetRealFinalVel(TP_STRUCT const* const tp, TC_STRUCT const* const tc, TC_STRUCT const* const nexttc);
static int tpCheckEndCondition(TP_STRUCT const* const tp, TC_STRUCT* const tc, TC_STRUCT const* const nexttc);
static int tpCompleteSegment(TP_STRUCT* const tp, TC_STRUCT* const tc);
static void tpHandleEmptyQueue(TP_STRUCT* const tp);
static int tpSetCycleTime(TP_STRUCT* const tp, double secs);
static int tp_init();
static tp_err_t tpActiveSegment(TP_STRUCT* const tp, TC_STRUCT* const tc);
static double tpGetFeedScale(TP_STRUCT const* const tp, TC_STRUCT const* const tc);

static TP_STRUCT coord_tp;
static int mot_comp_id;
double trajCycleTime = 0.0001;
PmCartesian carte_pos_cmd = { 0 };
double currentvel=0.0;

double pmSqrt(double x) {
	if (x > 0.0) {
		return sqrt(x);
	}
	if (x > SQRT_FUZZ) {
		return 0.0;
	}
	return 0.0;
}

double saturate(double x, double max) {
	if (x > max) {
		return max;
	}
	else if (x < (-max)) {
		return -max;
	}
	else {
		return x;
	}
}

int tcInit(TC_STRUCT* const tc, double cycle_time) {
	tc->cycle_time = cycle_time;
	/*TODO：*/
	//tc->tolerance = 0.001;
	return 0;
}

int tcInitKinkProperties(TC_STRUCT* tc) {
	tc->kink_vel = -1.0;
	tc->kink_accel_reduce_prev = 0.0;
	tc->kink_accel_reduce = 0.0;
	return 0;
}

int tcSetupMotion(TC_STRUCT* const tc, double vel, double ini_maxvel, double acc) {
	tc->maxaccel = acc;
	tc->maxvel = ini_maxvel;
	tc->reqvel = vel;
	tcInitKinkProperties(tc);
	return 0;
}

int pmCartCartSub(PmCartesian const* const v1, PmCartesian const* const v2, PmCartesian* const vout) {
	vout->x = v1->x - v2->x;
	vout->y = v1->y - v2->y;
	return 1;
}

int pmCartMag(PmCartesian const* const v, double* d) {
	*d = pmSqrt(pmSq(v->x) + pmSq(v->y));
	return 1;
}

int pmCartUnitEq(PmCartesian* const v) {
	double size = pmSqrt(pmSq(v->x) + pmSq(v->y));
	v->x /= size;
	v->y /= size;
	return 1;
}

double pmLineTarget(PmCartLine* const line, PmCartesian const*const start, PmCartesian const* const end) {
	line->start = *start;
	line->end = *end;
	pmCartCartSub(end, start, &line->uVec);
	line->tmag = pmSqrt(pmSq(line->uVec.x) + pmSq(line->uVec.y));
	pmCartMag(&line->uVec, &line->tmag);
	pmCartUnitEq(&line->uVec);
	return line->tmag;
}

int PmCartScalMultEq(PmCartesian* const v, double d) {
	v->x *= d;
	v->y *= d;
	return 1;
}

int pmCartScalMult(PmCartesian const* const v1, double d, PmCartesian* const vout) {
	if (v1 != vout) {
		*vout = *v1;
	}
	return PmCartScalMultEq(vout, d);
}

int pmCartCartAdd(PmCartesian const* const v1, PmCartesian const* const v2, PmCartesian* const vout) {
	vout->x = v1->x + v2->x;
	vout->y = v1->y + v2->y;
	return 1;
}

int pmCartLinePoint(PmCartLine const* const line, double len, PmCartesian* const point) {
	pmCartScalMult(&line->uVec, len, point);
	pmCartCartAdd(&line->start, point, point);
	return 1;
}

int tcGetPosReal(TC_STRUCT const* const tc, int of_point, PmCartesian* const pos) {
	PmCartesian xyz;
	double progress = 0.0;
	switch (of_point) {
	case TC_GET_PROGRESS:
		progress = tc->progress;
		break;
	case TC_GET_ENDPOINT:
		progress = tc->target;
		break;
	case TC_GET_STARTPOINT:
		progress = 0.0;
		break;
	}
	double angle = 0.0;		
	pmCartLinePoint(&tc->coords.line, progress * tc->coords.line.tmag / tc->target, &xyz);
	pos->x = xyz.x;
	pos->y = xyz.y;
	return 0;
}

int tcGetEndPoint(TC_STRUCT const* const tc, PmCartesian* const out) {
	tcGetPosReal(tc, TC_GET_ENDPOINT, out);
	return 0;
}

int tcqPut(TC_QUEUE_STRUCT* const tcq, TC_STRUCT const* const tc) {
	tcq->queue[tcq->end] = *tc;
	tcq->_len++;
	tcq->end = (tcq->end + 1) % tcq->size;
	return 0;
}

inline int tpAddSegmentToQueue(TP_STRUCT* const tp, TC_STRUCT* const tc, int inc_id) {
	tc->id = tp->nextId;
	tcqPut(&tp->queue, tc);
	if (inc_id) {
		tp->nextId++;
	}
	tcGetEndPoint(tc, &tp->goalPos);
	tp->done;
	return 0;
}

int tpAddLine(TP_STRUCT* const tp,PmCartesian end,double vel,double ini_maxvel,double acc) {
	//初始化一个段结构体
	TC_STRUCT tc = { 0 };
	tcInit(&tc, tp->cycletime);
	tcSetupMotion(&tc, vel, ini_maxvel, acc);
	tc.target = pmLineTarget(&tc.coords.line, &tp->goalPos, &end);
	tpAddSegmentToQueue(tp, &tc, TRUE);
	return 0;
}

int pmCartCartDot(PmCartesian const* const v1, PmCartesian const* const v2, double* d) {
	*d = v1->x * v2->x + v1->y * v2->y;
	return 0;
}

static inline int tcqCheck(TC_QUEUE_STRUCT const* const tcq) {
	if ((0 == tcq) || (0 == tcq->queue)) {
		return -1;
	}
	return 0;
}

TC_STRUCT* tcqItem(TC_QUEUE_STRUCT const* const tcq, int n) {
	if (tcqCheck(tcq) || (n < 0) || (n >= tcq->_len)) return NULL;
	return &(tcq->queue[(tcq->start + n) % tcq->size]);
}

int tcGetPos(TC_STRUCT const* const tc, PmCartesian* const out) {
	tcGetPosReal(tc, TC_GET_PROGRESS, out);
	return 0;
}

void tpCalculateTrapezoidalAccel(TP_STRUCT const* const tp, TC_STRUCT* const tc, TC_STRUCT const * const nexttc,
	double* const acc, double* const vel_desired) {
	double tc_target_vel = tc->reqvel;
	double tc_finalvel = tc->finalvel;
	double dx = tc->target - tc->progress;
	double maxaccel = tc->maxaccel;
	double discr_term1 = pmSq(tc_finalvel);
	double discr_term2 = maxaccel * (2.0 * dx - tc->currentvel * tc->cycle_time);
	double tmp_adt = maxaccel * tc->cycle_time * 0.5;
	double discr_term3 = pmSq(tmp_adt);
	double discr = discr_term1 + discr_term2 + discr_term3;
	double maxnewvel = -tmp_adt;
	if (discr > discr_term3) {
		maxnewvel += pmSqrt(discr);
	}
	double newvel = saturate(maxnewvel, tc_target_vel);
	double dt = fmax(tc->cycle_time, TP_TIME_EPSILON);
	double maxnewaccel = (newvel - tc->currentvel) / dt;
	*acc = saturate(maxnewaccel, maxaccel);
	*vel_desired = maxnewvel;
}

int tcUpdateDistFromAccel(TC_STRUCT* const tc, double acc, double vel_desired) {
	double v_next = tc->currentvel + acc * tc->cycle_time;
	if (v_next < 0.0) {
		v_next = 0.0;
		if ((tc->target - tc->progress) < (tc->cycle_time * tc->currentvel)) {
			tc->progress = tc->target;
		}
	}
	else {
		double displacement = (v_next + tc->currentvel) * tc->cycle_time * 0.5;
		tc->progress += displacement;
	}
	tc->currentvel = v_next;
	tc->on_final_decel = (fabs(vel_desired - tc->currentvel) < TP_VEL_EPSILON) && (acc < 0.0);
	return 0;
}

int emcPoseAdd(PmCartesian const* const p1, PmCartesian const* const p2, PmCartesian* const out) {
	out->x = p1->x + p2->x;
	out->y = p1->y + p2->y;
	return 1;
}

int emcPoseSelfAdd(PmCartesian* const self, PmCartesian const* const p2) {
	emcPoseAdd(self, p2, self);
	return 1;
}

int tpAddCurrentPos(TP_STRUCT* const tp, PmCartesian const* const disp) {
	emcPoseSelfAdd(&tp->currentPos, disp);
	return 1;
}

static inline double tpGetRealTargetVel(TP_STRUCT const* const tp, TC_STRUCT const* const tc) {
	if (!tc) {
		return 0.0;
	}
	return tc->reqvel;
}

static inline double tpGetRealFinalVel(TP_STRUCT const* const tp, TC_STRUCT const* const tc, TC_STRUCT const* const nexttc) {
	double v_target_this = tpGetRealTargetVel(tp, tc);
	double v_target_next = 0.0;
	if (nexttc) {
		v_target_next = tpGetRealTargetVel(tp, nexttc);
	}
	double v_target = fmin(v_target_this, v_target_next);
	return fmin(tc->finalvel, v_target);
}

int sat_inplace(double* const x, double max) {
	if (*x > max) {
		*x = max;
		return 1;
	}
	else if (*x < -max) {
		*x = -max;
		return -1;
	}
	return 0;
}

static int tpCheckEndCondition(TP_STRUCT const* const tp, TC_STRUCT* const tc, TC_STRUCT const* const nexttc) {
	tc->cycle_time = tp->cycletime;
	double dx = tc->target - tc->progress;
	if (dx <= TP_POS_EPSILON) {
		tc->progress = tc->target;
		tc->remove = 1;
	}
	double v_f = tpGetRealFinalVel(tp, tc, nexttc);
	double v_avg = (tc->currentvel + v_f) / 2.0;
	double dt = TP_TIME_EPSILON / 2.0;
	if (v_avg > TP_VEL_EPSILON) {
		dt = fmax(dt, dx / v_avg);
	}
	else {
		if (dx > (v_avg * tp->cycletime) && dx > TP_POS_EPSILON) {
			return 1;
		}
	}
	double dv = v_f - tc->currentvel;
	double a_f = dv / dt;
	double a_max = tc->maxaccel;
	double a = a_f;
	int recalc = sat_inplace(&a, a_max);
	if (recalc) {
		double disc = pmSq(tc->currentvel / a) + 2.0 / a * dx;
		if (disc < 0) {
			return 1;
		}
		if (disc < TP_TIME_EPSILON * TP_TIME_EPSILON) {
			dt = -tc->currentvel / a;
		}
		else if (a > 0) {
			dt = -tc->currentvel / a + pmSqrt(disc);
		}
		else {
			dt = -tc->currentvel / a - pmSqrt(disc);
		}
		v_f = tc->currentvel + dt * a;
	}
	if (dt < TP_TIME_EPSILON) {
		tc->progress = tc->target;
		tc->remove = 1;
	}
	return 1;
}

int	emcPoseSub(PmCartesian const* const p1, PmCartesian const* const p2, PmCartesian* const out) {
	out->x = p1->x - p2->x;
	out->y = p1->y - p2->y;
	return 1;
}

int emcPoseSelfSub(PmCartesian* const self, PmCartesian const* const p2) {
	return emcPoseSub(self, p2, self);
}

int tpUpdateCycle(TP_STRUCT* const tp, TC_STRUCT* const tc, TC_STRUCT const* const nexttc) {
	PmCartesian before;
	tcGetPos(tc, &before);
	if (!tc->blending_next) {
		tc->vel_at_blend_start = tc->currentvel;
	}
	int res_accel = 1;
	double acc = 0, vel_desired = 0;
	tpCalculateTrapezoidalAccel(tp, tc, nexttc, &acc, &vel_desired);
	tcUpdateDistFromAccel(tc, acc, vel_desired);
	tp->currentVel = tc->currentvel;
	tpCheckEndCondition(tp, tc, nexttc);
	PmCartesian displacement;
	tcGetPos(tc, &displacement);
	emcPoseSelfSub(&displacement, &before);
	tpAddCurrentPos(tp, &displacement);
	return 0;
}

double tpCalculateTriangleVel(TC_STRUCT const* tc) {
	double acc_scaled = tc->maxaccel;
	double length = tc->target;
	return pmSqrt(acc_scaled * length);
}

int findIntersectionAngel(PmCartesian const* const u1, PmCartesian const* const u2, double* const theta) {
	double dot;
	pmCartCartDot(u1, u2, &dot);	//计算向量的内积
	if (dot > 1.0 || dot < -1.0) {
		sat_inplace(&dot, 1.0);
	}
	*theta = acos(-dot) / 2.0;
	return 0;
}

int tcGetEndAccelUnitVector(TC_STRUCT const* const tc, PmCartesian* const out) {
	*out = tc->coords.line.uVec;
	return 1;
}

int tcGetStartAccelUnitVector(TC_STRUCT const* const tc, PmCartesian* const out) {
	*out = tc->coords.line.uVec;
	return 1;
}

int tpComputeBlendVelocity(TC_STRUCT const* tc, TC_STRUCT const* nexttc,
	double target_vel_this, double target_vel_next,
	double* v_blend_this, double* v_blend_next, double* v_blend_net) {
	if (!nexttc || !tc || !v_blend_this || !v_blend_next) {
		return 0;
	}
	double acc_this = tc->maxaccel;
	double acc_next = nexttc->maxaccel;
	double v_reachable_this = fmin(tpCalculateTriangleVel(tc), target_vel_this);
	double v_reachable_next = fmin(tpCalculateTriangleVel(nexttc), target_vel_next);
	double t_max_this = tc->target / v_reachable_this;
	double t_max_next = nexttc->target / v_reachable_next;
	double t_max_reachable = fmin(t_max_this, t_max_next);
	double t_min_blend_this = v_reachable_this / acc_this;
	double t_min_blend_next = v_reachable_next / acc_next;
	double t_max_blend = fmax(t_min_blend_this,t_min_blend_next);
	double t_blend = fmin(t_max_reachable, t_max_blend);
	*v_blend_this = fmin(v_reachable_this, t_blend * acc_this);
	*v_blend_next = fmin(v_reachable_next, t_blend * acc_next);
	double theta;
	PmCartesian v1, v2;
	tcGetEndAccelUnitVector(tc, &v1);
	tcGetStartAccelUnitVector(nexttc, &v2);
	findIntersectionAngel(&v1, &v2, &theta);
	double cos_theta = cos(theta);

	if (tc->tolerance > 0) {
		double tblend_vel;
		const double min_cos_theta = cos(PM_PI / 2.0 - TP_MIN_ARC_ANGLE);
		if (cos_theta > min_cos_theta) {
			tblend_vel = 2.0 * (acc_this * tc->tolerance / cos_theta);
			*v_blend_this = fmin(*v_blend_this, tblend_vel);
			*v_blend_next = fmin(*v_blend_next, tblend_vel);
		}
	}
	if (v_blend_net) {
		*v_blend_net = sin(theta) * (*v_blend_this + *v_blend_next) / 2.0;
	}
	return 1;
}

static double tpGetFeedScale(TP_STRUCT const* const tp, TC_STRUCT const* const tc) {
	if (!tc) {
		return 0.0;
	}
	bool pausing = tp->pausing;
	if (pausing) {
		return 0.0;
	}
	else if (tc->is_blending) {
		return fmin(0.1, 1.0);
	}
	return 0.1;
}

static void tpUpdateBlend(TP_STRUCT* const tp, TC_STRUCT* const tc, TC_STRUCT* const nexttc) {
	if (!nexttc) {
		return;
	}
	double save_vel = nexttc->target_vel;
	if (tpGetFeedScale(tp,nexttc) > TP_VEL_EPSILON) {
		double dv = tc->vel_at_blend_start - tc->currentvel;
		double vel_start = fmax(tc->vel_at_blend_start, TP_VEL_EPSILON);
		double blend_progress = fmax(fmin(dv / vel_start, 1.0), 0.0);
		double blend_scale = tc->vel_at_blend_start / tc->blend_vel;
		nexttc->target_vel = blend_progress * nexttc->blend_vel * blend_scale;
		nexttc->is_blending = true;
	}
	else {
		nexttc->target_vel = 0.0;
	}
	tpUpdateCycle(tp, nexttc, NULL);
	nexttc->target_vel = save_vel;
}

int tpUpdateMovementStatus(TP_STRUCT* const tp, TC_STRUCT* const tc) {
	PmCartesian tc_pos;
	tcGetEndPoint(tc, &tc_pos);
	tp->execId = tc->id;
	currentvel = tc->currentvel;
	return 0;
}

int tpDoParabolicBlending(TP_STRUCT* const tp, TC_STRUCT* const tc, TC_STRUCT* const nexttc) {
	tpUpdateBlend(tp, tc, nexttc);
	if (tc->currentvel > nexttc->currentvel) {
		tpUpdateMovementStatus(tp, tc);
	}
	else {
		tpUpdateMovementStatus(tp, nexttc);
	}
	return TP_ERR_OK;
}

int tcIsBlending(TC_STRUCT* const tc) {
	int is_blending_next = tc->on_final_decel && (tc->currentvel < tc->blend_vel);
	tc->blending_next |= is_blending_next;
	return tc->blending_next;
}

int tpHandleRegularCycle(TP_STRUCT* const tp, TC_STRUCT* const tc, TC_STRUCT* const nexttc) {
	if (tc->remove) {
		return -1;
	}
	tc->cycle_time = tp->cycletime;
	tpUpdateCycle(tp, tc, nexttc);
	/* 抛物线速度混合 */
	double v_this=0.0, v_next = 0.0;
	double target_vel_this = tpGetRealTargetVel(tp,tc);
	double target_vel_next = tpGetRealTargetVel(tp, nexttc);
	tpComputeBlendVelocity(tc, nexttc, target_vel_this, target_vel_next, &v_this, &v_next, NULL);
	tc->blend_vel = v_this;
	if (nexttc) {
		nexttc->blend_vel = v_next;
	}
	if (nexttc && tcIsBlending(tc)) {
		tpDoParabolicBlending(tp, tc, nexttc);
	}
	else {
		tpUpdateMovementStatus(tp, tc);
	}
	return 0;
}

int tcqPop(TC_QUEUE_STRUCT* const tcq) {
	if (tcqCheck(tcq)) {
		return -1;
	}
	tcq->start = (tcq->start + 1) % tcq->size;
	tcq->_len--;
	return 1;
}

static int tpCompleteSegment(TP_STRUCT* const tp, TC_STRUCT* const tc) {
	tc->remove = 0;
	tc->cycle_time = tp->cycletime;
	tc->currentvel = 0.0;
	tc->is_blending = 0;
	tcqPop(&tp->queue);
	return 1;
}

int tcqInit(TC_QUEUE_STRUCT* const tcq) {
	tcq->_len = 0;
	tcq->start = tcq->end = 0;
	return 0;
}

static void tpHandleEmptyQueue(TP_STRUCT* const tp) {
	tcqInit(&tp->queue);
	tp->goalPos = tp->currentPos;
	tp->done = 1;
	tp->execId = 0;
}

int tpRunCycle(TP_STRUCT* const tp, long period) {
	TC_STRUCT* tc;
	TC_STRUCT* nexttc;
	tc = tcqItem(&tp->queue, 0);
	nexttc = tcqItem(&tp->queue, 1);
	if (!tc) {
		tpHandleEmptyQueue(tp);
		return 0;
	}
	int res_activate = tpActiveSegment(tp, tc);
	tpHandleRegularCycle(tp, tc, nexttc);
	if (tc->remove) {
		tpCompleteSegment(tp, tc);
	}
	return 0;
}

int tcqCreate(TC_QUEUE_STRUCT* const tcq, int _size, TC_STRUCT* const tcSpace) {
	tcq->queue = tcSpace;
	tcq->size = _size;
	tcqInit(tcq);
	return 0;
}

TC_STRUCT queueTcSpace[DEFAULT_TC_QUEUE_SIZE + 10];

int tpClear(TP_STRUCT* const tp) {
	tcqInit(&tp->queue);
	tp->queueSize = 0;
	tp->goalPos = tp->currentPos;
	tp->nextId = 0;
	tp->execId = 0;
	return 0;
}

int tpInit(TP_STRUCT* const tp) {
	tp->cycletime = 0.0;
	ZERO_EMC_POSE(tp->currentPos);
	return tpClear(tp);
}

int tpCreat(TP_STRUCT* const tp, int _queueSize, int id) {
	if (0 == tp) {
		return 0;
	}
	if (_queueSize <= 0) {
		tp->queueSize = DEFAULT_TC_QUEUE_SIZE;
	}
	else {
		tp->queueSize = _queueSize;
	}
	TC_STRUCT* const tcSpace = queueTcSpace;
	tcqCreate(&tp->queue, tp->queueSize, tcSpace);
	return tpInit(tp);
}

static int tpSetCycleTime(TP_STRUCT* const tp, double secs) {
	tp->cycletime = secs;
	return 0;
}

int tpSetPos(TP_STRUCT* const tp, PmCartesian const* const pos) {
	tp->currentPos = *pos;
	tp->goalPos = *pos;
	return 0;
}

static int tp_init() {
	tpCreat(&coord_tp, DEFAULT_TC_QUEUE_SIZE, mot_comp_id);
	tpSetCycleTime(&coord_tp, trajCycleTime);
	tpSetPos(&coord_tp, &carte_pos_cmd);
	return 0;
}

static tp_err_t tpActiveSegment(TP_STRUCT* const tp, TC_STRUCT* const tc) {
	tc->on_final_decel = 0;
	tc->blending_next = 0;
	return TP_ERR_OK;
}

int main() {
	std::vector<double> record_speed;
	std::vector<double> record_progress;
	std::vector<double> record_acc;
	tp_init();
	PmCartesian pos = { 10,10 };
	PmCartesian pos_next = { 20,20 };
	tpAddLine(&coord_tp, pos, 250, 20, 10000);
	tpAddLine(&coord_tp, pos_next, 250, 20, 10000);
	int times = 0;
	while (1) {
		times++;
		tpRunCycle(&coord_tp, 0.0001);
		record_speed.push_back(currentvel);
		record_progress.push_back(coord_tp.currentPos.x);
		if (coord_tp.done==1) {
			break;
		}
	}
	initgraph(1300, 500);
	for (int i = 0; i < record_progress.size(); ++i) {
		putpixel(i*0.8, record_progress[i]*10, 0x0000AA);
		putpixel(i * 0.8, record_speed[i], 0xFF0000);
	}
	_getch();
	closegraph();
	return 1;
}
