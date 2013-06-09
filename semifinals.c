/*
 *
 * Written by:
 * Flavio 'darkjoker' Giobergia
 * Marco 'The Nano' Pomponio
 *
 */

#include <math.h>
#include <string.h>
#include "ZR_API.h"
#include "ZRGame.h"
#include "math_matrix.h"
#include "spheres_constants.h"

#ifdef ZRSIMULATION
extern void _init(void);
void *_start = &_init;
#endif

static int st; //DECL::VAR::st
static char ok; //DECL::VAR::ok
static char sp; //DECL::VAR::sp
static float a; //DECL::VAR::a
static int ast; //DECL::VAR::ast
static int act; //DECL::VAR::act
static int m; //DECL::VAR::m
static void mathVecMul (float *v, float *a, float n); //DECL::PROC::mathVecMul
static void funz (float a[3], float q[4], float v[3]); //DECL::PROC::funz
static void SetPosFaster (float c[3], float state[12]); //DECL::PROC::SetPosFaster
static int getInTime (float cother[3], float other[12], float c[3], float state[12]); //DECL::PROC::getInTime
static void spin (float c[3], float mystate[12]); //DECL::PROC::spin
static float timeStation (float c[3], float state[12]); //DECL::PROC::timeStation
static void revolve (float c[3], float mystate[12]); //DECL::PROC::revolve
static float checkTarget (float other[12], float c[3]); //DECL::PROC::checkTarget
static void meltIce (float state[12], float time); //DECL::PROC::meltIce

void ZRUser01(float *myState, float *otherState, float time)
{
//BEGIN::PROC::ZRUser


// These are just some costants, don't pay to much attention to them :)

#define	 TIME_OPP_STAT   140	// Will only start considering about going
				// to a station after TIME_OPP_STAT seconds.

#define	 ACCURACY_STAT   1.97f  // Number between 1 and 2, the higher it is,
				// the more accurate in deciding whether the 
				// opponent is headed towards the station.

#define	 FUEL_STAT	   25	// The minimum amount of fuel required to go
				// to a mining station (below this value, 
				// the sphere's not going anywhere


// msg  : contains the messages received from the opponent
// out  : contains the message we're sending our opponent
unsigned short msg,out=0;
// avoid: is true or false whether the collision avoidance system is activated or not
unsigned char avoid;
// Don't worry about x and v ^^
float x[2];
float v[3];
// stat : these are some coordinates of stuff, like stations, asteroids and
//		a couple of corners of the playground
float stat[6][3] =  {   {-0.5f,+0.31f,-0.55f},
			{+0.5f,-0.31f,+0.55f},
			{0,-0.35f,-0.2f},
			{0,+0.35f,+0.2f},
			{-0.5f,+0.65f,-0.55f},
			{+0.5f,-0.65f,+0.55f}};
msg=PgetMessage();
avoid = PisAvoidingCollision();

/*
 * y0b0tics protocol will be used only when 
 * we're controlling the SPH # 1
 */
 
if (!time) {
	if (myState[0]>0) { // SPH1
		a=0.4; //  a is the X coordinate of the nearest laser
		out=0x0400; // bit 11
		sp=1;
	}
	else
		a=-0.4;
}

/*
 * If the opponent is willing to work
 * on Indigens, we're going to Indigens
 * as well, but we'll try to revolve.
 * Otherwise we keep sending our
 * opponent a "I'll revolve on Opulens"
 * message.
 */
if (time<60 && time) {
	if ((msg>3 && msg<6) || ast==1) {
		out=5;
		m=0;
		act=2; // Revolve
		ast=1; // Indigens
	}
	else {
		out=3;
	}
}
if (time>45 && (msg==3 || msg==5) && PisRevolving(otherState) && act==2) {
	act=3;
	out=msg-1;
	m=(out==2);
}

/*
 * If we're headed to a mining station,
 * we're telling our opponents which one
 * it is
 */
else if (act==4 && st<2)
	out=7-st;
PsendMessage(out);

/*
 * First of all, the Laser is taken
 * (btw, there's a "minor" bug in this
 * condition (&& time<60) is needed
 * (will add that as soon as some space
 * is available)
 */
if (!PhaveLaser()) {
	/*
	 * If the opponent is trying to "steal"
	 * our laser, we will get the other one
	 */
	if (avoid && myState[0]*a>0) 
		a*=-1;
	/* v contains the
	 * laser coordinates
	 */
	v[0]=a;
	v[1]=0;
	v[2]=0;
	ZRSetPositionTarget(v);
}
else {
	if (time>TIME_OPP_STAT && st<2) {
		/*
		 * x[0] is given the value of the distance
		 * between us and our opponent
		 */
		mathVecSubtract (v,myState,otherState,3);
		x[0]=mathVecMagnitude (v,3);
		/*
		 * if our opponent is going to a station
		 * we decide to go to the other one, but
		 * only if we're able to arrive first
		 */
		if (checkTarget(otherState,stat[0])>ACCURACY_STAT) {
			/*
			 * Even if we're headed towards a station,
			 * if our opponent is going to
			 * the same station as us, we stop going
			 * there and, instead, we move to the
			 * closest corner of the playground,
			 * so that, even if we don't gain anything,
			 * at least we won't activate the collision
			 * avoidance system
			 */
			if (!st && x[0]<0.65f) {
					st=4;
			}
			else if (st==-1) {
				/*
				 * We go to the other station only if
				 * we can reach it first
				 */
				if (getInTime(stat[0],otherState,stat[1],myState)) {
					st=1;
					act=4;
				}
				else
					st=-2;
			}
		}
		/*
		 * Same stuff as before, but now
		 * the other station is checked.
		 */
		else if (checkTarget(otherState,stat[1])>ACCURACY_STAT) {
			if (st==1 && x[0]<0.65f) {
					st=5;
			}
			else if (st==-1) {
				if (getInTime(stat[1],otherState,stat[0],myState)) {
					st=0;
					act=4;
				}
				else
					st=-2;
			}
		}
		/*
		 * although our opponent is not going to 
		 * a station, we leave the asteroid we're
		 * working on to go to a station.
		 * The time we're leaving the asteroid is
		 * not fixed, but is calculated every time
		 * in order to leave the asteroid as late as
		 * possible
		 */
		else if (PgetPercentFuelRemaining()>FUEL_STAT && st==-1) {
			/*
			 * x[0] contains the value of time
			 * at which we have to leave the asteroid
			 * to reach a station
			 * 180-(time needed to reach a station)
			 */
			x[0]=timeStation(stat[0],myState);
			x[1]=timeStation(stat[1],myState);
				/*
				 * It is then decided which station is
				 * better, keeping in account that, if
				 * we're spinning, our opponent could be
				 * on our way to the station.
				 */
				if ((time>=x[0]) && (x[0]>x[1]) && !(act==3 && otherState[0]<0 && otherState[1]>0 && otherState[2]<0.2)) {
					act=4;
					st=0;
				}
				if ((x[0]<x[1]) && (time>=x[1]) && !(act==3 && otherState[0]>0 && otherState[1]<0 && otherState[2]>-0.2)) {
					act=4;
					st=1;
				}
		}
	}
	/*
	 * act (action) contains a value
	 * related to the action we're
	 * doing.
	 */
	switch (act) {
		/* act = 2 -> revolve */
		case 2:
			if (m) {
				meltIce (myState,time);
				if (PiceMelted()) 
					m=0;
			}
			else 
				revolve (stat[ast],myState);
			/*
			 * if the opponent is revolving as
			 * well, we keep revolving as well
			 * for WAIT_TIME seconds. After that
			 * seconds, if the opponent is 
			 * still revolving, we start spinning
			 */
			break;
		/* act = 3 -> spinning */
		case 3:
			if (m) {
				meltIce (myState,time);
				if (PiceMelted()) 
					m=0;
			}
			else 
				spin (stat[ast],myState);

			/*
			 * if the avoiding collision system
			 * is working, it probably means that
			 * our opponent wants to spin and he's
			 * close to us so, after a few seconds,
			 * we let him spin and we start revolving ^^
			 */
			break;
		/* act = 4 -> going somewhere */
		case 4:
			SetPosFaster(stat[st],myState);
			break;
		default:
			break;
	}
}
//END::PROC::ZRUser
}
void ZRInit01()
{
//BEGIN::PROC::ZRInit
st = -1;
ok = 0;
sp = -1;
a = 0.0f;
ast = 2;
act = 2;
m = 1;
//END::PROC::ZRInit
}
//User-defined procedures
static void mathVecMul (float *v, float *a, float n)
{
//BEGIN::PROC::mathVecMul
/*
 * The code is self-explanatory
 */
int i;
for (i=0;i<3;i++)
	v[i]=a[i]*n;
//END::PROC::mathVecMul
}
static void funz (float a[3], float q[4], float v[3])
{
//BEGIN::PROC::funz
/*
 * Magic happens here :)
 */
a[0] =  q[0] * v[3] + q[1] * v[2] - q[2] * v[1] + q[3] * v[0];
a[1] = -q[0] * v[2] + q[1] * v[3] + q[2] * v[0] + q[3] * v[1];
a[2] =  q[0] * v[1] - q[1] * v[0] + q[2] * v[3] + q[3] * v[2];
a[3] = -q[0] * v[0] - q[1] * v[1] - q[2] * v[2] + q[3] * v[3];

//END::PROC::funz
}
static void SetPosFaster (float c[3], float state[12])
{
//BEGIN::PROC::SetPosFaster
/*
 * SetPosFaster, as the name says,
 * is used to reach a point in
 * a shorter amount of seconds
 * than ZRSetPositionTarget ().
 */
float v[3];

mathVecSubtract(v,c,state,3);
mathVecMul (v,v,1.315f);
mathVecAdd (v,v,state,3);
ZRSetPositionTarget(v);
//END::PROC::SetPosFaster
}
static int getInTime (float cother[3], float other[12], float c[3], float state[12])
{
//BEGIN::PROC::getInTime
/*
 * getInTime () returns a boolean
 * value, which will be false if our
 * opponent is able to reach his mining
 * station before us, true otherwise
 */
float targetother[3];
float target[3];
float v1[3];
float v2[3];

mathVecSubtract(targetother,cother,other,3);
mathVecSubtract(target,c,state,3);

mathVecMul (other+3,other+3,2.05f);
mathVecMul (state+3,other+3,2.05f);
mathVecSubtract (v1,targetother,other+3,3);
mathVecSubtract (v2,target,state+3,3);

return ((mathVecMagnitude(v1, 3)+ 0.14f)>(mathVecMagnitude(v2, 3)));
//END::PROC::getInTime
}
static void spin (float c[3], float mystate[12])
{
//BEGIN::PROC::spin
/*
 * mmh.. Do I really need to
 * exaplain what this function
 * is for?
 */
float v[3];
char i;
float cosalfa=0;

if (ok<3) {
	PgetAsteroidNormal(v);
	
	for(i=0;i<3;i++)
		cosalfa+=v[i]*mystate[i+6];
		
	if(cosalfa<0)
		mathVecMul (v,v,-1);
			
	ZRSetAttitudeTarget(v);
	
	if (fabs(cosalfa)>0.986f)
		ok++;
}
else {
	v[0]=(0.5236f - mystate[9])/50;
	v[1]=0;
	v[2]=0;
	ZRSetTorques(v);
}
 
ZRSetPositionTarget(c); 
//END::PROC::spin
}
static float timeStation (float c[3], float state[12])
{
//BEGIN::PROC::timeStation
/*
 * timeStation return the second
 * at which we can leave the asteroid
 * to go to the mining station and
 * reach it before the end of
 * phase 3
 */
float target[3];
float v2[3],d;

mathVecSubtract(target,c,state,3);

mathVecMul (state+3,state+3,2.85);
mathVecSubtract (v2,target,state+3,3);

d=mathVecMagnitude (v2,3);
d*=26;
return 180-(d+6);
//END::PROC::timeStation
}
static void revolve (float c[3], float mystate[12])
{
//BEGIN::PROC::revolve
/*
 * Revolving function
 */
float v[4];
float q[4];
float a[4];
float a2[4];
float asse[3];
float d;
float b;
float b2;


PgetAsteroidNormal(asse);

v[3]=0;

mathVecSubtract (v,c,mystate,3);

d=mathVecMagnitude(v, 3);
if(d >= 0.3f) {
	d=0.2355f/d; // (pi/2)*(0.3/distance)
}
else {
	d= 1.57f - 2.615f*d; // pi - (pi/2)*(distance/0.25)
}

// This creates the quaternion
// (yeah, a quaternion, that's cool right?)
b=d-1.57f;
b2=b*b*b;
q[3]=(-b + b2/6 - b2*b*b/120);
b2=d*d*d;
mathVecMul (q,asse,sp*(d - b2/6 + b2*d*d/120)); 

// v[] is rotated according to the quaternion
funz(a2,q,v);
mathVecMul(q,q,-1);
funz(a,a2,q);

//then it's normalized and converted in velocity
mathVecNormalize(a, 3);
mathVecMul(a,a,0.0232f);
ZRSetVelocityTarget(a);


q[3]=1;
mathVecMul(q,asse,-sp*0.043619f);

funz(a2,q,v);
mathVecMul(q,q,-1);
funz(a,a2,q);

ZRSetAttitudeTarget(a);

//END::PROC::revolve
}
static float checkTarget (float other[12], float c[3])
{
//BEGIN::PROC::checkTarget
/*
 * checkTarget returns a value 
 * between 1 and 2: an high 
 * value means that the player
 * "other" is headed towards
 * the point "c".
 */
float v[3]={0,0,0};
float a[3];

mathVecAdd (v,v,other+3,3);
mathVecSubtract (a,c,other,3);
mathVecNormalize(v, 3);
mathVecNormalize(a, 3);
mathVecAdd (a,a,v,3);

return mathVecMagnitude (a,3);

//END::PROC::checkTarget
}
static void meltIce (float state[12], float time)
{
//BEGIN::PROC::meltIce
/*
 * meltIce melts the
 * ice sheet on Opulens
 * (no shit!)
 */float asse[4];
float q[4];
float a2[4]; 
float x[3]={0,-0.35,-0.2};
int i;
asse[3]=0;
if (act==2) {
	if (PiceHits()+PotherIceHits()>15 || time<47) {
		revolve (x,state);
	}
	else {
		ZRSetPositionTarget (state);
		mathVecSubtract(x,x,state,3);
		ZRSetAttitudeTarget(x);	
	}
}
else {
	mathVecSubtract (asse,x,state,3);
	mathVecNormalize(asse, 3);
	mathVecCross(q, state+6, asse); 
	q[3]=0;
	for(i=0;i<3;i++)
		q[3]+=state[i+6]*asse[i];
	if(q[3]>0.51f) {
		funz(a2,q,asse);
		mathVecMul(q,q,-1);
		funz(asse,a2,q);
	}   
	ZRSetAttitudeTarget(asse);
	PgetAsteroidNormal(asse);
	mathVecMul(asse, asse, 0.01f*(80-time));
	mathVecAdd(asse, x, asse, 3);
	SetPosFaster(asse,state);
}
if (time>59)
	Plaser();
//END::PROC::meltIce
}
