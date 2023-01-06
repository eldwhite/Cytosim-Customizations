// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "assert_macro.h"
#include "my_treadmilling_fiber.h"
#include "my_treadmilling_fiber_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"
#include "space.h"


//------------------------------------------------------------------------------
MyTreadmillingFiber::MyTreadmillingFiber(MyTreadmillingFiberProp const* p) : Fiber(p), prop(p)
{

//    /// set PLUS_END as growing
//    unitP[0] = 1;
//    unitP[1] = 1;
//    unitP[2] = 1;
    mStateP = STATE_GREEN; 
    mGrowthP = 0; 

    treadmillOnsetTime = RNG.exponential()+7;

//    /// set MINUS_END as static
//    unitM[0] = 0;
//    unitM[1] = 0;
//    unitM[2] = 0;
    mStateM = STATE_WHITE; 
    mGrowthM = 0; 

    /// get time to start treadmilling, and start disassembly	
    ///treadmillOnsetTime = RNG.exponential()+7;
    disassemblyOnsetTime = RNG.exponential()+7;
}

MyTreadmillingFiber::~MyTreadmillingFiber()
{
    prop = nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -

/// Minus end state
state_t MyTreadmillingFiber::calculateStateM() const
{
    if ( treadmillOnsetTime < 0 )
	{
        return STATE_RED;
	}
	else
	{
        return STATE_WHITE;
	}
}

state_t MyTreadmillingFiber::dynamicStateM() const
{
//    calculateStateM();
//    assert_true( mStateM == calculateStateM());
    return mStateM;
}

/// Check for err minus end
void MyTreadmillingFiber::setDynamicStateM(state_t s)
{

    if ( s == STATE_WHITE || s == STATE_RED )
        mStateM = s;

//    else if ( s != mStateM)
//    {
//        mStateM     = s;
//        unitM[1]    = (4 - s ) / 2;
//        unitM[0]    = (4 - s - 2*unitM[1] );
//        assert_true( 0==unitP[0] || unitP[0]==1 );
//        assert_true( 0==unitP[1] || unitP[1]==1 );
//        assert_true( mStateM == calculateStateM() );
//    }
//    }
    else
        throw InvalidParameter("Invalid State for MyTreadmillingFiber MINUS_END");
}

/// Update
int MyTreadmillingFiber::stepMinusEnd()
{
//    int res = 0;
	constexpr int P = 0, M = 1;
    ///Approach treadmilling onset
    treadmillOnsetTime -= 0.001;
    mStateM = calculateStateM();
    return treadmillOnsetTime;
//
//	if (mStateM == STATE_WHITE )
//	{
//		mGrowthM = 0;
//	}
//	else if ( mStateM == STATE_RED)
//	{
//		mGrowthM = prop->shrinking_speed_dt[M];
//	}
////    return res;
} 

//-----------------------------------------------------------------------------

/// Plus end state - initiate stopping assembly
state_t MyTreadmillingFiber::calculateStateP() const
{ 

	if ( disassemblyOnsetTime <= 0 )
	{
        return STATE_WHITE;
	}
	else
	{
        return STATE_GREEN;
	} 
}

state_t MyTreadmillingFiber::dynamicStateP() const
{
//    mStateP = calculateStateP();
//	assert_true( mStateP == calculateStateP() );
	return mStateP;
}

/// Check for err plus end
void MyTreadmillingFiber::setDynamicStateP(state_t s)
{
    
    if ( s == STATE_WHITE || s == STATE_GREEN )
        mStateP = s;
//    else if ( s ! mStateP)
//        mStateP     = s;
//        assert_true( mStateP == calculateStateP() );
    else
        throw InvalidParameter("Invalid State for MyTreadmillingFiber PLUS_END");
}

/// Update
int MyTreadmillingFiber::stepPlusEnd()
{
//    res = 0;
	constexpr int P = 0, M = 1;
	
    disassemblyOnsetTime -= 0.001;
    mStateP = calculateStateP();
    
    return disassemblyOnsetTime;
    
//	///If it's growing - prior to disassembly onset
//	if ( mStateP == STATE_GREEN )
//	{
//		real forceP = projectedForceEndP();
//		mGrowthP = prop->growing_speed_dt[P] * prop->free_polymer;
//
//        assert_true(mGrowthP >=0);
//
//        if ( forceP < 0 && prop->growing_force[P] < INFINITY )
//			mGrowthP *= exp(forceP/prop->growing_force[P]);
//	}
//    else if (mStateP == STATE_RED)
//    {
//        mGrowthP = prop ->shrinking_speed_dt[P];
//    }
//	else
//	{
//		mGrowthP = 0;
//	}
////    return res;
}

//------------------------------------------------------------------------------
void MyTreadmillingFiber::step()
{
    constexpr int P = 0, M = 1; 
    int incP = stepPlusEnd(); 
    int incM = stepMinusEnd();
    
    if ( mStateP == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceP = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        mGrowthP = prop->growing_speed_dt[P] * prop->free_polymer;
        
        assert_true(mGrowthP>=0);
        // antagonistic force (< 0) decreases assembly rate exponentially
        if ( forceP < 0  &&  prop->growing_force[P] < INFINITY )
            mGrowthP *= exp(forceP/prop->growing_force[P]);
    }
    else if ( mStateP == STATE_RED )
    {
        mGrowthP = prop->shrinking_speed_dt[P];
    }
    else
    {
        mGrowthP = 0;
    }
    
    
    if ( mStateM == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceM = projectedForceEndM();
        
        // growth is reduced if free monomers are scarce:
        mGrowthM = prop->growing_speed_dt[M] * prop->free_polymer;

        assert_true(mGrowthM>=0);
        // antagonistic force (< 0) decreases assembly rate exponentially
        if ( forceM < 0  &&  prop->growing_force[M] < INFINITY  )
            mGrowthM *= exp(forceM/prop->growing_force[M]);
    }
    else if ( mStateM == STATE_RED )
    {
        mGrowthM = prop->shrinking_speed_dt[M];
    }
    else
    {
        mGrowthM = 0;
    }
    
    
    real len = length();
    real inc = mGrowthP + mGrowthM;
    if ( inc < 0  &&  len + inc < prop->min_length )
    {
        // the fiber is too short, we delete it:
        delete(this);
        return;
    }
    else if ( len + inc < prop->max_length )
    {
        if ( mGrowthM ) growM(mGrowthM);
        if ( mGrowthP ) growP(mGrowthP);
    }
    else if ( len < prop->max_length )
    {
        // the remaining possible growth is distributed to the two ends:
        inc = ( prop->max_length - len ) / inc;
        if ( mGrowthM ) growM(inc*mGrowthM);
        if ( mGrowthP ) growP(inc*mGrowthP);
    }


//    mGrowthP = incP * prop->unit_length;
//    mGrowthM = incM * prop->unit_length;
//
//    if ( incM || incP )
//    {
//        if ( length() + mGrowthM + mGrowthP < prop->min_length )
//        {
//            if ( !prop->persistent )
//            {
//                delete(this);
//                return;
//            }
//        }
//        else if ( length() + mGrowthM + mGrowthP < prop->max_length )
//        {
//            if ( incP ) growP(mGrowthP);
//            if ( incM ) growM(mGrowthM);
//        }
//    }

    Fiber::step(); 

}

                  
//------------------------------------------------------------------------------
#pragma mark -


void MyTreadmillingFiber::write(Outputter& out) const
{
    Fiber::write(out);

    /// write variables describing the dynamic state of the ends:
    writeHeader(out, TAG_DYNAMIC);
    out.writeUInt16(mStateM); 
    out.writeUInt16(mStateP); 
//    out.writeUInt16(unitM[0]);
//    out.writeUInt16(unitM[1]);
//    out.writeUInt16(unitP[0]);
//    out.writeUInt16(unitP[1]);
}

void MyTreadmillingFiber::read(Inputter& in, Simul& sim, ObjectTag tag)
{
#ifdef BACKWARD_COMPATIBILITY
    if ( tag == TAG_DYNAMIC || ( tag == TAG && in.formatID() < 44 ) )
#else
    if ( tag == TAG_DYNAMIC )
#endif
    {
//        mStateM = calculateStateM();
//        mStateP = calculateStateP();
//	unitM[0] = in.readUInt16();
//	unitM[1] = in.readUInt16();  
	mStateM = in.readUInt16();

//	unitP[0] = in.readUInt16();
//    unitP[1] = in.readUInt16();
	mStateP = in.readUInt16();
    } 
        
#ifdef BACKWARD_COMPATIBILITY
    if ( tag != TAG_DYNAMIC || in.formatID() < 44 )
#else
    else
#endif
        Fiber::read(in, sim, tag);
}

