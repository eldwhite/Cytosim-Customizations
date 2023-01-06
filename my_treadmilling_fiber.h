// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef MY_TREADMILLING_FIBER_H
#define MY_TREADMILLING_FIBER_H

#include "sim.h"
#include "vector.h"
#include "node_list.h"
#include "fiber.h"

class MyTreadmillingFiberProp;


/// A Fiber that grows for a period of time before assembly stops and disassembly begins 
/**
 This is not documented yet!

 See the @ref MyTreadmillingFiberPar.

 @todo Document MyTreadmillingFiber
 @ingroup FiberGroup
 */
class MyTreadmillingFiber : public Fiber
{   
private:
    /// assembly during last time-step 
    real 	mGrowthM;
    real 	mGrowthP;

    /// countdown timers for MINUS_END
    real 	    treadmillOnsetTime;
    /// countdown timers for PLUS_END
    real 	    disassemblyOnsetTime;
    
    unsigned    unitP[3];

    /// dynamic state of PLUS_END
    state_t     mStateP;
    /// assembly during last time-step
    /// real        mGrowthP;
    
    unsigned    unitM[3];
    
    /// dynamic state of MINUS_END
    state_t     mStateM;
    /// assembly during last time-step
    /// real        mGrowthM;
 
    /// calculate state at MINUS_END
    state_t	calculateStateM() const;

    /// calculate state at PLUS_END
    state_t	calculateStateP() const;
public:
    
    /// the Property of this object
    MyTreadmillingFiberProp const* prop;
  
    /// constructor
    MyTreadmillingFiber(MyTreadmillingFiberProp const*);

    /// destructor
    virtual ~MyTreadmillingFiber();
        
    //--------------------------------------------------------------------------
    
    /// return assembly/disassembly state of MINUS_END
    state_t     dynamicStateM() const; ///{ return mStateM; }
    
    /// return assembly/disassembly state of PLUS_END
    state_t     dynamicStateP() const; ///{ return mStateP; }
    
    /// change state of MINUS_END
    void        setDynamicStateM(state_t s);
    
    /// change state of PLUS_END
    void        setDynamicStateP(state_t s);
    
    /// length increment at MINUS_END during last time-step
    real        freshAssemblyM() const { return mGrowthM; }
    
    /// length increment at PLUS_END during last time-step
    real        freshAssemblyP() const { return mGrowthP; }
    
    //--------------------------------------------------------------------------
    
    /// simulat instability of PLUS_END
    int 	stepPlusEnd(); 
    /// simulate instability of MINUS_END
    int 	stepMinusEnd(); 

    /// monte-carlo step
    void        step();
    
    //--------------------------------------------------------------------------
    
    /// write to Outputter
    void        write(Outputter&) const;

    /// read from Inputter
    void        read(Inputter&, Simul&, ObjectTag);
    
};


#endif
