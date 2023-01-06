// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef MY_TREADMILLING_FIBER_PROP
#define MY_TREADMILLING_FIBER_PROP

#include "fiber_prop.h"


/// additional Property for MyTreadmillingFiber
/**
 @ingroup Properties
 Assembly is continuous, and can occur at both ends.
 */
class MyTreadmillingFiberProp : public FiberProp
{
    friend class MyTreadmillingFiber;
    
public:
    
    /**
     @defgroup MyTreadmillingFiberPar Parameters of MyTreadmillingFiber
     @ingroup Parameters
     Inherits @ref FiberPar.
     @{
     */
    
    /// see @ref MyTreadmillingFiber
    
    /// Length of discrete units of assembly/disassembly
    real    unit_length;
    
    /// Speed of assembly
    real    growing_speed[2];

    /// Characteristic force for polymer assembly
    real    growing_force[2];
    
    /// Speed of disassembly
    real    shrinking_speed[2];
    
    /// @}
    
private:
    
    /// derived variable:
    real    growing_speed_dt[2];
    /// derived variable:
    real    shrinking_speed_dt[2];
    
public:
    
    /// constructor
    MyTreadmillingFiberProp(const std::string& n) : FiberProp(n) { clear(); }

    /// destructor
    ~MyTreadmillingFiberProp() { }
    
    /// return a Fiber with this property
    Fiber* newFiber() const;
    
    /// set default values
    void clear();
       
    /// set using a Glossary
    void read(Glossary&);
   
    /// check and derive parameter values
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new MyTreadmillingFiberProp(*this); }

    /// write
    void write_values(std::ostream&) const;

    /// print predicted length adn time
    void splash(std::ostream&, real, real) const; 
};

#endif

