﻿Module CS_Model_Environment !Module of environment
    type DailyEnvironment_type
        
        real :: TMin_ = 0 !
        real :: TMax_ = 0 ! 
        integer :: Hour_ = 0 ! 
        
    contains
    
        procedure, pass (this) :: getTMin
        procedure, pass (this) :: setTMin
        procedure, pass (this) :: getTMax
        procedure, pass (this) :: setTMax
        procedure, pass (this) :: getHour
        procedure, pass (this) :: setHour
        procedure, pass (this) :: fetchTemperature
    
    end Type DailyEnvironment_type
    
    ! interface to reference the constructor
    interface DailyEnvironment_type
        module procedure DailyEnvironment_type_constructor
    end interface DailyEnvironment_type
    
    contains
    
    ! constructor for the type
    type (DailyEnvironment_type) function DailyEnvironment_type_constructor(TMin, TMax, Hour)
        implicit none
        real, intent (in) :: TMin, TMax
        integer, intent (in) :: Hour
        DailyEnvironment_type_constructor%TMin_ = TMin
        DailyEnvironment_type_constructor%TMax_ = TMax
        DailyEnvironment_type_constructor%Hour_ = Hour
    end function DailyEnvironment_type_constructor    

    ! get TMin
    real function getTMin(this)
        implicit none
        class (DailyEnvironment_type), intent(in) :: this
        
        getTMin = this%TMin_
    end function getTMin
    
    ! set TMin    
    subroutine setTMin(this, TMin)
        implicit none
        class (DailyEnvironment_type), intent(inout) :: this
        real, intent (in) :: TMin
        
        this%TMin_ = TMin
    end subroutine setTMin
    
    ! get TMax
    real function getTMax(this)
        implicit none
        class (DailyEnvironment_type), intent(in) :: this
        
        getTMax = this%TMax_
    end function getTMax
    
    
    ! set TMax    
    subroutine setTMax(this, TMax)
        implicit none
        class (DailyEnvironment_type), intent(inout) :: this
        real, intent (in) :: TMax
        
        this%TMax_ = TMax
    end subroutine setTMax
    
    ! get Hour
     integer function getHour(this)
         implicit none
         class (DailyEnvironment_type), intent(in) :: this
         
         getHour = this%Hour_
     end function getHour
          
     ! set Hour    
     subroutine setHour(this, Hour)
         implicit none
         class (DailyEnvironment_type), intent(inout) :: this
         integer, intent (in) :: Hour
         
         this%Hour_ = Hour
     end subroutine setHour
    
    ! obtain the temperature accrding to the hour of the day
    ! T(Hour) = Amplitude*sin[w(t - a)] + C.
    ! Amplitude is called the amplitude the height of each peak above the baseline
    ! Hod is Hours Of Day, the period or wavelength (the length of each cycle) 
    ! a  is the phase shift (the horizontal offset of the basepoint; where the curve crosses the baseline as it ascends)
    ! C is average temperature,  the vertical offset (height of the baseline) 
    ! w is the angular frequency, given by w = 2PI/hod 
    real function fetchTemperature(this)
        implicit none
        class (DailyEnvironment_type), intent(in) :: this
        REAL :: Amplitude, C, hod, w, a, g, pi, dawn

        dawn=5
        pi=  4 * atan (1.0_8)
        a = 12                                              ! 12 hours to shift
        hod = 24                                            ! 24 hours a day
        w = (2*pi)/hod
        Amplitude = ((this%TMax_ - this%TMin_)/2)           ! half distance between temperatures
        C = (this%TMin_ + this%TMax_)/2                     ! mean temperature
        g = w*(this%Hour_ - a)
        
        !if (this%Hour_ > dawn) then                         ! if it is later dawn time
            !fetchTemperature = Amplitude*SIN(g)+C           ! calculate temperature acording to the current time
        !else                                                ! else
            !g = w*(dawn - a)
            !fetchTemperature = Amplitude*SIN(g)+C           !calculate temperature like if it was dawn time
        !end if
        
        fetchTemperature = Amplitude*SIN(g)+C           ! calculate temperature acording to the current time

    end function fetchTemperature
    
    ! obtain the Saturation Vapour Pressure (pascals)
    real function fetchSVP(this)
        implicit none
        class (DailyEnvironment_type), intent(in) :: this
        
        fetchSVP = 610.78 * exp( fetchTemperature(this) / ( fetchTemperature(this)  + 238.3 ) * 17.2694 )        !  DA Saturation vapour pressure in pascals: svp = 610.78 *exp( t / ( t + 238.3 ) *17.2694 ) 

    end function fetchSVP
    
    ! obtain the water holding capacity of the air (kg/m3)
    real function fetchWaterHoldingCapacity(this)
        implicit none
        class (DailyEnvironment_type), intent(in) :: this
        
        fetchWaterHoldingCapacity = 0.002166 * fetchSVP(this) / ( fetchTemperature(this) + 273.16 )                                !  DA water holding capacity of the air WHC = 0.002166 * SVP / ( t + 273.16 )   

    end function fetchWaterHoldingCapacity
    
    ! obtain the incoming radiation at the given hour
    !real function fetchIncomingRadiation(this)
    !    implicit none
    !    class (DailyEnvironment_type), intent(in) :: this
    !    
    !    fetchIncomingRadiation = Radiation * (Tmax_- Tmin_)
    !
    !end function fetchIncomingRadiation

        
END Module CS_Model_Environment    
    