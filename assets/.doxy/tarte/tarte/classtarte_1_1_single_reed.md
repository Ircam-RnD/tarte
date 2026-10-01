

# Class tarte::SingleReed

**template &lt;typename ftype&gt;**



[**ClassList**](annotated.md) **>** [**tarte**](namespacetarte.md) **>** [**SingleReed**](classtarte_1_1_single_reed.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**DspSetup**](#function-dspsetup) (float samplerate) <br> |
|  void | [**Process**](#function-process) (float Pin) <br> |
|  ftype | [**ReadAuxiliaryVariable**](#function-readauxiliaryvariable) () <br> |
|  ftype | [**ReadDisplacement**](#function-readdisplacement) () <br> |
|  ftype | [**ReadEffectiveOpening**](#function-readeffectiveopening) () <br> |
|  ftype | [**ReadMeanMouthpieceFlow**](#function-readmeanmouthpieceflow) () <br> |
|  ftype | [**ReadMouthpieceFlow**](#function-readmouthpieceflow) () <br> |
|  ftype | [**ReadPressureDrop**](#function-readpressuredrop) () <br> |
|  ftype | [**ReadRadiatedPressure**](#function-readradiatedpressure) () <br> |
|   | [**SingleReed**](#function-singlereed) (float samplerate) <br> |
|  std::tuple&lt; ftype, ftype, ftype &gt; | [**getCurrentDissipatedPowers**](#function-getcurrentdissipatedpowers) () <br> |
|  std::tuple&lt; ftype, ftype, ftype &gt; | [**getCurrentExchangedPowers**](#function-getcurrentexchangedpowers) () <br> |
|  ftype | [**getCurrentState**](#function-getcurrentstate) () <br> |
|  std::tuple&lt; ftype, ftype, ftype &gt; | [**getCurrentStoredPowers**](#function-getcurrentstoredpowers) () <br> |
|  ftype | [**get\_contact\_stiffness**](#function-get_contact_stiffness) () <br> |
|  ftype | [**get\_damping**](#function-get_damping) () <br> |
|  ftype | [**get\_epsilon\_smooth**](#function-get_epsilon_smooth) () <br> |
|  ftype | [**get\_lay\_position**](#function-get_lay_position) () <br> |
|  ftype | [**get\_mass**](#function-get_mass) () <br> |
|  ftype | [**get\_noise\_ratio**](#function-get_noise_ratio) () <br> |
|  std::shared\_ptr&lt; [**WebsterFDTD**](classtarte_1_1_webster_f_d_t_d.md)&lt; ftype, 300 &gt; &gt; | [**get\_resonator**](#function-get_resonator) () <br> |
|  ftype | [**get\_stiffness**](#function-get_stiffness) () <br> |
|  ftype | [**get\_surface**](#function-get_surface) () <br> |
|  ftype | [**get\_width**](#function-get_width) () <br> |
|  void | [**set\_contact\_stiffness**](#function-set_contact_stiffness) (const ftype & contact\_stiffness) <br> |
|  void | [**set\_damping**](#function-set_damping) (const ftype & damping) <br> |
|  void | [**set\_epsilon\_smooth**](#function-set_epsilon_smooth) (const ftype & epsilon\_smooth) <br> |
|  void | [**set\_lay\_position**](#function-set_lay_position) (const ftype & lay\_position) <br> |
|  void | [**set\_mass**](#function-set_mass) (const ftype & mass) <br> |
|  void | [**set\_noise\_ratio**](#function-set_noise_ratio) (const ftype & noise\_ratio) <br> |
|  void | [**set\_stiffness**](#function-set_stiffness) (const ftype & stiffness) <br> |
|  void | [**set\_surface**](#function-set_surface) (const ftype & surface) <br> |
|  void | [**set\_width**](#function-set_width) (const ftype & width) <br> |




























## Public Functions Documentation




### function DspSetup 

```C++
void tarte::SingleReed::DspSetup (
    float samplerate
) 
```




<hr>



### function Process 

```C++
void tarte::SingleReed::Process (
    float Pin
) 
```




<hr>



### function ReadAuxiliaryVariable 

```C++
inline ftype tarte::SingleReed::ReadAuxiliaryVariable () 
```




<hr>



### function ReadDisplacement 

```C++
inline ftype tarte::SingleReed::ReadDisplacement () 
```




<hr>



### function ReadEffectiveOpening 

```C++
inline ftype tarte::SingleReed::ReadEffectiveOpening () 
```




<hr>



### function ReadMeanMouthpieceFlow 

```C++
inline ftype tarte::SingleReed::ReadMeanMouthpieceFlow () 
```




<hr>



### function ReadMouthpieceFlow 

```C++
inline ftype tarte::SingleReed::ReadMouthpieceFlow () 
```




<hr>



### function ReadPressureDrop 

```C++
inline ftype tarte::SingleReed::ReadPressureDrop () 
```




<hr>



### function ReadRadiatedPressure 

```C++
inline ftype tarte::SingleReed::ReadRadiatedPressure () 
```




<hr>



### function SingleReed 

```C++
tarte::SingleReed::SingleReed (
    float samplerate
) 
```




<hr>



### function getCurrentDissipatedPowers 

```C++
inline std::tuple< ftype, ftype, ftype > tarte::SingleReed::getCurrentDissipatedPowers () 
```




<hr>



### function getCurrentExchangedPowers 

```C++
inline std::tuple< ftype, ftype, ftype > tarte::SingleReed::getCurrentExchangedPowers () 
```




<hr>



### function getCurrentState 

```C++
inline ftype tarte::SingleReed::getCurrentState () 
```




<hr>



### function getCurrentStoredPowers 

```C++
inline std::tuple< ftype, ftype, ftype > tarte::SingleReed::getCurrentStoredPowers () 
```




<hr>



### function get\_contact\_stiffness 

```C++
inline ftype tarte::SingleReed::get_contact_stiffness () 
```




<hr>



### function get\_damping 

```C++
inline ftype tarte::SingleReed::get_damping () 
```




<hr>



### function get\_epsilon\_smooth 

```C++
inline ftype tarte::SingleReed::get_epsilon_smooth () 
```




<hr>



### function get\_lay\_position 

```C++
inline ftype tarte::SingleReed::get_lay_position () 
```




<hr>



### function get\_mass 

```C++
inline ftype tarte::SingleReed::get_mass () 
```




<hr>



### function get\_noise\_ratio 

```C++
inline ftype tarte::SingleReed::get_noise_ratio () 
```




<hr>



### function get\_resonator 

```C++
inline std::shared_ptr< WebsterFDTD < ftype, 300 > > tarte::SingleReed::get_resonator () 
```




<hr>



### function get\_stiffness 

```C++
inline ftype tarte::SingleReed::get_stiffness () 
```




<hr>



### function get\_surface 

```C++
inline ftype tarte::SingleReed::get_surface () 
```




<hr>



### function get\_width 

```C++
inline ftype tarte::SingleReed::get_width () 
```




<hr>



### function set\_contact\_stiffness 

```C++
inline void tarte::SingleReed::set_contact_stiffness (
    const ftype & contact_stiffness
) 
```




<hr>



### function set\_damping 

```C++
inline void tarte::SingleReed::set_damping (
    const ftype & damping
) 
```




<hr>



### function set\_epsilon\_smooth 

```C++
inline void tarte::SingleReed::set_epsilon_smooth (
    const ftype & epsilon_smooth
) 
```




<hr>



### function set\_lay\_position 

```C++
inline void tarte::SingleReed::set_lay_position (
    const ftype & lay_position
) 
```




<hr>



### function set\_mass 

```C++
inline void tarte::SingleReed::set_mass (
    const ftype & mass
) 
```




<hr>



### function set\_noise\_ratio 

```C++
inline void tarte::SingleReed::set_noise_ratio (
    const ftype & noise_ratio
) 
```




<hr>



### function set\_stiffness 

```C++
inline void tarte::SingleReed::set_stiffness (
    const ftype & stiffness
) 
```




<hr>



### function set\_surface 

```C++
inline void tarte::SingleReed::set_surface (
    const ftype & surface
) 
```




<hr>



### function set\_width 

```C++
inline void tarte::SingleReed::set_width (
    const ftype & width
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `include/single_reed.h`

