

# Class tarte::WebsterFDTD

**template &lt;typename ftype, int kMaxN&gt;**



[**ClassList**](annotated.md) **>** [**tarte**](namespacetarte.md) **>** [**WebsterFDTD**](classtarte_1_1_webster_f_d_t_d.md)




















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**FrequencyResponse**](structtarte_1_1_webster_f_d_t_d_1_1_frequency_response.md) <br> |
| struct | [**PolesResidues**](structtarte_1_1_webster_f_d_t_d_1_1_poles_residues.md) <br> |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  std::vector&lt; [**FrequencyResponse**](structtarte_1_1_webster_f_d_t_d_1_1_frequency_response.md) &gt; | [**ComputeFrequencyResponse**](#function-computefrequencyresponse) (const std::vector&lt; double &gt; & frequenciesHz) <br> |
|  [**PolesResidues**](structtarte_1_1_webster_f_d_t_d_1_1_poles_residues.md) | [**ComputePolesResidues**](#function-computepolesresidues) () <br> |
|  void | [**DspSetup**](#function-dspsetup) (ftype sampleRate, [**Articulation**](classtarte_1_1_articulation.md) \* art=nullptr) <br> |
|  std::tuple&lt; ftype, ftype &gt; | [**GetIOLinearDependencyCoefficients**](#function-getiolineardependencycoefficients) () <br> |
|  void | [**Process**](#function-process) (ftype inputFlow, ftype outputFlow=0) <br> |
|  ArrayN | [**ReadCurrentDensityDistribution**](#function-readcurrentdensitydistribution) () <br> |
|  ArrayNm1 | [**ReadCurrentVelocityDistribution**](#function-readcurrentvelocitydistribution) () <br> |
|  ftype | [**ReadInputPressure**](#function-readinputpressure) () <br> |
|  ftype | [**ReadPowerDissipatedRadiation**](#function-readpowerdissipatedradiation) () <br> |
|  ftype | [**ReadPowerDissipatedWalls**](#function-readpowerdissipatedwalls) () <br> |
|  ftype | [**ReadPowerExchanged**](#function-readpowerexchanged) () <br> |
|  ftype | [**ReadPowerStoredFluid**](#function-readpowerstoredfluid) () <br> |
|  ftype | [**ReadPowerStoredFluidKinetic**](#function-readpowerstoredfluidkinetic) () <br> |
|  ftype | [**ReadPowerStoredFluidPotential**](#function-readpowerstoredfluidpotential) () <br> |
|  ftype | [**ReadPowerStoredRadiation**](#function-readpowerstoredradiation) () <br> |
|  ftype | [**ReadPowerStoredWalls**](#function-readpowerstoredwalls) () <br> |
|  ftype | [**ReadPowerTotal**](#function-readpowertotal) () <br> |
|  ftype | [**ReadRadiatedPressure**](#function-readradiatedpressure) () <br> |
|  void | [**SetConstantSection**](#function-setconstantsection) (ftype section) <br> |
|  void | [**SetTargetGeometry**](#function-settargetgeometry) (intype const \* in, std::size\_t const size) <br> |
|  void | [**SetTargetGeometryFromArticulation**](#function-settargetgeometryfromarticulation) ([**Articulation**](classtarte_1_1_articulation.md) articulation, bool force\_direct=false) <br> |
|   | [**WebsterFDTD**](#function-websterfdtd) (ftype sampleRate, ftype length=ftype(17e-2), [**Articulation**](classtarte_1_1_articulation.md) \* art=nullptr) <br> |
|  void | [**getTargetGeometry**](#function-gettargetgeometry) (ftype \* out, const size\_t N) <br> |
|  std::size\_t | [**get\_N**](#function-get_n) () <br> |
|  ftype | [**get\_c0**](#function-get_c0) () <br> |
|  ftype | [**get\_l0**](#function-get_l0) () <br> |
|  float | [**get\_lpf\_frequency**](#function-get_lpf_frequency) () <br> |
|  ftype | [**get\_rho0**](#function-get_rho0) () <br> |
|  ftype | [**get\_wall\_area\_damping**](#function-get_wall_area_damping) () <br> |
|  ftype | [**get\_wall\_area\_mass**](#function-get_wall_area_mass) () <br> |
|  void | [**initializeFilters**](#function-initializefilters) () <br> |
|  void | [**set\_N\_update\_geometry**](#function-set_n_update_geometry) (const int NUpdateGeometry) <br> |
|  void | [**set\_c0**](#function-set_c0) (const ftype sound\_velocity) <br> |
|  void | [**set\_compute\_powers**](#function-set_compute_powers) (const bool compute\_powers) <br> |
|  void | [**set\_l0**](#function-set_l0) (const ftype length) <br> |
|  void | [**set\_lp\_Q**](#function-set_lp_q) (int index, ftype Q) <br> |
|  void | [**set\_lp\_Qs**](#function-set_lp_qs) (ftype Q) <br> |
|  void | [**set\_lp\_frequencies**](#function-set_lp_frequencies) (ftype freq) <br> |
|  void | [**set\_lp\_frequency**](#function-set_lp_frequency) (int index, ftype freq) <br> |
|  void | [**set\_pumped\_flow**](#function-set_pumped_flow) (const bool pumpedFlow) <br> |
|  void | [**set\_radiation**](#function-set_radiation) (const bool isRadiating) <br> |
|  void | [**set\_rho0**](#function-set_rho0) (const ftype rest\_density) <br> |
|  void | [**set\_time\_varying\_geometry**](#function-set_time_varying_geometry) (const bool isVarying) <br> |
|  void | [**set\_yielding\_walls**](#function-set_yielding_walls) (const bool isYielding) <br> |




























## Public Functions Documentation




### function ComputeFrequencyResponse 

```C++
std::vector< FrequencyResponse > tarte::WebsterFDTD::ComputeFrequencyResponse (
    const std::vector< double > & frequenciesHz
) 
```




<hr>



### function ComputePolesResidues 

```C++
PolesResidues tarte::WebsterFDTD::ComputePolesResidues () 
```




<hr>



### function DspSetup 

```C++
void tarte::WebsterFDTD::DspSetup (
    ftype sampleRate,
    Articulation * art=nullptr
) 
```




<hr>



### function GetIOLinearDependencyCoefficients 

```C++
std::tuple< ftype, ftype > tarte::WebsterFDTD::GetIOLinearDependencyCoefficients () 
```




<hr>



### function Process 

```C++
void tarte::WebsterFDTD::Process (
    ftype inputFlow,
    ftype outputFlow=0
) 
```




<hr>



### function ReadCurrentDensityDistribution 

```C++
inline ArrayN tarte::WebsterFDTD::ReadCurrentDensityDistribution () 
```




<hr>



### function ReadCurrentVelocityDistribution 

```C++
inline ArrayNm1 tarte::WebsterFDTD::ReadCurrentVelocityDistribution () 
```




<hr>



### function ReadInputPressure 

```C++
inline ftype tarte::WebsterFDTD::ReadInputPressure () 
```




<hr>



### function ReadPowerDissipatedRadiation 

```C++
inline ftype tarte::WebsterFDTD::ReadPowerDissipatedRadiation () 
```




<hr>



### function ReadPowerDissipatedWalls 

```C++
inline ftype tarte::WebsterFDTD::ReadPowerDissipatedWalls () 
```




<hr>



### function ReadPowerExchanged 

```C++
inline ftype tarte::WebsterFDTD::ReadPowerExchanged () 
```




<hr>



### function ReadPowerStoredFluid 

```C++
inline ftype tarte::WebsterFDTD::ReadPowerStoredFluid () 
```




<hr>



### function ReadPowerStoredFluidKinetic 

```C++
inline ftype tarte::WebsterFDTD::ReadPowerStoredFluidKinetic () 
```




<hr>



### function ReadPowerStoredFluidPotential 

```C++
inline ftype tarte::WebsterFDTD::ReadPowerStoredFluidPotential () 
```




<hr>



### function ReadPowerStoredRadiation 

```C++
inline ftype tarte::WebsterFDTD::ReadPowerStoredRadiation () 
```




<hr>



### function ReadPowerStoredWalls 

```C++
inline ftype tarte::WebsterFDTD::ReadPowerStoredWalls () 
```




<hr>



### function ReadPowerTotal 

```C++
inline ftype tarte::WebsterFDTD::ReadPowerTotal () 
```




<hr>



### function ReadRadiatedPressure 

```C++
inline ftype tarte::WebsterFDTD::ReadRadiatedPressure () 
```




<hr>



### function SetConstantSection 

```C++
void tarte::WebsterFDTD::SetConstantSection (
    ftype section
) 
```




<hr>



### function SetTargetGeometry 

```C++
template<typename intype>
inline void tarte::WebsterFDTD::SetTargetGeometry (
    intype const * in,
    std::size_t const size
) 
```




<hr>



### function SetTargetGeometryFromArticulation 

```C++
void tarte::WebsterFDTD::SetTargetGeometryFromArticulation (
    Articulation articulation,
    bool force_direct=false
) 
```




<hr>



### function WebsterFDTD 

```C++
tarte::WebsterFDTD::WebsterFDTD (
    ftype sampleRate,
    ftype length=ftype(17e-2),
    Articulation * art=nullptr
) 
```




<hr>



### function getTargetGeometry 

```C++
inline void tarte::WebsterFDTD::getTargetGeometry (
    ftype * out,
    const size_t N
) 
```




<hr>



### function get\_N 

```C++
inline std::size_t tarte::WebsterFDTD::get_N () 
```




<hr>



### function get\_c0 

```C++
inline ftype tarte::WebsterFDTD::get_c0 () 
```




<hr>



### function get\_l0 

```C++
inline ftype tarte::WebsterFDTD::get_l0 () 
```




<hr>



### function get\_lpf\_frequency 

```C++
inline float tarte::WebsterFDTD::get_lpf_frequency () 
```




<hr>



### function get\_rho0 

```C++
inline ftype tarte::WebsterFDTD::get_rho0 () 
```




<hr>



### function get\_wall\_area\_damping 

```C++
inline ftype tarte::WebsterFDTD::get_wall_area_damping () 
```




<hr>



### function get\_wall\_area\_mass 

```C++
inline ftype tarte::WebsterFDTD::get_wall_area_mass () 
```




<hr>



### function initializeFilters 

```C++
void tarte::WebsterFDTD::initializeFilters () 
```




<hr>



### function set\_N\_update\_geometry 

```C++
inline void tarte::WebsterFDTD::set_N_update_geometry (
    const int NUpdateGeometry
) 
```




<hr>



### function set\_c0 

```C++
inline void tarte::WebsterFDTD::set_c0 (
    const ftype sound_velocity
) 
```




<hr>



### function set\_compute\_powers 

```C++
inline void tarte::WebsterFDTD::set_compute_powers (
    const bool compute_powers
) 
```




<hr>



### function set\_l0 

```C++
inline void tarte::WebsterFDTD::set_l0 (
    const ftype length
) 
```




<hr>



### function set\_lp\_Q 

```C++
void tarte::WebsterFDTD::set_lp_Q (
    int index,
    ftype Q
) 
```




<hr>



### function set\_lp\_Qs 

```C++
void tarte::WebsterFDTD::set_lp_Qs (
    ftype Q
) 
```




<hr>



### function set\_lp\_frequencies 

```C++
void tarte::WebsterFDTD::set_lp_frequencies (
    ftype freq
) 
```




<hr>



### function set\_lp\_frequency 

```C++
void tarte::WebsterFDTD::set_lp_frequency (
    int index,
    ftype freq
) 
```




<hr>



### function set\_pumped\_flow 

```C++
inline void tarte::WebsterFDTD::set_pumped_flow (
    const bool pumpedFlow
) 
```




<hr>



### function set\_radiation 

```C++
inline void tarte::WebsterFDTD::set_radiation (
    const bool isRadiating
) 
```




<hr>



### function set\_rho0 

```C++
inline void tarte::WebsterFDTD::set_rho0 (
    const ftype rest_density
) 
```




<hr>



### function set\_time\_varying\_geometry 

```C++
inline void tarte::WebsterFDTD::set_time_varying_geometry (
    const bool isVarying
) 
```




<hr>



### function set\_yielding\_walls 

```C++
inline void tarte::WebsterFDTD::set_yielding_walls (
    const bool isYielding
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `include/websterFDTD.h`

