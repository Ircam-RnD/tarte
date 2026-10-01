

# Class tarte::BodyCoverPair

**template &lt;typename ftype&gt;**



[**ClassList**](annotated.md) **>** [**tarte**](namespacetarte.md) **>** [**BodyCoverPair**](classtarte_1_1_body_cover_pair.md)






















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Eigen::Vector&lt; ftype, half\_N &gt; | [**half\_state\_type**](#typedef-half_state_type)  <br> |
| typedef ftype | [**scalar\_type**](#typedef-scalar_type)  <br> |
| typedef Eigen::Vector&lt; ftype, N &gt; | [**state\_type**](#typedef-state_type)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BodyCoverPair**](#function-bodycoverpair) () <br> |
|  ftype | [**ComputeFlow**](#function-computeflow) (const ftype & Psub, const ftype & Psup) <br> |
|  void | [**EffectiveAreas**](#function-effectiveareas) (state\_type & out\_area\_P\_sub, state\_type & out\_area\_P\_sup) <br> |
|  ftype | [**Enl**](#function-enl) (const state\_type & state\_q, bool recompute\_intermediary=false) <br> |
|  void | [**FillIntermediary**](#function-fillintermediary) (const state\_type & state\_q) <br> |
|  void | [**Fnl**](#function-fnl) (const state\_type & state\_q, state\_type & out, bool recompute\_intermediary=false) <br> |
|  void | [**KOp**](#function-kop) (const state\_type & state\_q, state\_type & out) <br> |
|  state\_type | [**MinvOp**](#function-minvop) (const state\_type & state\_p) <br> |
|  void | [**ROp**](#function-rop) (const state\_type & state\_p, state\_type & out) <br> |
|  half\_state\_type | [**ReadEffectiveOpenings**](#function-readeffectiveopenings) () <br> |
|  ftype | [**get\_alpha\_contact\_stiffness**](#function-get_alpha_contact_stiffness) () <br> |
|  ftype | [**get\_c0**](#function-get_c0) () <br> |
|  ftype | [**get\_contact\_stiffness**](#function-get_contact_stiffness) () <br> |
|  ftype | [**get\_epsilon\_smooth**](#function-get_epsilon_smooth) () <br> |
|  ftype | [**get\_eta\_stiffness**](#function-get_eta_stiffness) (FoldIdentifier fold\_id=kBoth) <br> |
|  Gender | [**get\_gender**](#function-get_gender) () <br> |
|  ftype | [**get\_kt**](#function-get_kt) () <br> |
|  half\_state\_type | [**get\_lengths**](#function-get_lengths) (FoldIdentifier fold\_id=kBoth) <br> |
|  half\_state\_type | [**get\_masses**](#function-get_masses) (FoldIdentifier fold\_id=kBoth) <br> |
|  ftype | [**get\_noise\_ratio**](#function-get_noise_ratio) () <br> |
|  half\_state\_type | [**get\_rest\_positions**](#function-get_rest_positions) (FoldIdentifier fold\_id=kBoth) <br> |
|  ftype | [**get\_rho0**](#function-get_rho0) () <br> |
|  Eigen::Vector&lt; ftype, 4 &gt; | [**get\_stiffnesses**](#function-get_stiffnesses) (FoldIdentifier fold\_id=kBoth) <br> |
|  half\_state\_type | [**get\_thicknesses**](#function-get_thicknesses) (FoldIdentifier fold\_id=kBoth) <br> |
|  void | [**set\_alpha\_contact\_stiffness**](#function-set_alpha_contact_stiffness) (const ftype & alpha\_contact\_stiffness) <br> |
|  void | [**set\_c0**](#function-set_c0) (const ftype & c0) <br> |
|  void | [**set\_contact\_stiffness**](#function-set_contact_stiffness) (const ftype & contact\_stiffness) <br> |
|  void | [**set\_epsilon\_smooth**](#function-set_epsilon_smooth) (const ftype & epsilon\_smooth) <br> |
|  void | [**set\_eta\_stiffness**](#function-set_eta_stiffness) (const ftype & eta\_stiffness, FoldIdentifier fold\_id=kBoth) <br> |
|  void | [**set\_gender**](#function-set_gender) (Gender gender) <br> |
|  void | [**set\_kt**](#function-set_kt) (const ftype & kt) <br> |
|  void | [**set\_length**](#function-set_length) (const ftype & length, FoldIdentifier fold\_id=kBoth) <br> |
|  void | [**set\_masses**](#function-set_masses-12) (const half\_state\_type & masses, FoldIdentifier fold\_id=kBoth) <br> |
|  void | [**set\_masses**](#function-set_masses-22) (const ftype & lower, const ftype & upper, const ftype & body, FoldIdentifier fold\_id=kBoth) <br> |
|  void | [**set\_muscles\_activation**](#function-set_muscles_activation) (const ftype & ct\_activity, const ftype & ta\_activity, const ftype & lc\_activity, FoldIdentifier fold\_id=kBoth) <br> |
|  void | [**set\_noise\_ratio**](#function-set_noise_ratio) (const ftype & noise\_ratio) <br> |
|  void | [**set\_rest\_positions**](#function-set_rest_positions-12) (const half\_state\_type & rest\_positions, FoldIdentifier fold\_id=kBoth) <br>_--------------------------------_  _Control functions --------------------------------------_ _///_ |
|  void | [**set\_rest\_positions**](#function-set_rest_positions-22) (const ftype & lower, const ftype & upper, const ftype & body, FoldIdentifier fold\_id=kBoth) <br> |
|  void | [**set\_rho0**](#function-set_rho0) (const ftype & rho0) <br> |
|  void | [**set\_stiffnesses**](#function-set_stiffnesses-12) (const Eigen::Vector&lt; ftype, 4 &gt; & stiffnesses, FoldIdentifier fold\_id=kBoth) <br> |
|  void | [**set\_stiffnesses**](#function-set_stiffnesses-22) (const ftype & lower, const ftype & upper, const ftype & body, const ftype & coupling, FoldIdentifier fold\_id=kBoth) <br> |
|  void | [**set\_thicknesses**](#function-set_thicknesses-12) (const half\_state\_type & thicknesses, FoldIdentifier fold\_id=kBoth) <br> |
|  void | [**set\_thicknesses**](#function-set_thicknesses-22) (const ftype & lower, const ftype & upper, const ftype & body, FoldIdentifier fold\_id=kBoth) <br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  int | [**get\_N**](#function-get_n) () <br> |
|  int | [**get\_half\_N**](#function-get_half_n) () <br> |


























## Public Types Documentation




### typedef half\_state\_type 

```C++
using tarte::BodyCoverPair< ftype >::half_state_type = Eigen::Vector<ftype, half_N>;
```




<hr>



### typedef scalar\_type 

```C++
using tarte::BodyCoverPair< ftype >::scalar_type = ftype;
```




<hr>



### typedef state\_type 

```C++
using tarte::BodyCoverPair< ftype >::state_type = Eigen::Vector<ftype, N>;
```




<hr>
## Public Functions Documentation




### function BodyCoverPair 

```C++
inline tarte::BodyCoverPair::BodyCoverPair () 
```




<hr>



### function ComputeFlow 

```C++
inline ftype tarte::BodyCoverPair::ComputeFlow (
    const ftype & Psub,
    const ftype & Psup
) 
```




<hr>



### function EffectiveAreas 

```C++
inline void tarte::BodyCoverPair::EffectiveAreas (
    state_type & out_area_P_sub,
    state_type & out_area_P_sup
) 
```




<hr>



### function Enl 

```C++
inline ftype tarte::BodyCoverPair::Enl (
    const state_type & state_q,
    bool recompute_intermediary=false
) 
```




<hr>



### function FillIntermediary 

```C++
inline void tarte::BodyCoverPair::FillIntermediary (
    const state_type & state_q
) 
```




<hr>



### function Fnl 

```C++
inline void tarte::BodyCoverPair::Fnl (
    const state_type & state_q,
    state_type & out,
    bool recompute_intermediary=false
) 
```




<hr>



### function KOp 

```C++
inline void tarte::BodyCoverPair::KOp (
    const state_type & state_q,
    state_type & out
) 
```




<hr>



### function MinvOp 

```C++
inline state_type tarte::BodyCoverPair::MinvOp (
    const state_type & state_p
) 
```




<hr>



### function ROp 

```C++
inline void tarte::BodyCoverPair::ROp (
    const state_type & state_p,
    state_type & out
) 
```




<hr>



### function ReadEffectiveOpenings 

```C++
inline half_state_type tarte::BodyCoverPair::ReadEffectiveOpenings () 
```




<hr>



### function get\_alpha\_contact\_stiffness 

```C++
inline ftype tarte::BodyCoverPair::get_alpha_contact_stiffness () 
```




<hr>



### function get\_c0 

```C++
inline ftype tarte::BodyCoverPair::get_c0 () 
```




<hr>



### function get\_contact\_stiffness 

```C++
inline ftype tarte::BodyCoverPair::get_contact_stiffness () 
```




<hr>



### function get\_epsilon\_smooth 

```C++
inline ftype tarte::BodyCoverPair::get_epsilon_smooth () 
```




<hr>



### function get\_eta\_stiffness 

```C++
inline ftype tarte::BodyCoverPair::get_eta_stiffness (
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function get\_gender 

```C++
inline Gender tarte::BodyCoverPair::get_gender () 
```




<hr>



### function get\_kt 

```C++
inline ftype tarte::BodyCoverPair::get_kt () 
```




<hr>



### function get\_lengths 

```C++
inline half_state_type tarte::BodyCoverPair::get_lengths (
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function get\_masses 

```C++
inline half_state_type tarte::BodyCoverPair::get_masses (
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function get\_noise\_ratio 

```C++
inline ftype tarte::BodyCoverPair::get_noise_ratio () 
```




<hr>



### function get\_rest\_positions 

```C++
inline half_state_type tarte::BodyCoverPair::get_rest_positions (
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function get\_rho0 

```C++
inline ftype tarte::BodyCoverPair::get_rho0 () 
```




<hr>



### function get\_stiffnesses 

```C++
inline Eigen::Vector< ftype, 4 > tarte::BodyCoverPair::get_stiffnesses (
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function get\_thicknesses 

```C++
inline half_state_type tarte::BodyCoverPair::get_thicknesses (
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function set\_alpha\_contact\_stiffness 

```C++
inline void tarte::BodyCoverPair::set_alpha_contact_stiffness (
    const ftype & alpha_contact_stiffness
) 
```




<hr>



### function set\_c0 

```C++
inline void tarte::BodyCoverPair::set_c0 (
    const ftype & c0
) 
```




<hr>



### function set\_contact\_stiffness 

```C++
inline void tarte::BodyCoverPair::set_contact_stiffness (
    const ftype & contact_stiffness
) 
```




<hr>



### function set\_epsilon\_smooth 

```C++
inline void tarte::BodyCoverPair::set_epsilon_smooth (
    const ftype & epsilon_smooth
) 
```




<hr>



### function set\_eta\_stiffness 

```C++
inline void tarte::BodyCoverPair::set_eta_stiffness (
    const ftype & eta_stiffness,
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function set\_gender 

```C++
inline void tarte::BodyCoverPair::set_gender (
    Gender gender
) 
```




<hr>



### function set\_kt 

```C++
inline void tarte::BodyCoverPair::set_kt (
    const ftype & kt
) 
```




<hr>



### function set\_length 

```C++
inline void tarte::BodyCoverPair::set_length (
    const ftype & length,
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function set\_masses [1/2]

```C++
inline void tarte::BodyCoverPair::set_masses (
    const half_state_type & masses,
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function set\_masses [2/2]

```C++
inline void tarte::BodyCoverPair::set_masses (
    const ftype & lower,
    const ftype & upper,
    const ftype & body,
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function set\_muscles\_activation 

```C++
inline void tarte::BodyCoverPair::set_muscles_activation (
    const ftype & ct_activity,
    const ftype & ta_activity,
    const ftype & lc_activity,
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function set\_noise\_ratio 

```C++
inline void tarte::BodyCoverPair::set_noise_ratio (
    const ftype & noise_ratio
) 
```




<hr>



### function set\_rest\_positions [1/2]

_--------------------------------_  _Control functions --------------------------------------_ _///_
```C++
inline void tarte::BodyCoverPair::set_rest_positions (
    const half_state_type & rest_positions,
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function set\_rest\_positions [2/2]

```C++
inline void tarte::BodyCoverPair::set_rest_positions (
    const ftype & lower,
    const ftype & upper,
    const ftype & body,
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function set\_rho0 

```C++
inline void tarte::BodyCoverPair::set_rho0 (
    const ftype & rho0
) 
```




<hr>



### function set\_stiffnesses [1/2]

```C++
inline void tarte::BodyCoverPair::set_stiffnesses (
    const Eigen::Vector< ftype, 4 > & stiffnesses,
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function set\_stiffnesses [2/2]

```C++
inline void tarte::BodyCoverPair::set_stiffnesses (
    const ftype & lower,
    const ftype & upper,
    const ftype & body,
    const ftype & coupling,
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function set\_thicknesses [1/2]

```C++
inline void tarte::BodyCoverPair::set_thicknesses (
    const half_state_type & thicknesses,
    FoldIdentifier fold_id=kBoth
) 
```




<hr>



### function set\_thicknesses [2/2]

```C++
inline void tarte::BodyCoverPair::set_thicknesses (
    const ftype & lower,
    const ftype & upper,
    const ftype & body,
    FoldIdentifier fold_id=kBoth
) 
```




<hr>
## Public Static Functions Documentation




### function get\_N 

```C++
static inline int tarte::BodyCoverPair::get_N () 
```




<hr>



### function get\_half\_N 

```C++
static inline int tarte::BodyCoverPair::get_half_N () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `include/vocal_folds/body_cover_pair.h`

