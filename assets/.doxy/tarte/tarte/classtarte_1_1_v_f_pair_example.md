

# Class tarte::VFPairExample

**template &lt;typename ftype&gt;**



[**ClassList**](annotated.md) **>** [**tarte**](namespacetarte.md) **>** [**VFPairExample**](classtarte_1_1_v_f_pair_example.md)






















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Eigen::Vector&lt; ftype, half\_N &gt; | [**half\_state\_type**](#typedef-half_state_type)  <br> |
| typedef ftype | [**scalar\_type**](#typedef-scalar_type)  <br> |
| typedef Eigen::Vector&lt; ftype, N &gt; | [**state\_type**](#typedef-state_type)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  ftype | [**ComputeFlow**](#function-computeflow) (const ftype & Psub, const ftype & Psup) <br> |
|  void | [**EffectiveAreas**](#function-effectiveareas) (state\_type & out\_area\_P\_sub, state\_type & out\_area\_P\_sup) <br> |
|  ftype | [**Enl**](#function-enl) (const state\_type & state\_q, bool recompute\_intermediary=false) <br> |
|  void | [**FillIntermediary**](#function-fillintermediary) (const state\_type & state\_q) <br> |
|  void | [**Fnl**](#function-fnl) (const state\_type & state\_q, state\_type & out, bool recompute\_intermediary=false) <br> |
|  void | [**KOp**](#function-kop) (const state\_type & state\_q, state\_type & out) <br> |
|  state\_type | [**MinvOp**](#function-minvop) (const state\_type & state\_p) <br> |
|  void | [**ROp**](#function-rop) (const state\_type & state\_p, state\_type & out) <br> |
|  half\_state\_type | [**ReadEffectiveOpenings**](#function-readeffectiveopenings) () <br> |
|   | [**VFPairExample**](#function-vfpairexample) () <br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  int | [**get\_N**](#function-get_n) () <br> |
|  int | [**get\_half\_N**](#function-get_half_n) () <br> |


























## Public Types Documentation




### typedef half\_state\_type 

```C++
using tarte::VFPairExample< ftype >::half_state_type = Eigen::Vector<ftype, half_N>;
```




<hr>



### typedef scalar\_type 

```C++
using tarte::VFPairExample< ftype >::scalar_type = ftype;
```




<hr>



### typedef state\_type 

```C++
using tarte::VFPairExample< ftype >::state_type = Eigen::Vector<ftype, N>;
```




<hr>
## Public Functions Documentation




### function ComputeFlow 

```C++
inline ftype tarte::VFPairExample::ComputeFlow (
    const ftype & Psub,
    const ftype & Psup
) 
```




<hr>



### function EffectiveAreas 

```C++
inline void tarte::VFPairExample::EffectiveAreas (
    state_type & out_area_P_sub,
    state_type & out_area_P_sup
) 
```




<hr>



### function Enl 

```C++
inline ftype tarte::VFPairExample::Enl (
    const state_type & state_q,
    bool recompute_intermediary=false
) 
```




<hr>



### function FillIntermediary 

```C++
inline void tarte::VFPairExample::FillIntermediary (
    const state_type & state_q
) 
```




<hr>



### function Fnl 

```C++
inline void tarte::VFPairExample::Fnl (
    const state_type & state_q,
    state_type & out,
    bool recompute_intermediary=false
) 
```




<hr>



### function KOp 

```C++
inline void tarte::VFPairExample::KOp (
    const state_type & state_q,
    state_type & out
) 
```




<hr>



### function MinvOp 

```C++
inline state_type tarte::VFPairExample::MinvOp (
    const state_type & state_p
) 
```




<hr>



### function ROp 

```C++
inline void tarte::VFPairExample::ROp (
    const state_type & state_p,
    state_type & out
) 
```




<hr>



### function ReadEffectiveOpenings 

```C++
inline half_state_type tarte::VFPairExample::ReadEffectiveOpenings () 
```




<hr>



### function VFPairExample 

```C++
inline tarte::VFPairExample::VFPairExample () 
```




<hr>
## Public Static Functions Documentation




### function get\_N 

```C++
static inline int tarte::VFPairExample::get_N () 
```




<hr>



### function get\_half\_N 

```C++
static inline int tarte::VFPairExample::get_half_N () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `include/vocal_folds/vf_pair_model_example.h`

