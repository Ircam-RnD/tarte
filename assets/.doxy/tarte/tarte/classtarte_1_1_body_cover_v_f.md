

# Class tarte::BodyCoverVF

**template &lt;typename T&gt;**



[**ClassList**](annotated.md) **>** [**tarte**](namespacetarte.md) **>** [**BodyCoverVF**](classtarte_1_1_body_cover_v_f.md)




























## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  const Eigen::Matrix&lt; T, 4, 3 &gt; | [**elongation\_matrix\_**](#variable-elongation_matrix_)   = `/* multi line expression */`<br> |














## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BodyCoverVF**](#function-bodycovervf) () <br> |
|  void | [**ComputeParametersFromMuscularActivity**](#function-computeparametersfrommuscularactivity) () <br> |
|  void | [**FillDissipationCoefficients**](#function-filldissipationcoefficients) () <br> |
|  void | [**SetCricoarytenoidActivity**](#function-setcricoarytenoidactivity) (const T & activity) <br> |
|  void | [**SetCricothyroidActivity**](#function-setcricothyroidactivity) (const T & activity) <br> |
|  void | [**SetThyroarytenoidActivity**](#function-setthyroarytenoidactivity) (const T & activity) <br> |
|  Eigen::DiagonalMatrix&lt; T, 4 &gt; | [**dissipation\_coefficients**](#function-dissipation_coefficients) () <br> |
|  T | [**eta\_stiffness**](#function-eta_stiffness) () <br> |
|  Gender | [**get\_gender**](#function-get_gender) () <br> |
|  Eigen::Vector&lt; T, 3 &gt; | [**lengths**](#function-lengths) () <br> |
|  Eigen::DiagonalMatrix&lt; T, 3 &gt; | [**masses**](#function-masses) () <br> |
|  Eigen::Vector&lt; T, 3 &gt; | [**rest\_positions**](#function-rest_positions) () <br> |
|  void | [**set\_eta\_stiffness**](#function-set_eta_stiffness) (const T & eta\_stiffness) <br> |
|  void | [**set\_gender**](#function-set_gender) (Gender gender) <br> |
|  void | [**set\_length**](#function-set_length) (const T & length) <br> |
|  void | [**set\_masses**](#function-set_masses) (const Eigen::Vector&lt; T, 3 &gt; & masses) <br> |
|  void | [**set\_rest\_positions**](#function-set_rest_positions) (const Eigen::Vector&lt; T, 3 &gt; & rest\_positions) <br> |
|  void | [**set\_stiffnesses**](#function-set_stiffnesses) (const Eigen::Vector&lt; T, 4 &gt; & stiffnesses) <br> |
|  void | [**set\_thicknesses**](#function-set_thicknesses) (const Eigen::Vector&lt; T, 3 &gt; & thicknesses) <br> |
|  Eigen::DiagonalMatrix&lt; T, 4 &gt; | [**stiffnesses**](#function-stiffnesses) () <br> |
|  Eigen::Vector&lt; T, 3 &gt; | [**thicknesses**](#function-thicknesses) () <br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  int | [**getN**](#function-getn) () <br> |


























## Public Static Attributes Documentation




### variable elongation\_matrix\_ 

```C++
const Eigen::Matrix<T, 4, 3> tarte::BodyCoverVF< T >::elongation_matrix_;
```




<hr>
## Public Functions Documentation




### function BodyCoverVF 

```C++
inline tarte::BodyCoverVF::BodyCoverVF () 
```




<hr>



### function ComputeParametersFromMuscularActivity 

```C++
inline void tarte::BodyCoverVF::ComputeParametersFromMuscularActivity () 
```




<hr>



### function FillDissipationCoefficients 

```C++
inline void tarte::BodyCoverVF::FillDissipationCoefficients () 
```




<hr>



### function SetCricoarytenoidActivity 

```C++
inline void tarte::BodyCoverVF::SetCricoarytenoidActivity (
    const T & activity
) 
```




<hr>



### function SetCricothyroidActivity 

```C++
inline void tarte::BodyCoverVF::SetCricothyroidActivity (
    const T & activity
) 
```




<hr>



### function SetThyroarytenoidActivity 

```C++
inline void tarte::BodyCoverVF::SetThyroarytenoidActivity (
    const T & activity
) 
```




<hr>



### function dissipation\_coefficients 

```C++
inline Eigen::DiagonalMatrix< T, 4 > tarte::BodyCoverVF::dissipation_coefficients () 
```




<hr>



### function eta\_stiffness 

```C++
inline T tarte::BodyCoverVF::eta_stiffness () 
```




<hr>



### function get\_gender 

```C++
inline Gender tarte::BodyCoverVF::get_gender () 
```




<hr>



### function lengths 

```C++
inline Eigen::Vector< T, 3 > tarte::BodyCoverVF::lengths () 
```




<hr>



### function masses 

```C++
inline Eigen::DiagonalMatrix< T, 3 > tarte::BodyCoverVF::masses () 
```




<hr>



### function rest\_positions 

```C++
inline Eigen::Vector< T, 3 > tarte::BodyCoverVF::rest_positions () 
```




<hr>



### function set\_eta\_stiffness 

```C++
inline void tarte::BodyCoverVF::set_eta_stiffness (
    const T & eta_stiffness
) 
```




<hr>



### function set\_gender 

```C++
inline void tarte::BodyCoverVF::set_gender (
    Gender gender
) 
```




<hr>



### function set\_length 

```C++
inline void tarte::BodyCoverVF::set_length (
    const T & length
) 
```




<hr>



### function set\_masses 

```C++
inline void tarte::BodyCoverVF::set_masses (
    const Eigen::Vector< T, 3 > & masses
) 
```




<hr>



### function set\_rest\_positions 

```C++
inline void tarte::BodyCoverVF::set_rest_positions (
    const Eigen::Vector< T, 3 > & rest_positions
) 
```




<hr>



### function set\_stiffnesses 

```C++
inline void tarte::BodyCoverVF::set_stiffnesses (
    const Eigen::Vector< T, 4 > & stiffnesses
) 
```




<hr>



### function set\_thicknesses 

```C++
inline void tarte::BodyCoverVF::set_thicknesses (
    const Eigen::Vector< T, 3 > & thicknesses
) 
```




<hr>



### function stiffnesses 

```C++
inline Eigen::DiagonalMatrix< T, 4 > tarte::BodyCoverVF::stiffnesses () 
```




<hr>



### function thicknesses 

```C++
inline Eigen::Vector< T, 3 > tarte::BodyCoverVF::thicknesses () 
```




<hr>
## Public Static Functions Documentation




### function getN 

```C++
static inline int tarte::BodyCoverVF::getN () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `include/vocal_folds/body_cover_vf.h`

