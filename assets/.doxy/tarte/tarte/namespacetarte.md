

# Namespace tarte



[**Namespace List**](namespaces.md) **>** [**tarte**](namespacetarte.md)


















## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**vowels**](namespacetarte_1_1vowels.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**Articulation**](classtarte_1_1_articulation.md) <br> |
| class | [**Biquad**](classtarte_1_1_biquad.md) <br> |
| class | [**BodyCoverPair**](classtarte_1_1_body_cover_pair.md) &lt;typename ftype&gt;<br> |
| class | [**BodyCoverVF**](classtarte_1_1_body_cover_v_f.md) &lt;typename T&gt;<br> |
| class | [**Bow**](classtarte_1_1_bow.md) &lt;class T&gt;<br> |
| class | [**Duffing**](classtarte_1_1_duffing.md) &lt;typename T&gt;<br> |
| class | [**NoiseGenerator**](classtarte_1_1_noise_generator.md) <br> |
| class | [**SingleReed**](classtarte_1_1_single_reed.md) &lt;typename ftype&gt;<br> |
| class | [**VFPairExample**](classtarte_1_1_v_f_pair_example.md) &lt;typename ftype&gt;<br> |
| class | [**Voice**](classtarte_1_1_voice.md) &lt;vf\_pair, typename ftype&gt;<br> |
| class | [**WebsterFDTD**](classtarte_1_1_webster_f_d_t_d.md) &lt;typename ftype, kMaxN&gt;<br> |


## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**BiquadMode**](#enum-biquadmode)  <br> |
| enum  | [**BowMode**](#enum-bowmode)  <br> |
| enum  | [**FoldIdentifier**](#enum-foldidentifier)  <br> |
| enum  | [**Gender**](#enum-gender)  <br> |
| enum  | [**NoiseColor**](#enum-noisecolor)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  auto | [**ClipEigen**](#function-clipeigen) (Eigen::Ref&lt; const Eigen::ArrayX&lt; Derived &gt; &gt; array, const ftype & min, const ftype & max) <br> |
|  auto | [**ClipEigen**](#function-clipeigen) (const Eigen::ArrayBase&lt; Derived &gt; & array, const ftype & min, const ftype & max) <br> |
|  auto | [**SafeSetEigen**](#function-safeseteigen) (Eigen::ArrayBase&lt; Derived1 &gt; & array, Eigen::Ref&lt; const Eigen::ArrayX&lt; Derived2 &gt; &gt; array2) <br> |
|  auto | [**SafeSetEigen**](#function-safeseteigen) (Eigen::ArrayBase&lt; Derived1 &gt; & array, const Eigen::ArrayBase&lt; Derived2 &gt; & array2) <br> |
|  void | [**WriteWav**](#function-writewav) (const std::string & path, const std::vector&lt; float &gt; & samples, int sampleRate) <br> |
|  void | [**normalize\_vector**](#function-normalize_vector) (std::vector&lt; T &gt; & vec) <br> |
|  int | [**sgn**](#function-sgn) (T val) <br> |
|  auto | [**smoothRampMatrix**](#function-smoothrampmatrix) (const Eigen::MatrixBase&lt; derived &gt; & x, ftype epsilon=1) <br> |
|  ftype | [**softplus**](#function-softplus) (const ftype & x, ftype epsilon=1, ftype threshold=40) <br> |
|  ftype | [**softplusDerivative**](#function-softplusderivative) (const ftype & x, ftype epsilon=1, ftype threshold=40) <br> |
|  auto | [**softplusDerivativeMatrix**](#function-softplusderivativematrix) (const Eigen::MatrixBase&lt; derived &gt; & x, ftype epsilon=1, ftype threshold=40) <br> |
|  auto | [**softplusMatrix**](#function-softplusmatrix) (const Eigen::MatrixBase&lt; derived &gt; & x, ftype epsilon=1, ftype threshold=40) <br> |




























## Public Types Documentation




### enum BiquadMode 

```C++
enum tarte::BiquadMode {
    kLowPass,
    kHighPass,
    kBandPass,
    kNotch,
    kPeak,
    kLowShelf,
    kHighShelf
};
```




<hr>



### enum BowMode 

```C++
enum tarte::BowMode {
    kMatusiak,
    kVigue,
    kTerrien
};
```




<hr>



### enum FoldIdentifier 

```C++
enum tarte::FoldIdentifier {
    kRight,
    kLeft,
    kBoth
};
```




<hr>



### enum Gender 

```C++
enum tarte::Gender {
    MALE,
    FEMALE
};
```




<hr>



### enum NoiseColor 

```C++
enum tarte::NoiseColor {
    White,
    Pink
};
```




<hr>
## Public Functions Documentation




### function ClipEigen 

```C++
template<typename Derived, typename ftype>
auto tarte::ClipEigen (
    Eigen::Ref< const Eigen::ArrayX< Derived > > array,
    const ftype & min,
    const ftype & max
) 
```




<hr>



### function ClipEigen 

```C++
template<typename Derived, typename ftype>
auto tarte::ClipEigen (
    const Eigen::ArrayBase< Derived > & array,
    const ftype & min,
    const ftype & max
) 
```




<hr>



### function SafeSetEigen 

```C++
template<typename Derived1, typename Derived2>
auto tarte::SafeSetEigen (
    Eigen::ArrayBase< Derived1 > & array,
    Eigen::Ref< const Eigen::ArrayX< Derived2 > > array2
) 
```




<hr>



### function SafeSetEigen 

```C++
template<typename Derived1, typename Derived2>
auto tarte::SafeSetEigen (
    Eigen::ArrayBase< Derived1 > & array,
    const Eigen::ArrayBase< Derived2 > & array2
) 
```




<hr>



### function WriteWav 

```C++
void tarte::WriteWav (
    const std::string & path,
    const std::vector< float > & samples,
    int sampleRate
) 
```




<hr>



### function normalize\_vector 

```C++
template<typename T>
void tarte::normalize_vector (
    std::vector< T > & vec
) 
```




<hr>



### function sgn 

```C++
template<typename T>
int tarte::sgn (
    T val
) 
```




<hr>



### function smoothRampMatrix 

```C++
template<typename derived, typename ftype>
auto tarte::smoothRampMatrix (
    const Eigen::MatrixBase< derived > & x,
    ftype epsilon=1
) 
```




<hr>



### function softplus 

```C++
template<typename ftype>
ftype tarte::softplus (
    const ftype & x,
    ftype epsilon=1,
    ftype threshold=40
) 
```




<hr>



### function softplusDerivative 

```C++
template<typename ftype>
ftype tarte::softplusDerivative (
    const ftype & x,
    ftype epsilon=1,
    ftype threshold=40
) 
```




<hr>



### function softplusDerivativeMatrix 

```C++
template<typename derived, typename ftype>
auto tarte::softplusDerivativeMatrix (
    const Eigen::MatrixBase< derived > & x,
    ftype epsilon=1,
    ftype threshold=40
) 
```




<hr>



### function softplusMatrix 

```C++
template<typename derived, typename ftype>
auto tarte::softplusMatrix (
    const Eigen::MatrixBase< derived > & x,
    ftype epsilon=1,
    ftype threshold=40
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/duffing.cpp`

