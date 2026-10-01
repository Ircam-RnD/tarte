

# Class tarte::Biquad



[**ClassList**](annotated.md) **>** [**tarte**](namespacetarte.md) **>** [**Biquad**](classtarte_1_1_biquad.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Biquad**](#function-biquad) (float sample\_rate=44100, BiquadMode mode=kLowPass, float frequency=100, float dB\_gain=0, float Q=0.7) <br> |
|  void | [**InitializeState**](#function-initializestate) (double initial\_value) <br> |
|  double | [**Process**](#function-process) (double input) <br> |
|  void | [**set\_Q**](#function-set_q) (float Q\_) <br> |
|  void | [**set\_freq**](#function-set_freq) (float freq\_) <br> |
|  void | [**set\_gain**](#function-set_gain) (float dB\_gain\_) <br> |
|  void | [**set\_mode**](#function-set_mode) (BiquadMode mode) <br> |




























## Public Functions Documentation




### function Biquad 

```C++
inline tarte::Biquad::Biquad (
    float sample_rate=44100,
    BiquadMode mode=kLowPass,
    float frequency=100,
    float dB_gain=0,
    float Q=0.7
) 
```




<hr>



### function InitializeState 

```C++
inline void tarte::Biquad::InitializeState (
    double initial_value
) 
```




<hr>



### function Process 

```C++
inline double tarte::Biquad::Process (
    double input
) 
```




<hr>



### function set\_Q 

```C++
inline void tarte::Biquad::set_Q (
    float Q_
) 
```




<hr>



### function set\_freq 

```C++
inline void tarte::Biquad::set_freq (
    float freq_
) 
```




<hr>



### function set\_gain 

```C++
inline void tarte::Biquad::set_gain (
    float dB_gain_
) 
```




<hr>



### function set\_mode 

```C++
inline void tarte::Biquad::set_mode (
    BiquadMode mode
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `include/utility/biquad.h`

