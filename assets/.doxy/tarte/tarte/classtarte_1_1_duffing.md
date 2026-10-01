

# Class tarte::Duffing

**template &lt;typename T&gt;**



[**ClassList**](annotated.md) **>** [**tarte**](namespacetarte.md) **>** [**Duffing**](classtarte_1_1_duffing.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Duffing**](#function-duffing) (float sample\_rate) <br> |
|  void | [**Process**](#function-process) (T input\_force) <br> |
|  T | [**ReadDisplacement**](#function-readdisplacement) () <br> |
|  T | [**ReadVelocity**](#function-readvelocity) () <br> |
|  void | [**ReinitDsp**](#function-reinitdsp) (float sample\_rate) <br> |
|  void | [**set\_amplitude**](#function-set_amplitude) (T amplitude) <br> |
|  void | [**set\_decay\_time**](#function-set_decay_time) (T decay\_time) <br> |
|  void | [**set\_frequency**](#function-set_frequency) (T frequency) <br> |
|  void | [**set\_lambda\_ctrl**](#function-set_lambda_ctrl) (T lambda\_ctrl) <br> |
|  void | [**set\_linear\_parameters**](#function-set_linear_parameters) (T amplitude, T pulsation, T decay\_time) <br> |
|  void | [**set\_nonlinearity**](#function-set_nonlinearity) (T eta) <br> |
|  void | [**set\_physical\_parameters**](#function-set_physical_parameters) (T mass, T stiffness, T dissipation, T eta\_nl=0) <br> |




























## Public Functions Documentation




### function Duffing 

```C++
tarte::Duffing::Duffing (
    float sample_rate
) 
```




<hr>



### function Process 

```C++
void tarte::Duffing::Process (
    T input_force
) 
```




<hr>



### function ReadDisplacement 

```C++
inline T tarte::Duffing::ReadDisplacement () 
```




<hr>



### function ReadVelocity 

```C++
inline T tarte::Duffing::ReadVelocity () 
```




<hr>



### function ReinitDsp 

```C++
void tarte::Duffing::ReinitDsp (
    float sample_rate
) 
```




<hr>



### function set\_amplitude 

```C++
void tarte::Duffing::set_amplitude (
    T amplitude
) 
```




<hr>



### function set\_decay\_time 

```C++
void tarte::Duffing::set_decay_time (
    T decay_time
) 
```




<hr>



### function set\_frequency 

```C++
void tarte::Duffing::set_frequency (
    T frequency
) 
```




<hr>



### function set\_lambda\_ctrl 

```C++
inline void tarte::Duffing::set_lambda_ctrl (
    T lambda_ctrl
) 
```




<hr>



### function set\_linear\_parameters 

```C++
void tarte::Duffing::set_linear_parameters (
    T amplitude,
    T pulsation,
    T decay_time
) 
```




<hr>



### function set\_nonlinearity 

```C++
inline void tarte::Duffing::set_nonlinearity (
    T eta
) 
```




<hr>



### function set\_physical\_parameters 

```C++
void tarte::Duffing::set_physical_parameters (
    T mass,
    T stiffness,
    T dissipation,
    T eta_nl=0
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `include/duffing.h`

