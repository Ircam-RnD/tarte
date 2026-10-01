

# Class tarte::Articulation



[**ClassList**](annotated.md) **>** [**tarte**](namespacetarte.md) **>** [**Articulation**](classtarte_1_1_articulation.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Articulation**](#function-articulation) ([**vowels::Vowel**](structtarte_1_1vowels_1_1_vowel.md) const & v=vowels::a) <br> |
|  void | [**FindClosestVowels**](#function-findclosestvowels) (float targetF0, float targetF1, size\_t N, std::vector&lt; [**vowels::Vowel**](structtarte_1_1vowels_1_1_vowel.md) &gt; & closest\_vowels, std::vector&lt; float &gt; & alphas) <br> |
|  void | [**Interpolate2Vowels**](#function-interpolate2vowels) ([**vowels::Vowel**](structtarte_1_1vowels_1_1_vowel.md) const & v1, [**vowels::Vowel**](structtarte_1_1vowels_1_1_vowel.md) const & v2, float alpha) <br> |
|  void | [**InterpolateNVowels**](#function-interpolatenvowels) ([**vowels::Vowel**](structtarte_1_1vowels_1_1_vowel.md) const \* vs, float const \* alphas, size\_t const n) <br> |
|  void | [**SetFromFormants**](#function-setfromformants) (float targetF0, float targetF1, size\_t N=3) <br> |
|  void | [**SetFromVowel**](#function-setfromvowel) ([**vowels::Vowel**](structtarte_1_1vowels_1_1_vowel.md) const & v) <br> |
|  void | [**getAreas**](#function-getareas) (T const \* evaluation\_positions, T \* out, std::size\_t const size) <br> |




























## Public Functions Documentation




### function Articulation 

```C++
inline tarte::Articulation::Articulation (
    vowels::Vowel const & v=vowels::a
) 
```




<hr>



### function FindClosestVowels 

```C++
inline void tarte::Articulation::FindClosestVowels (
    float targetF0,
    float targetF1,
    size_t N,
    std::vector< vowels::Vowel > & closest_vowels,
    std::vector< float > & alphas
) 
```




<hr>



### function Interpolate2Vowels 

```C++
inline void tarte::Articulation::Interpolate2Vowels (
    vowels::Vowel const & v1,
    vowels::Vowel const & v2,
    float alpha
) 
```




<hr>



### function InterpolateNVowels 

```C++
inline void tarte::Articulation::InterpolateNVowels (
    vowels::Vowel const * vs,
    float const * alphas,
    size_t const n
) 
```




<hr>



### function SetFromFormants 

```C++
inline void tarte::Articulation::SetFromFormants (
    float targetF0,
    float targetF1,
    size_t N=3
) 
```




<hr>



### function SetFromVowel 

```C++
inline void tarte::Articulation::SetFromVowel (
    vowels::Vowel const & v
) 
```




<hr>



### function getAreas 

```C++
template<typename T>
inline void tarte::Articulation::getAreas (
    T const * evaluation_positions,
    T * out,
    std::size_t const size
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `include/articulation.h`

