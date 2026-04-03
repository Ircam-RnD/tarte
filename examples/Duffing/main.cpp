#include <duffing.h>

#include <iostream>

int main(int, char*[])
{
    float samplerate = 44100;
    tarte::Duffing<float> duffingOscillator(samplerate);
    return 0;
}