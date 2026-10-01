

# File audiowrite.h

[**File List**](files.md) **>** [**include**](dir_d44c64559bbebec7f509842c48db8b23.md) **>** [**utility**](dir_3a4d35156996fb70540e167b84a39bd1.md) **>** [**audiowrite.h**](audiowrite_8h.md)

[Go to the documentation of this file](audiowrite_8h.md)


```C++
#pragma once
#include "external/miniaudio.h"
#include <stdexcept>
#include <string>
#include <vector>

namespace tarte {

void WriteWav(const std::string& path, const std::vector<float>& samples, int sampleRate)
{
    ma_encoder_config config = ma_encoder_config_init(ma_encoding_format_wav, // format
                                                      ma_format_f32,          // sample type
                                                      1,                      // channels (mono)
                                                      sampleRate);

    ma_encoder encoder;
    if (ma_encoder_init_file(path.c_str(), &config, &encoder) != MA_SUCCESS)
        throw std::runtime_error("Failed to open file for writing: " + path);

    ma_uint64 framesWritten;
    ma_encoder_write_pcm_frames(&encoder, samples.data(), samples.size(), &framesWritten);

    ma_encoder_uninit(&encoder);
}

} // namespace tarte
```


