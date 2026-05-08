The RIFF container, utilizing the ATRAC3 encoder, supports bitrates of 66, 105, and 132 kbps.

This work is derived from the project at https://github.com/hilman2/atrac3-rs with modifications. Special thanks to hilman2.

usage：

at3enc -e input.wav(default 132kb)

at3enc -e input.wav -b 66

at3enc -e input.wav -b 105


The official ATRAC3 encoder only supports 44.1 kHz WAV input, which is very unfriendly. To address this, I have added a resampling feature: if the input WAV is not at 44.1 kHz, it will be automatically resampled before encoding.
