===============Variables================
ImageSize: Size of the binary image transmitted
NoNoise_RxSignal: Received signal without noise
LowNoise_RxSignal: Received signal with noise (variance = -25 dB)
HighNoise_RxSignal: Received signal with noise (variance = -10 dB)


===============Setting=======================
To generate data, we use binary images as the transmitted data. Every pixel of the binary image is represented by a 0 (black) or a 1 (white). Quadrature Phase Shift Keying(QPSK) mapping is used for the data bits, and the QPSK symbol stream is padded with zeros to make its size a multiple of the ODFM symbol length. We use common pilot symbol for all the OFDM pilot indices.  The cyclic prefixed sequence is convolved with the channel vector to obtain the noiseless received signal. For simplicity, we have avoided error correction coding, passband transmission, pulse shaping, etc.  

===============Hints======================
For each received signal, do the following:
 1. Remove cyclic prefixing (length of cyclic prefixing = channel length -1)
 2. Extract the pilot part of the received signal in the frequency domain
 3. Estimate the channel using a filtering/smoothing technique
 4. Extract the data symbols and remove the zero padding
 5. Use QAM demapper to reconstruct transmitted data bits
 6. Display the reconstructed image using the following code where '#' should be replaced with the recovered data bits: 
	DecodeImage = cast(255*reshape(#,ImageSize),'uint8');
        imshow(DecodeImage);
