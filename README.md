# Audio Compression Using MDCT in Scilab

A clean, hands-on simulation of an audio compression pipeline built from scratch in **Scilab**. This project demonstrates how modern audio codecs (like MP3, AAC, and Ogg Vorbis) compress sound using the **Modified Discrete Cosine Transform (MDCT)** and psychoacoustic quantization.

---

## 🎵 What is MDCT and Why Use It?
The **Modified Discrete Cosine Transform (MDCT)** is a specialized mathematical transform widely used in audio coding. It uses **50% overlapping frames** (lapped transforms), which provides two huge benefits:
* **No Blocking Artifacts:** Standard block transforms (like regular DCT or FFT) create harsh clicking or buzzing sounds at frame boundaries. MDCT eliminates this smoothly.
* **Efficient Compression:** It converts time-domain audio samples into frequency coefficients while maintaining critical sampling.

---

## 🚀 How the Code / Pipeline Works
1. **Framing & Windowing:** Splits the audio into overlapping blocks and applies a sine window.
2. **Forward MDCT:** Converts the audio frames into frequency coefficients.
3. **Perceptual Quantization:** Applies coarser quantization steps to high-frequency bands (since human ears are less sensitive to high-frequency details), driving many small coefficients down to zero.
4. **Entropy Estimation:** Measures the theoretical bits needed to store the remaining values to calculate the compression ratio and bitrate.
5. **Reconstruction (IMDCT):** Dequantizes the coefficients, runs the Inverse MDCT, and uses overlap-add synthesis to rebuild the audio signal.

---

## 📂 Repository Structure
* `mdct_audio_compression.sce` — The main Scilab simulation script.
* `my_audio.wav` — The input test audio file.
* `Reconstructed Audio/` — The output folder where the final compressed and reconstructed `.wav` file is automatically saved.

---

## ▶️ How to Run
1. Open **Scilab**.
2. Open the `mdct_audio_compression.sce` script in the SciNotes editor.
3. Run the script. Check the Scilab console for the compression ratio, SNR metrics, and listen to the reconstructed audio output!
